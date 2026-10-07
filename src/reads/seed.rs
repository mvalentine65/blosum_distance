//! Seed index for the ungapped adapter scan.
//!
//! Not from fastp — this replaces its brute-force offset loop with an
//! equivalent that tests far fewer offsets.
//!
//! Pass 1 accepts an offset when the adapter matches with at most
//! `cmplen / 8` mismatches. Cut the adapter into `allowed + 1` disjoint blocks:
//! if every block held a mismatch there would be `allowed + 1` of them, one too
//! many. So **at least one block must match exactly**, and only offsets where
//! some block lands exactly can possibly match. Finding those costs one rolling
//! k-mer pass over the read instead of a comparison at all 146 offsets.
//!
//! The filter never discards a match fastp would have found — it is the
//! pigeonhole principle, not a heuristic.

/// Block length. 6 bases pack into 12 bits, so the lookup table is 4096 entries.
pub const SEED_K: usize = 6;
/// Bytes a vector compare may read past the data it is given.
pub const PAD: usize = 32;
/// Read bytes one step of the packed scan covers.
#[cfg(target_arch = "x86_64")]
const SPAN: usize = 32 + SEED_K - 1;
const TABLE_SIZE: usize = 1 << (2 * SEED_K);

/// 2-bit code per base, 0xFF for anything that cannot start a seed. A table
/// keeps the rolling scan branchless — it runs once per base of every read.
const BASE_CODE: [u8; 256] = {
    let mut t = [0xFFu8; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t
};

#[inline]
fn base_code(b: u8) -> Option<u8> {
    let c = BASE_CODE[b as usize];
    if c == 0xFF {
        None
    } else {
        Some(c)
    }
}

/// An adapter with its seed table built once, ahead of the run.
#[derive(Debug, Clone)]
pub struct PreparedAdapter {
    pub seq: Vec<u8>,
    /// `seq` with PAD zero bytes after it, for compares that read past its end.
    pub padded: Vec<u8>,
    /// For each k-mer code, a bitmask of the blocks holding it.
    table: Vec<u16>,
    /// Number of blocks; 0 means the adapter is too short to filter and the
    /// caller must scan every offset.
    nblocks: usize,
    /// Each block in the table as its two packed halves and its offset in the adapter.
    halves: Vec<(u8, u8, usize)>,
}

/// A base in 2 bits: A 0, C 1, T 2, G 3. Any other byte lands on one of them.
#[inline]
const fn code2(b: u8) -> u8 {
    (b >> 1) & 3
}

/// Three bases packed into 6 bits.
fn pack3(bases: &[u8]) -> u8 {
    code2(bases[0]) | code2(bases[1]) << 2 | code2(bases[2]) << 4
}

impl PreparedAdapter {
    pub fn new(seq: &[u8]) -> Self {
        let alen = seq.len();
        // Blocks needed to guarantee a clean one, and blocks available.
        let allowed = alen / 8;
        let needed = allowed + 1;
        let available = alen / SEED_K;

        if alen < SEED_K || available < needed || needed > 16 {
            // Cannot guarantee the pigeonhole (or cannot fit the mask), so the
            // caller falls back to the exhaustive scan.
            return PreparedAdapter { seq: seq.to_vec(), padded: padded(seq), table: Vec::new(), nblocks: 0, halves: Vec::new() };
        }

        let nblocks = needed;
        let mut table = vec![0u16; TABLE_SIZE];
        let mut halves = Vec::new();
        for b in 0..nblocks {
            let start = b * SEED_K;
            let block = &seq[start..start + SEED_K];
            if let Some(code) = encode(block) {
                table[code as usize] |= 1 << b;
                halves.push((pack3(&block[..3]), pack3(&block[3..]), start));
            }
        }
        PreparedAdapter { seq: seq.to_vec(), padded: padded(seq), table, nblocks, halves }
    }

    /// The rolling scan's marks, 32 read positions a step: each position's six
    /// bases are packed as two halves of three and compared with every block's.
    /// False, with nothing marked, when the read holds a byte other than A, C,
    /// G or T, whose packed code would pass for a base.
    ///
    /// # Safety
    /// The CPU must support AVX2 (`seqops::avx2()`); `read` holds SPAN bytes or more.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn mark_packed_avx2(&self, read: &[u8], max_pos: usize, seen: &mut [u64]) -> bool {
        use std::arch::x86_64::*;
        let (len, p) = (read.len(), read.as_ptr());
        debug_assert!(len >= SPAN);
        // SAFETY: every load below starts 32 bytes or more before the end of `read`.
        unsafe {
            let three = _mm256_set1_epi8(3);
            let load = |at: usize| _mm256_loadu_si256(p.add(at) as *const __m256i);
            let codes = |at: usize| _mm256_and_si256(_mm256_srli_epi16::<1>(load(at)), three);

            let letters = _mm256_broadcastsi128_si256(_mm_setr_epi8(
                b'A' as i8, b'C' as i8, b'T' as i8, b'G' as i8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ));
            let mut i = 0;
            loop {
                let at = i.min(len - 32);
                let x = load(at);
                let back = _mm256_shuffle_epi8(letters, _mm256_and_si256(_mm256_srli_epi16::<1>(x), three));
                if _mm256_movemask_epi8(_mm256_cmpeq_epi8(x, back)) != -1 {
                    return false;
                }
                if at == len - 32 {
                    break;
                }
                i += 32;
            }

            // an adapter has 16 blocks at most (PreparedAdapter::new)
            let mut blocks = [(_mm256_setzero_si256(), _mm256_setzero_si256(), 0usize); 16];
            for (slot, &(lo, hi, shift)) in blocks.iter_mut().zip(&self.halves) {
                *slot = (_mm256_set1_epi8(lo as i8), _mm256_set1_epi8(hi as i8), shift);
            }
            let blocks = &blocks[..self.halves.len()];
            let far = self.halves.last().map_or(0, |h| h.2);
            let last = (len - SEED_K).min(max_pos + far);
            let mut i = 0;
            loop {
                // the final step is pulled back to end on the read's last base
                let at = i.min(len - SPAN);
                let pack = |from: usize| {
                    _mm256_or_si256(
                        _mm256_or_si256(codes(from), _mm256_slli_epi16::<2>(codes(from + 1))),
                        _mm256_slli_epi16::<4>(codes(from + 2)),
                    )
                };
                let (lo, hi) = (pack(at), pack(at + 3));
                for &(block_lo, block_hi, shift) in blocks {
                    let hit = _mm256_and_si256(_mm256_cmpeq_epi8(lo, block_lo), _mm256_cmpeq_epi8(hi, block_hi));
                    let mut hits = _mm256_movemask_epi8(hit) as u32;
                    while hits != 0 {
                        let start = at + hits.trailing_zeros() as usize;
                        hits &= hits - 1;
                        if start >= shift && start - shift <= max_pos {
                            let pos = start - shift;
                            seen[pos / 64] |= 1 << (pos % 64);
                        }
                    }
                }
                if at == len - SPAN || at + 32 > last {
                    break;
                }
                i += 32;
            }
        }
        true
    }

    #[inline]
    pub fn can_filter(&self) -> bool {
        self.nblocks > 0
    }

    /// Offsets in `0..=max_pos` where some block lands exactly, ascending.
    ///
    /// Every offset fastp's pass 1 would accept appears here; the caller still
    /// verifies each one, so extra candidates only cost time, never accuracy.
    pub fn candidates(
        &self,
        read: &[u8],
        max_pos: usize,
        seen: &mut Vec<u64>,
        out: &mut Vec<usize>,
    ) {
        out.clear();
        if self.nblocks == 0 || read.len() < SEED_K {
            return;
        }
        // A bitset keeps the offsets sorted and deduplicated for free. It
        // lives in the caller's scratch: allocating it per read cost more than
        // the filter saved.
        let words = max_pos / 64 + 1;
        seen.clear();
        seen.resize(words, 0);

        #[cfg(target_arch = "x86_64")]
        // SAFETY: AVX2 is checked and the read is long enough for a step.
        if read.len() >= SPAN && crate::reads::seqops::avx2() && unsafe { self.mark_packed_avx2(read, max_pos, seen) } {
            push_marked(seen, out);
            return;
        }

        let mut code: u32 = 0;
        let mut valid = 0usize;
        const MASK: u32 = (1 << (2 * SEED_K)) - 1;
        for (i, &b) in read.iter().enumerate() {
            let c = BASE_CODE[b as usize];
            if c == 0xFF {
                valid = 0;
                continue;
            }
            code = ((code << 2) | u32::from(c)) & MASK;
            valid += 1;
            if valid < SEED_K {
                continue;
            }
            let start = i + 1 - SEED_K;
            let mut mask = self.table[code as usize];
            while mask != 0 {
                let b = mask.trailing_zeros() as usize;
                mask &= mask - 1;
                let shift = b * SEED_K;
                if start >= shift {
                    let pos = start - shift;
                    if pos <= max_pos {
                        seen[pos / 64] |= 1 << (pos % 64);
                    }
                }
            }
        }

        push_marked(seen, out);
    }
}

/// The marked offsets, ascending.
fn push_marked(seen: &[u64], out: &mut Vec<usize>) {
    for (w, &word) in seen.iter().enumerate() {
        let mut bits = word;
        while bits != 0 {
            let b = bits.trailing_zeros() as usize;
            bits &= bits - 1;
            out.push(w * 64 + b);
        }
    }
}

fn padded(seq: &[u8]) -> Vec<u8> {
    let mut out = seq.to_vec();
    out.resize(seq.len() + PAD, 0);
    out
}

fn encode(kmer: &[u8]) -> Option<u32> {
    let mut code = 0u32;
    for &b in kmer {
        code = (code << 2) | u32::from(base_code(b)?);
    }
    Some(code)
}

#[cfg(test)]
#[path = "../tests/seed.rs"]
mod tests;
