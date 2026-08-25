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
    /// For each k-mer code, a bitmask of the blocks holding it.
    table: Vec<u16>,
    /// Number of blocks; 0 means the adapter is too short to filter and the
    /// caller must scan every offset.
    nblocks: usize,
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
            return PreparedAdapter { seq: seq.to_vec(), table: Vec::new(), nblocks: 0 };
        }

        let nblocks = needed;
        let mut table = vec![0u16; TABLE_SIZE];
        for b in 0..nblocks {
            let start = b * SEED_K;
            if let Some(code) = encode(&seq[start..start + SEED_K]) {
                table[code as usize] |= 1 << b;
            }
        }
        PreparedAdapter { seq: seq.to_vec(), table, nblocks }
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

        for (w, &word) in seen.iter().enumerate() {
            let mut bits = word;
            while bits != 0 {
                let b = bits.trailing_zeros() as usize;
                bits &= bits - 1;
                out.push(w * 64 + b);
            }
        }
    }
}

fn encode(kmer: &[u8]) -> Option<u32> {
    let mut code = 0u32;
    for &b in kmer {
        code = (code << 2) | u32::from(base_code(b)?);
    }
    Some(code)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn truseq_is_filterable_and_finds_its_own_offset() {
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let prep = PreparedAdapter::new(adapter);
        assert!(prep.can_filter());

        let mut read = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let at = read.len();
        read.extend_from_slice(adapter);

        let mut out = Vec::new();
        prep.candidates(&read, read.len() - adapter.len(), &mut Vec::<u64>::new(), &mut out);
        assert!(out.contains(&at), "offset {at} missing from {out:?}");
    }

    /// The guarantee that matters: an adapter carrying the maximum tolerated
    /// number of mismatches is still proposed as a candidate.
    #[test]
    fn a_match_at_the_mismatch_limit_survives_the_filter() {
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let prep = PreparedAdapter::new(adapter);
        let allowed = adapter.len() / 8; // 4

        let mut damaged = adapter.to_vec();
        // Spread the mismatches so they land in distinct blocks.
        for k in 0..allowed {
            let i = k * SEED_K + 2;
            damaged[i] = if damaged[i] == b'A' { b'C' } else { b'A' };
        }
        let mut read = b"TTTTTTTTTTTTTTTTTTTT".to_vec();
        let at = read.len();
        read.extend_from_slice(&damaged);

        let mut out = Vec::new();
        prep.candidates(&read, read.len() - adapter.len(), &mut Vec::<u64>::new(), &mut out);
        assert!(out.contains(&at), "pigeonhole broken: {at} not in {out:?}");
    }

    #[test]
    fn a_clean_read_yields_few_candidates() {
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let prep = PreparedAdapter::new(adapter);
        let read = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        let mut out = Vec::new();
        prep.candidates(read, read.len() - adapter.len(), &mut Vec::<u64>::new(), &mut out);
        assert!(out.is_empty(), "expected no candidates, got {out:?}");
    }

    #[test]
    fn ns_break_seeds_without_breaking_soundness() {
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let prep = PreparedAdapter::new(adapter);
        let mut read = b"ACGTACGTACGTACGTACGT".to_vec();
        let at = read.len();
        let mut damaged = adapter.to_vec();
        damaged[1] = b'N'; // kills block 0 only
        read.extend_from_slice(&damaged);
        let mut out = Vec::new();
        prep.candidates(&read, read.len() - adapter.len(), &mut Vec::<u64>::new(), &mut out);
        assert!(out.contains(&at), "later blocks should still anchor it");
    }

    #[test]
    fn a_short_adapter_disables_filtering() {
        let prep = PreparedAdapter::new(b"AGATCGGA");
        assert!(!prep.can_filter(), "8 bp cannot guarantee the pigeonhole");
    }

    #[test]
    fn candidates_come_back_sorted_and_unique() {
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let prep = PreparedAdapter::new(adapter);
        let mut read = Vec::new();
        for _ in 0..3 {
            read.extend_from_slice(adapter);
        }
        let mut out = Vec::new();
        prep.candidates(&read, read.len() - adapter.len(), &mut Vec::<u64>::new(), &mut out);
        let mut sorted = out.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(out, sorted);
    }
}
