//! Sequence primitives.
//!
//! Ported from fastp 1.3.6 `src/simd.cpp` (MIT, (c) 2016 OpenGene). fastp
//! dispatches these through Google Highway; we keep the scalar forms, which
//! LLVM vectorises on its own and which carry no build dependency. Semantics
//! match the Highway scalar tails exactly.

/// Complement lookup, built the way fastp's `kComplement` table works: A/a,
/// C/c, T/t, G/g map across cases and everything else falls to N. A table
/// beats a match here because this runs once per base of every mate.
const COMPLEMENT: [u8; 256] = {
    let mut t = [b'N'; 256];
    t[b'A' as usize] = b'T';
    t[b'a' as usize] = b'T';
    t[b'T' as usize] = b'A';
    t[b't' as usize] = b'A';
    t[b'C' as usize] = b'G';
    t[b'c' as usize] = b'G';
    t[b'G' as usize] = b'C';
    t[b'g' as usize] = b'C';
    t
};

#[inline]
pub fn complement(base: u8) -> u8 {
    COMPLEMENT[base as usize]
}

/// `dst` receives the reverse complement of `src`. `dst` must be `src.len()`.
pub fn reverse_complement_into(src: &[u8], dst: &mut [u8]) {
    debug_assert_eq!(src.len(), dst.len());
    // Walk both ends toward the middle so each iteration is a pair of
    // independent table lookups rather than one reversed dependent stride.
    let len = src.len();
    let (mut i, mut j) = (0usize, len);
    while i + 1 < j {
        j -= 1;
        dst[len - 1 - i] = COMPLEMENT[src[i] as usize];
        dst[len - 1 - j] = COMPLEMENT[src[j] as usize];
        i += 1;
    }
    if i < j {
        dst[len - 1 - i] = COMPLEMENT[src[i] as usize];
    }
}

/// Allocating form. Used by the unwired merge path and by tests; the hot
/// path uses `reverse_complement_into` with a reused buffer.
#[allow(dead_code)]
pub fn reverse_complement(src: &[u8]) -> Vec<u8> {
    let mut out = vec![0u8; src.len()];
    reverse_complement_into(src, &mut out);
    out
}

/// Total mismatches over the first `len` bytes of both slices.
#[inline]
pub fn count_mismatches(a: &[u8], b: &[u8], len: usize) -> usize {
    let n = len.min(a.len()).min(b.len());
    #[cfg(target_arch = "x86_64")]
    {
        // SSE2 is part of the x86_64 baseline, so no runtime probe is needed.
        return unsafe { count_mismatches_sse2(&a[..n], &b[..n]) };
    }
    #[cfg(not(target_arch = "x86_64"))]
    count_mismatches_scalar(&a[..n], &b[..n])
}

/// Mismatches over the first `len` bytes, abandoning the count once it passes
/// `limit`. The return value is exact whenever it is at or under `limit` —
/// which is all any caller relies on: they either test `<= limit`, or use the
/// value only on the accepting branch, where no early exit can have happened.
#[inline]
pub fn count_mismatches_bounded(a: &[u8], b: &[u8], len: usize, limit: usize) -> usize {
    let n = len.min(a.len()).min(b.len());
    #[cfg(target_arch = "x86_64")]
    {
        return unsafe { count_mismatches_bounded_sse2(&a[..n], &b[..n], limit) };
    }
    #[cfg(not(target_arch = "x86_64"))]
    count_mismatches_bounded_scalar(&a[..n], &b[..n], limit)
}

#[inline]
fn count_mismatches_scalar(a: &[u8], b: &[u8]) -> usize {
    let mut diff = 0;
    for i in 0..a.len() {
        if a[i] != b[i] {
            diff += 1;
        }
    }
    diff
}

/// Fallback for non-x86_64, and the reference the SIMD path is tested against.
#[allow(dead_code)]
#[inline]
fn count_mismatches_bounded_scalar(a: &[u8], b: &[u8], limit: usize) -> usize {
    let mut diff = 0;
    for i in 0..a.len() {
        if a[i] != b[i] {
            diff += 1;
            if diff > limit {
                return diff;
            }
        }
    }
    diff
}

#[cfg(target_arch = "x86_64")]
#[inline]
unsafe fn count_mismatches_sse2(a: &[u8], b: &[u8]) -> usize {
    use std::arch::x86_64::*;
    let n = a.len();
    let mut i = 0;
    let mut diff = 0usize;
    while i + 16 <= n {
        let va = _mm_loadu_si128(a.as_ptr().add(i).cast());
        let vb = _mm_loadu_si128(b.as_ptr().add(i).cast());
        // movemask sets a bit per *equal* byte, so the mismatches are the
        // zero bits within the low 16.
        let eq = _mm_movemask_epi8(_mm_cmpeq_epi8(va, vb)) as u32;
        diff += 16 - (eq & 0xFFFF).count_ones() as usize;
        i += 16;
    }
    diff + count_mismatches_scalar(&a[i..], &b[i..])
}

#[cfg(target_arch = "x86_64")]
#[inline]
unsafe fn count_mismatches_bounded_sse2(a: &[u8], b: &[u8], limit: usize) -> usize {
    use std::arch::x86_64::*;
    let n = a.len();
    let mut i = 0;
    let mut diff = 0usize;
    while i + 16 <= n {
        let va = _mm_loadu_si128(a.as_ptr().add(i).cast());
        let vb = _mm_loadu_si128(b.as_ptr().add(i).cast());
        let eq = _mm_movemask_epi8(_mm_cmpeq_epi8(va, vb)) as u32;
        diff += 16 - (eq & 0xFFFF).count_ones() as usize;
        if diff > limit {
            // Overshoots the scalar count, but every caller only asks whether
            // the limit was blown.
            return diff;
        }
        i += 16;
    }
    for j in i..n {
        if a[j] != b[j] {
            diff += 1;
            if diff > limit {
                return diff;
            }
        }
    }
    diff
}

/// Count of positions where a base differs from the one after it. fastp uses
/// this over the whole read as its complexity measure.
#[inline]
pub fn count_adjacent_diffs(data: &[u8]) -> usize {
    if data.len() <= 1 {
        return 0;
    }
    let mut diff = 0;
    for i in 0..data.len() - 1 {
        if data[i] != data[i + 1] {
            diff += 1;
        }
    }
    diff
}

/// Per-read quality tallies, all three in one pass.
///
/// * `low_qual` — bases scoring under `qualified_qual` (a raw phred+33 byte)
/// * `n_bases`  — literal `N` calls
/// * `total_qual` — sum of phred scores, i.e. each byte less 33
pub struct QualityMetrics {
    pub low_qual: usize,
    pub n_bases: usize,
    pub total_qual: i64,
}

pub fn count_quality_metrics(qual: &[u8], seq: &[u8], qualified_qual: u8) -> QualityMetrics {
    let n = qual.len().min(seq.len());
    let mut m = QualityMetrics { low_qual: 0, n_bases: 0, total_qual: 0 };
    for i in 0..n {
        let q = qual[i];
        m.total_qual += i64::from(q) - 33;
        if q < qualified_qual {
            m.low_qual += 1;
        }
        if seq[i] == b'N' {
            m.n_bases += 1;
        }
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_matches_fastp_table() {
        assert_eq!(reverse_complement(b"ATCGN"), b"NCGAT".to_vec());
        // lowercase shares the nibble, unknown bases fall to N
        assert_eq!(reverse_complement(b"atcg"), b"CGAT".to_vec());
        assert_eq!(reverse_complement(b"AXA"), b"TNT".to_vec());
    }

    #[test]
    fn mismatch_counters() {
        assert_eq!(count_mismatches(b"AAAA", b"AATA", 4), 1);
        assert!(count_mismatches_bounded(b"AAAA", b"TTTT", 4, 1) > 1);
        assert!(count_mismatches_bounded(b"AAAA", b"TTTT", 4, 0) > 0);
        assert_eq!(count_mismatches_bounded(b"AAAA", b"AAAA", 4, 0), 0);
    }

    /// The vector path must agree with the scalar one everywhere the result is
    /// load-bearing: exact totals, and exact counts up to the limit. Lengths
    /// span the 16-byte block boundary in both directions.
    #[test]
    fn simd_agrees_with_scalar() {
        let mut seed = 0x12345678u64;
        let mut rand = move || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            (seed >> 33) as usize
        };
        let bases = b"ACGT";
        for len in 0..80usize {
            for _ in 0..20 {
                let a: Vec<u8> = (0..len).map(|_| bases[rand() % 4]).collect();
                let mut b = a.clone();
                for _ in 0..rand() % 6 {
                    if len > 0 {
                        let i = rand() % len;
                        b[i] = bases[rand() % 4];
                    }
                }
                let want = count_mismatches_scalar(&a, &b);
                assert_eq!(count_mismatches(&a, &b, len), want, "len {len}");
                for limit in 0..6usize {
                    let got = count_mismatches_bounded(&a, &b, len, limit);
                    let want_b = count_mismatches_bounded_scalar(&a, &b, limit);
                    if want_b <= limit {
                        assert_eq!(got, want_b, "exact under limit, len {len} limit {limit}");
                    } else {
                        assert!(got > limit, "must still reject, len {len} limit {limit}");
                    }
                }
            }
        }
    }

    #[test]
    fn adjacent_diffs_is_the_complexity_measure() {
        assert_eq!(count_adjacent_diffs(b"AAAA"), 0);
        assert_eq!(count_adjacent_diffs(b"ATAT"), 3);
        assert_eq!(count_adjacent_diffs(b""), 0);
        assert_eq!(count_adjacent_diffs(b"A"), 0);
    }

    #[test]
    fn quality_metrics_sum_phred_not_bytes() {
        // 'I' is phred 40, '#' is phred 2
        let m = count_quality_metrics(b"II#I", b"ACNG", b'5');
        assert_eq!(m.total_qual, 40 + 40 + 2 + 40);
        assert_eq!(m.low_qual, 1);
        assert_eq!(m.n_bases, 1);
    }
}
