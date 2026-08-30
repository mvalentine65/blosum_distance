//! Sequence primitives.
//!
//! Ported from fastp 1.3.6 `src/simd.cpp` (MIT, (c) 2016 OpenGene). fastp
//! dispatches these through Google Highway; we stay in safe Rust — bio's
//! vectorised `hamming` for the exact count, a block-at-a-time loop LLVM
//! vectorises for the bounded one. Semantics match the Highway scalar tails
//! exactly.

use bio::alignment::distance::simd::hamming;

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
    // bio's vectorised hamming; the clamp above gives it the equal lengths it
    // requires.
    hamming(&a[..n], &b[..n]) as usize
}

/// Mismatches over the first `len` bytes, abandoning the count once it passes
/// `limit`. The return value is exact whenever it is at or under `limit` —
/// which is all any caller relies on: they either test `<= limit`, or use the
/// value only on the accepting branch, where no early exit can have happened.
#[inline]
pub fn count_mismatches_bounded(a: &[u8], b: &[u8], len: usize, limit: usize) -> usize {
    /// Bases compared before the budget is re-checked.
    const BLOCK: usize = 16;
    let n = len.min(a.len()).min(b.len());
    let (a, b) = (&a[..n], &b[..n]);
    let mut diff = 0usize;
    let mut i = 0usize;
    // A block at a time so the compare still vectorises; overshooting `limit`
    // within a block is what the contract above already allows.
    while i + BLOCK <= n {
        diff += count_mismatches_scalar(&a[i..i + BLOCK], &b[i..i + BLOCK]);
        if diff > limit {
            return diff;
        }
        i += BLOCK;
    }
    diff + count_mismatches_bounded_scalar(&a[i..], &b[i..], limit - diff)
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

/// The tail of the bounded scan, and the reference it is tested against.
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
/// * `low_qual` — bases scoring under `qualified_qual` (a raw phred+33 value,
///   widened so the caller's `+ 33` cannot wrap)
/// * `n_bases`  — literal `N` calls
/// * `total_qual` — sum of phred scores, i.e. each byte less 33
pub struct QualityMetrics {
    pub low_qual: usize,
    pub n_bases: usize,
    pub total_qual: i64,
}

pub fn count_quality_metrics(qual: &[u8], seq: &[u8], qualified_qual: i64) -> QualityMetrics {
    let n = qual.len().min(seq.len());
    let mut m = QualityMetrics { low_qual: 0, n_bases: 0, total_qual: 0 };
    for i in 0..n {
        let q = qual[i];
        m.total_qual += i64::from(q) - 33;
        if i64::from(q) < qualified_qual {
            m.low_qual += 1;
        }
        if seq[i] == b'N' {
            m.n_bases += 1;
        }
    }
    m
}

#[cfg(test)]
#[path = "../tests/seqops.rs"]
mod tests;
