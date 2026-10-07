//! Sequence primitives.
//!
//! Ported from fastp 1.3.6 `src/simd.cpp` (MIT, (c) 2016 OpenGene). fastp
//! dispatches these through Google Highway; we stay in safe Rust — triple_accel's
//! vectorised `hamming` for the exact count, a block-at-a-time loop LLVM
//! vectorises for the bounded one. Semantics match the Highway scalar tails
//! exactly.

use crate::hamming;

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
    #[cfg(target_arch = "x86_64")]
    if src.len() >= 32 && avx2() {
        // SAFETY: AVX2 was just checked; the slices are equally long.
        unsafe { reverse_complement_avx2(src, dst) };
        return;
    }
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

/// `reverse_complement_into`, 32 bases a step.
///
/// # Safety
/// The CPU must support AVX2 (`avx2()`); `dst` is as long as `src`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn reverse_complement_avx2(src: &[u8], dst: &mut [u8]) {
    use std::arch::x86_64::*;
    let len = src.len();
    let mut i = 0;
    // SAFETY: a step reads 32 bytes at `i` and writes the 32 that mirror them.
    unsafe {
        let upper = _mm256_set1_epi8(0xDFu8 as i8);
        let (a, c, g, t) = (
            _mm256_set1_epi8(b'A' as i8),
            _mm256_set1_epi8(b'C' as i8),
            _mm256_set1_epi8(b'G' as i8),
            _mm256_set1_epi8(b'T' as i8),
        );
        let backwards = _mm256_setr_epi8(
            15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
            15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
        );
        while i + 32 <= len {
            let x = _mm256_and_si256(_mm256_loadu_si256(src.as_ptr().add(i) as *const __m256i), upper);
            let mut y = _mm256_set1_epi8(b'N' as i8);
            y = _mm256_blendv_epi8(y, t, _mm256_cmpeq_epi8(x, a));
            y = _mm256_blendv_epi8(y, a, _mm256_cmpeq_epi8(x, t));
            y = _mm256_blendv_epi8(y, g, _mm256_cmpeq_epi8(x, c));
            y = _mm256_blendv_epi8(y, c, _mm256_cmpeq_epi8(x, g));
            let y = _mm256_shuffle_epi8(y, backwards);
            let y = _mm256_permute2x128_si256::<1>(y, y);
            _mm256_storeu_si256(dst.as_mut_ptr().add(len - i - 32) as *mut __m256i, y);
            i += 32;
        }
    }
    for k in i..len {
        dst[len - 1 - k] = COMPLEMENT[src[k] as usize];
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
    // vectorised hamming; the clamp above gives it the equal lengths it
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
    #[cfg(target_arch = "x86_64")]
    if avx2() {
        // SAFETY: AVX2 and POPCNT were just checked; the slices are equally long.
        return unsafe { count_mismatches_bounded_avx2(a, b, limit) };
    }
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

/// Whether the AVX2 compare below may run.
#[inline]
pub fn avx2() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        is_x86_feature_detected!("avx2") && is_x86_feature_detected!("popcnt")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// Mismatches over the first `n` bytes of both slices, 32 bases a step.
///
/// A last short step is one masked compare when both slices can be read 32
/// bytes on from it, so a caller that pads its buffers never reaches the byte
/// loop.
///
/// # Safety
/// The CPU must support AVX2 and POPCNT (`avx2()`); both slices hold `n` bytes.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
pub unsafe fn count_mismatches_avx2(a: &[u8], b: &[u8], n: usize) -> usize {
    use std::arch::x86_64::{__m256i, _mm256_cmpeq_epi8, _mm256_loadu_si256, _mm256_movemask_epi8};
    debug_assert!(a.len() >= n && b.len() >= n);
    let same_at = |i: usize| -> u32 {
        // SAFETY: the callers below only pass an `i` with 32 readable bytes in both.
        unsafe {
            let x = _mm256_loadu_si256(a.as_ptr().add(i) as *const __m256i);
            let y = _mm256_loadu_si256(b.as_ptr().add(i) as *const __m256i);
            _mm256_movemask_epi8(_mm256_cmpeq_epi8(x, y)) as u32
        }
    };
    let (mut i, mut same) = (0usize, 0u32);
    while i + 32 <= n {
        same += same_at(i).count_ones();
        i += 32;
    }
    let rest = n - i;
    if rest > 0 {
        if a.len() >= i + 32 && b.len() >= i + 32 {
            same += (same_at(i) & ((1u32 << rest) - 1)).count_ones();
        } else {
            same += (rest - count_mismatches_scalar(&a[i..n], &b[i..n])) as u32;
        }
    }
    n - same as usize
}

/// `count_mismatches_bounded` over two equally long slices, 32 bases a step
/// and the budget checked after each.
///
/// # Safety
/// The CPU must support AVX2 and POPCNT (`avx2()`); `b` is as long as `a`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
unsafe fn count_mismatches_bounded_avx2(a: &[u8], b: &[u8], limit: usize) -> usize {
    use std::arch::x86_64::*;
    let n = a.len();
    let (mut i, mut diff) = (0usize, 0usize);
    // SAFETY: each step loads bytes that lie inside both slices.
    unsafe {
        while i + 32 <= n {
            let x = _mm256_loadu_si256(a.as_ptr().add(i) as *const __m256i);
            let y = _mm256_loadu_si256(b.as_ptr().add(i) as *const __m256i);
            diff += 32 - (_mm256_movemask_epi8(_mm256_cmpeq_epi8(x, y)) as u32).count_ones() as usize;
            if diff > limit {
                return diff;
            }
            i += 32;
        }
        if i + 16 <= n {
            let x = _mm_loadu_si128(a.as_ptr().add(i) as *const __m128i);
            let y = _mm_loadu_si128(b.as_ptr().add(i) as *const __m128i);
            diff += 16 - (_mm_movemask_epi8(_mm_cmpeq_epi8(x, y)) as u32).count_ones() as usize;
            if diff > limit {
                return diff;
            }
            i += 16;
        }
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

/// `count_quality_metrics` over the whole 32-base steps of the first `n`
/// bases, added to `m`. Returns the bases it covered.
///
/// # Safety
/// The CPU must support AVX2 and POPCNT (`avx2()`); both slices hold `n` bytes.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
unsafe fn quality_metrics_avx2(qual: &[u8], seq: &[u8], n: usize, qualified_qual: i64, m: &mut QualityMetrics) -> usize {
    use std::arch::x86_64::*;
    // a byte is under the threshold never when it is 0 or less, always when it is over 255
    let threshold = qualified_qual.clamp(0, 256);
    let (mut i, mut at_or_over, mut n_bases) = (0usize, 0u32, 0u32);
    let mut sums = [0u64; 4];
    // SAFETY: a step loads the 32 bytes at `i` of both slices, `i + 32 <= n`.
    unsafe {
        let floor = _mm256_set1_epi8(threshold.min(255) as u8 as i8);
        let (zero, base_n) = (_mm256_setzero_si256(), _mm256_set1_epi8(b'N' as i8));
        let mut sum = zero;
        while i + 32 <= n {
            let q = _mm256_loadu_si256(qual.as_ptr().add(i) as *const __m256i);
            let s = _mm256_loadu_si256(seq.as_ptr().add(i) as *const __m256i);
            sum = _mm256_add_epi64(sum, _mm256_sad_epu8(q, zero));
            at_or_over += (_mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_max_epu8(q, floor), q)) as u32).count_ones();
            n_bases += (_mm256_movemask_epi8(_mm256_cmpeq_epi8(s, base_n)) as u32).count_ones();
            i += 32;
        }
        _mm256_storeu_si256(sums.as_mut_ptr() as *mut __m256i, sum);
    }
    m.total_qual += sums.iter().sum::<u64>() as i64 - 33 * i as i64;
    m.n_bases += n_bases as usize;
    m.low_qual += match threshold {
        0 => 0,
        256 => i,
        _ => i - at_or_over as usize,
    };
    i
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
    let mut done = 0;
    #[cfg(target_arch = "x86_64")]
    if n >= 32 && avx2() {
        // SAFETY: AVX2 and POPCNT were just checked; both slices hold `n` bytes.
        done = unsafe { quality_metrics_avx2(qual, seq, n, qualified_qual, &mut m) };
    }
    for i in done..n {
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
