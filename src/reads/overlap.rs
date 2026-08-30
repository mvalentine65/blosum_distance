//! Paired-end overlap analysis, and adapter trimming derived from it.
//!
//! Ported from fastp 1.3.6 `src/overlapanalysis.cpp` (MIT, (c) 2016 OpenGene),
//! which credits AfterQC for the original Python.
//!
//! When the insert is shorter than the read length both mates run off the end
//! of the fragment and into the adapter. Sliding R1 against the reverse
//! complement of R2 finds where the fragment actually stops, which locates the
//! adapter without needing to know its sequence — the reason this is the
//! better trimmer on a library with a short insert.
//!
//! Measured note on optimisation: the pigeonhole seed filter that pays off in
//! `seed.rs` does **not** transfer here. The slide tests 162 offsets per pair
//! on this library and returns at the answer, while a filter would have to
//! sweep ~150 positions of r1 for the forward direction plus ~150 of the
//! reverse complement — and with the inner comparison vectorised, one offset
//! test costs about what one filter position costs. The filter's fixed
//! full-pass cost cancels its benefit, so the slide stays as it is.

use crate::reads::matcher::{diff_with_one_insertion, MatchScratch};
use crate::reads::read::ReadRec;
use crate::reads::seqops::{count_mismatches, count_mismatches_bounded, reverse_complement_into};

/// fastp's `complete_compare_require`: past this many bases the mismatch limit
/// stops gating acceptance. 1.3.6 behaviour, and the reason it accepts some
/// overlaps 0.23.4 rejects.
const COMPLETE_COMPARE_REQUIRE: usize = 50;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct OverlapResult {
    pub overlapped: bool,
    /// How far R1 is shifted against the reverse complement of R2. Negative
    /// means the fragment is shorter than the read, i.e. adapter was sequenced.
    pub offset: isize,
    pub overlap_len: usize,
    pub diff: usize,
    pub has_gap: bool,
}

/// Reusable scratch so the reverse complement is not reallocated per pair.
#[derive(Default)]
pub struct OverlapScratch {
    rcr2: Vec<u8>,
    matcher: MatchScratch,
}

/// Accept-or-not for one candidate offset, fastp's `acceptNoGapOverlap`.
///
/// The mismatch limit is enforced only over the first `COMPLETE_COMPARE_REQUIRE`
/// bases; beyond that the overlap is accepted and the remaining mismatches are
/// merely counted. fastp's own unit test asserts an overlap with 30 mismatches
/// is accepted this way.
fn accept_no_gap(a: &[u8], b: &[u8], len: usize, diff_limit: usize) -> Option<usize> {
    let len = len.min(a.len()).min(b.len());
    let protected_prefix = len.min(COMPLETE_COMPARE_REQUIRE);
    let mut mismatch = count_mismatches_bounded(a, b, protected_prefix, diff_limit);
    if mismatch > diff_limit {
        return None;
    }
    if len > protected_prefix {
        // The prefix count is exact on this branch, so only the rest is left.
        mismatch += count_mismatches(
            &a[protected_prefix..],
            &b[protected_prefix..],
            len - protected_prefix,
        );
    }
    Some(mismatch)
}

/// Slide R1 against the reverse complement of R2 and report the best overlap.
pub fn analyze(
    r1: &[u8],
    r2: &[u8],
    diff_limit: usize,
    overlap_require: usize,
    diff_percent_limit: f64,
    allow_gap: bool,
    scratch: &mut OverlapScratch,
) -> OverlapResult {
    let len1 = r1.len();
    let len2 = r2.len();
    if len1 == 0 || len2 == 0 {
        return OverlapResult::default();
    }
    // revcomp overwrites every byte it is given, so the buffer only ever grows.
    if scratch.rcr2.len() < len2 {
        scratch.rcr2.resize(len2, 0);
    }
    reverse_complement_into(r2, &mut scratch.rcr2[..len2]);
    let str2 = &scratch.rcr2[..len2];

    let limit_for = |overlap_len: usize| -> usize {
        diff_limit.min((overlap_len as f64 * diff_percent_limit) as usize)
    };

    // Forward, no gap: the fragment is at least as long as the read.
    let mut offset: isize = 0;
    while offset < len1 as isize - overlap_require as isize {
        let overlap_len = (len1 - offset as usize).min(len2);
        let lim = limit_for(overlap_len);
        if let Some(diff) = accept_no_gap(&r1[offset as usize..], str2, overlap_len, lim) {
            return OverlapResult {
                overlapped: true,
                offset,
                overlap_len,
                diff,
                has_gap: false,
            };
        }
        offset += 1;
    }

    // Reverse, no gap: the fragment is shorter than the read, so adapter was
    // sequenced on both mates.
    offset = 0;
    while offset > -(len2 as isize - overlap_require as isize) {
        let overlap_len = len1.min(len2 - offset.unsigned_abs());
        let lim = limit_for(overlap_len);
        if let Some(diff) = accept_no_gap(r1, &str2[offset.unsigned_abs()..], overlap_len, lim) {
            return OverlapResult {
                overlapped: true,
                offset,
                overlap_len,
                diff,
                has_gap: false,
            };
        }
        offset -= 1;
    }

    if allow_gap {
        // Forward with one gap.
        offset = 0;
        while offset < len1 as isize - overlap_require as isize {
            let overlap_len = (len1 - offset as usize).min(len2);
            let lim = limit_for(overlap_len) as i32;
            let a = &r1[offset as usize..];
            let mut diff = diff_with_one_insertion(&mut scratch.matcher, a, str2, overlap_len - 1, lim);
            if diff < 0 || diff > lim {
                diff = diff_with_one_insertion(&mut scratch.matcher, str2, a, overlap_len - 1, lim);
            }
            if diff >= 0 && diff <= lim {
                return OverlapResult {
                    overlapped: true,
                    offset,
                    overlap_len,
                    diff: diff as usize,
                    has_gap: true,
                };
            }
            offset += 1;
        }

        // Reverse with one gap.
        offset = 0;
        while offset > -(len2 as isize - overlap_require as isize) {
            let overlap_len = len1.min(len2 - offset.unsigned_abs());
            let lim = limit_for(overlap_len) as i32;
            let b = &str2[offset.unsigned_abs()..];
            let mut diff = diff_with_one_insertion(&mut scratch.matcher, r1, b, overlap_len - 1, lim);
            if diff < 0 || diff > lim {
                diff = diff_with_one_insertion(&mut scratch.matcher, b, r1, overlap_len - 1, lim);
            }
            if diff >= 0 && diff <= lim {
                return OverlapResult {
                    overlapped: true,
                    offset,
                    overlap_len,
                    diff: diff as usize,
                    has_gap: true,
                };
            }
            offset -= 1;
        }
    }

    OverlapResult::default()
}

/// Trim both mates back to the fragment the overlap identified.
///
/// Only a negative offset means adapter was sequenced. Each mate is clipped
/// through a `min`, so a mate already at or inside its boundary is left alone
/// — which is why a pair often shows a change on one side only.
///
/// Returns `(r1 removed, r2 removed)`.
pub fn trim_by_overlap(
    r1: &mut ReadRec<'_>,
    r2: &mut ReadRec<'_>,
    ov: OverlapResult,
    front_trimmed1: usize,
    front_trimmed2: usize,
) -> Option<(usize, usize)> {
    if !ov.overlapped || ov.offset >= 0 {
        return None;
    }
    let ol = ov.overlap_len;
    let len1 = r1.len().min(ol + front_trimmed2);
    let len2 = r2.len().min(ol + front_trimmed1);
    let removed1 = r1.len() - len1;
    let removed2 = r2.len() - len2;
    r1.resize(len1);
    r2.resize(len2);
    Some((removed1, removed2))
}

#[cfg(test)]
#[path = "../tests/overlap.rs"]
mod tests;
