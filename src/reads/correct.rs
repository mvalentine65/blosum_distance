#![allow(dead_code)]
//! Overlap base correction, and merging a pair into one record.
//!
//! **Not yet wired into the ingest chain.** This is the optional read-merging
//! half of the port: it works and is tested, but nothing in `build_table`
//! calls it, because merging changes what a read *is* — one record instead of
//! two, a different length distribution, and dupe counts that mean something
//! different downstream. Switching it on is a deliberate decision, not a
//! default. The `allow` above keeps that state from generating noise.
//!
//! Ported from fastp 1.3.6 `src/basecorrector.cpp` and the `merge` half of
//! `src/overlapanalysis.cpp` (MIT, (c) 2016 OpenGene).
//!
//! Where the mates disagree inside the overlap, one of them is usually right:
//! a high-quality call opposite a low-quality one wins, and the loser is
//! overwritten. Anything less lopsided is left alone.

use crate::reads::overlap::OverlapResult;
use crate::reads::seqops::complement;

/// fastp's thresholds: correct only when one side is clearly good and the
/// other clearly bad.
const GOOD_QUAL: u8 = 30 + 33;
const BAD_QUAL: u8 = 14 + 33;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CorrectionResult {
    pub corrected: usize,
    pub uncorrected: usize,
}

/// A pair whose bases may be rewritten. Correction mutates both mates, so
/// unlike the trimmers this stage needs owned buffers.
pub struct PairBuffers {
    pub seq1: Vec<u8>,
    pub qual1: Vec<u8>,
    pub seq2: Vec<u8>,
    pub qual2: Vec<u8>,
}

/// Fix disagreements inside the overlap, fastp's `BaseCorrector::correctByOverlapAnalysis`.
pub fn correct_by_overlap(p: &mut PairBuffers, ov: OverlapResult) -> CorrectionResult {
    let mut out = CorrectionResult::default();
    if !ov.overlapped || ov.diff == 0 || ov.has_gap {
        return out;
    }
    let ol = ov.overlap_len;
    let start1 = ov.offset.max(0) as usize;
    let len2 = p.seq2.len();
    let neg = (-ov.offset).max(0) as usize;
    if len2 < neg + 1 {
        return out;
    }
    let start2 = len2 - neg - 1;

    for i in 0..ol {
        let p1 = start1 + i;
        if p1 >= p.seq1.len() || i > start2 {
            break;
        }
        let p2 = start2 - i;

        if p.seq1[p1] != complement(p.seq2[p2]) {
            if p.qual1[p1] >= GOOD_QUAL && p.qual2[p2] <= BAD_QUAL {
                p.seq2[p2] = complement(p.seq1[p1]);
                p.qual2[p2] = p.qual1[p1];
                out.corrected += 1;
            } else if p.qual2[p2] >= GOOD_QUAL && p.qual1[p1] <= BAD_QUAL {
                p.seq1[p1] = complement(p.seq2[p2]);
                p.qual1[p1] = p.qual2[p2];
                out.corrected += 1;
            } else {
                out.uncorrected += 1;
            }
        }
    }
    out
}

/// One merged record.
pub struct Merged {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    /// Bases contributed by R1 and by the tail of R2, as fastp's name suffix.
    pub from_r1: usize,
    pub from_r2: usize,
}

/// Splice an overlapping pair into a single record, fastp's
/// `OverlapAnalysis::merge`. Returns `None` when the mates do not overlap.
pub fn merge(
    seq1: &[u8],
    qual1: &[u8],
    seq2: &[u8],
    qual2: &[u8],
    ov: OverlapResult,
) -> Option<Merged> {
    if !ov.overlapped {
        return None;
    }
    let ol = ov.overlap_len;
    let len1 = (ol as isize + ov.offset.max(0)).max(0) as usize;
    let len1 = len1.min(seq1.len());
    let len2 = if ov.offset > 0 { seq2.len().saturating_sub(ol) } else { 0 };

    let mut seq = seq1[..len1].to_vec();
    let mut qual = qual1[..len1.min(qual1.len())].to_vec();
    // Only a positive offset splices R2's tail on, so nothing is reversed until
    // it is known to be needed.
    if ov.offset > 0 && len2 > 0 && ol + len2 <= seq2.len() {
        let rc2 = crate::reads::seqops::reverse_complement(seq2);
        let rq2: Vec<u8> = qual2.iter().rev().copied().collect();
        seq.extend_from_slice(&rc2[ol..ol + len2]);
        qual.extend_from_slice(&rq2[ol..ol + len2]);
    }
    Some(Merged { seq, qual, from_r1: len1, from_r2: len2 })
}

#[cfg(test)]
#[path = "../tests/correct.rs"]
mod tests;
