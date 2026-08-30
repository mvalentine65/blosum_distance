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
mod tests {
    use super::*;
    use crate::reads::overlap::{analyze, OverlapScratch};
    use crate::reads::seqops::reverse_complement;

    fn pair_for(insert: &[u8]) -> (Vec<u8>, Vec<u8>) {
        (insert.to_vec(), reverse_complement(insert))
    }

    #[test]
    fn a_bad_base_is_overwritten_by_its_good_mate() {
        let insert = b"ACGTTGCAAGGTCCATTGACGATCGGATCCAAGGTTCCAAGGTTCCAAGGTTCCAA";
        let (s1, s2) = pair_for(insert);
        let mut p = PairBuffers {
            qual1: vec![b'I'; s1.len()],   // phred 40
            qual2: vec![b'I'; s2.len()],
            seq1: s1,
            seq2: s2,
        };
        // Corrupt one R2 base and mark it terrible.
        let bad = 10;
        p.seq2[bad] = if p.seq2[bad] == b'A' { b'C' } else { b'A' };
        p.qual2[bad] = b'#'; // phred 2

        let mut scratch = OverlapScratch::default();
        let ov = analyze(&p.seq1, &p.seq2, 5, 30, 0.2, false, &mut scratch);
        assert!(ov.overlapped);
        let res = correct_by_overlap(&mut p, ov);
        assert_eq!(res.corrected, 1);
        assert_eq!(p.seq2, reverse_complement(insert));
    }

    #[test]
    fn two_good_bases_disagreeing_are_left_alone() {
        let insert = b"ACGTTGCAAGGTCCATTGACGATCGGATCCAAGGTTCCAAGGTTCCAAGGTTCCAA";
        let (s1, s2) = pair_for(insert);
        let mut p = PairBuffers {
            qual1: vec![b'I'; s1.len()],
            qual2: vec![b'I'; s2.len()],
            seq1: s1,
            seq2: s2,
        };
        let bad = 10;
        let before = p.seq2.clone();
        p.seq2[bad] = if p.seq2[bad] == b'A' { b'C' } else { b'A' };
        let corrupted = p.seq2.clone();

        let mut scratch = OverlapScratch::default();
        let ov = analyze(&p.seq1, &p.seq2, 5, 30, 0.2, false, &mut scratch);
        let res = correct_by_overlap(&mut p, ov);
        assert_eq!(res.corrected, 0);
        assert_eq!(res.uncorrected, 1);
        assert_eq!(p.seq2, corrupted, "must not have been repaired");
        assert_ne!(p.seq2, before);
    }

    #[test]
    fn merging_a_fully_overlapping_pair_returns_the_insert() {
        let insert = b"ACGTTGCAAGGTCCATTGACGATCGGATCCAAGGTTCCAAGGTTCCAAGGTTCCAA";
        let (s1, s2) = pair_for(insert);
        let q1 = vec![b'I'; s1.len()];
        let q2 = vec![b'I'; s2.len()];
        let mut scratch = OverlapScratch::default();
        let ov = analyze(&s1, &s2, 5, 30, 0.2, false, &mut scratch);
        let m = merge(&s1, &q1, &s2, &q2, ov).expect("overlapping pair merges");
        assert_eq!(m.seq, insert.to_vec());
        assert_eq!(m.qual.len(), m.seq.len());
    }

    #[test]
    fn merging_a_staggered_pair_spans_the_whole_fragment() {
        //  fragment longer than either read, so R2 contributes a tail
        let fragment = b"TTTTTTTTTTACGTTGCAAGGTCCATTGACGATCGGATCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGGGGGGGGG";
        let r1 = &fragment[..60];
        let r2 = reverse_complement(&fragment[fragment.len() - 60..]);
        let q1 = vec![b'I'; r1.len()];
        let q2 = vec![b'I'; r2.len()];
        let mut scratch = OverlapScratch::default();
        let ov = analyze(r1, &r2, 5, 30, 0.2, false, &mut scratch);
        assert!(ov.overlapped && ov.offset > 0, "offset {}", ov.offset);
        let m = merge(r1, &q1, &r2, &q2, ov).expect("staggered pair merges");
        assert_eq!(m.seq, fragment.to_vec());
    }

    #[test]
    fn non_overlapping_pairs_do_not_merge() {
        let s1 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let s2 = b"TTTTTTTTTTGGGGGGGGGGCCCCCCCCCCAAAAAAAAAA";
        let q = vec![b'I'; 40];
        let mut scratch = OverlapScratch::default();
        let ov = analyze(s1, s2, 2, 30, 0.2, false, &mut scratch);
        assert!(merge(s1, &q, s2, &q, ov).is_none());
    }
}
