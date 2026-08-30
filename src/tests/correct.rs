//! Unit tests for [`correct`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

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
