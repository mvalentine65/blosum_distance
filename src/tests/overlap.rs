//! Unit tests for [`overlap`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;
use crate::reads::seqops::reverse_complement;

/// fastp's own `OverlapAnalysis::test()`, both assertions.
#[test]
fn fastp_overlap_test_case() {
    let r1 = b"CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGC";
    let r2 = b"ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGTCC";
    let mut scratch = OverlapScratch::default();
    let ov = analyze(r1, r2, 2, 30, 0.2, false, &mut scratch);
    assert!(ov.overlapped);
    assert_eq!(ov.offset, 10);
    assert_eq!(ov.overlap_len, 79);
    assert_eq!(ov.diff, 1);
}

/// The 1.3.6 late-mismatch case: mismatches past base 50 no longer reject.
#[test]
fn fastp_late_mismatch_case() {
    let mut r1 = vec![b'A'; 50];
    r1.extend(std::iter::repeat(b'C').take(30));
    let mut rcr2 = vec![b'A'; 50];
    rcr2.extend(std::iter::repeat(b'G').take(30));
    let r2 = reverse_complement(&rcr2);

    let mut scratch = OverlapScratch::default();
    let ov = analyze(&r1, &r2, 0, 30, 0.0, false, &mut scratch);
    assert!(ov.overlapped);
    assert_eq!(ov.offset, 0);
    assert_eq!(ov.overlap_len, 80);
    assert_eq!(ov.diff, 30);
}

#[test]
fn short_insert_gives_a_negative_offset_and_trims_both_mates() {
    // 60 bp fragment sequenced by 80 bp reads: 20 bp of adapter on each end
    let insert = b"ACGTTGCAAGGTCCATTGACGATCGGATCCAAGGTTCCAAGGTTCCAAGGTTCCAA";
    let adapter1 = b"AGATCGGAAGAGCACACGTCT";
    let adapter2 = b"AGATCGGAAGAGCGTCGTGTA";

    let mut read1 = insert.to_vec();
    read1.extend_from_slice(adapter1);
    let mut read2 = reverse_complement(insert);
    read2.extend_from_slice(adapter2);

    let mut scratch = OverlapScratch::default();
    let ov = analyze(&read1, &read2, 5, 30, 0.2, false, &mut scratch);
    assert!(ov.overlapped, "a short insert must be detected");
    assert!(ov.offset < 0, "offset was {}", ov.offset);
    assert_eq!(ov.overlap_len, insert.len());

    let mut rec1 = ReadRec::from_seq(&read1);
    let mut rec2 = ReadRec::from_seq(&read2);
    let removed = trim_by_overlap(&mut rec1, &mut rec2, ov, 0, 0);
    assert!(removed.is_some());
    assert_eq!(rec1.bases(), &insert[..]);
    assert_eq!(rec2.bases(), &reverse_complement(insert)[..]);
}

#[test]
fn unrelated_mates_do_not_overlap() {
    let r1 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let r2 = b"TTTTTTTTTTGGGGGGGGGGCCCCCCCCCCAAAAAAAAAA";
    let mut scratch = OverlapScratch::default();
    let ov = analyze(r1, r2, 2, 30, 0.2, false, &mut scratch);
    assert!(!ov.overlapped);
}

#[test]
fn positive_offset_never_trims() {
    let ov = OverlapResult { overlapped: true, offset: 10, overlap_len: 79, diff: 1, has_gap: false };
    let seq = b"ACGTACGT";
    let mut r1 = ReadRec::from_seq(seq);
    let mut r2 = ReadRec::from_seq(seq);
    assert!(trim_by_overlap(&mut r1, &mut r2, ov, 0, 0).is_none());
    assert_eq!(r1.bases(), seq);
}
