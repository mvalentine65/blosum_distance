//! Unit tests for [`read`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;

#[test]
fn window_trims_without_copying() {
    let seq = b"ACGTACGT";
    let qual = b"IIIIIIII";
    let mut r = ReadRec::new(seq, qual);
    r.resize(6);
    assert_eq!(r.bases(), b"ACGTAC");
    r.trim_front(2);
    assert_eq!(r.bases(), b"GTAC");
    assert_eq!(r.quals(), b"IIII");
    assert_eq!(r.front_trimmed(), 2);
    assert_eq!(r.len(), 4);
}

#[test]
fn resize_past_the_end_is_a_no_op() {
    let mut r = ReadRec::from_seq(b"ACGT");
    r.resize(99);
    assert_eq!(r.bases(), b"ACGT");
    r.resize(0);
    assert!(r.is_empty());
}

#[test]
fn fasta_reads_carry_no_quality() {
    let r = ReadRec::from_seq(b"ACGT");
    assert!(!r.has_qual());
    assert_eq!(r.quals(), b"");
}
