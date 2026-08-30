//! Unit tests for [`polyx`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;

/// fastp's own `PolyX::test()`, expectations included.
#[test]
fn fastp_polyx_test_case() {
    let seq = b"ATTTTAAAAAAAAAATAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAT";
    let mut r = ReadRec::from_seq(seq);
    let trimmed = trim_poly_x(&mut r, 10);
    assert_eq!(r.bases(), b"ATTTT");
    let (base, removed) = trimmed.expect("a polyA tail should be found");
    assert_eq!(base, b'A');
    assert_eq!(removed, 51);
}

/// The cut lands on the leftmost G that the mismatch budget can still
/// reach, not on the first G of the visible run: walking back over the 12
/// Gs costs nothing, the `T` at 15 is the one allowed mismatch, and the
/// `G` at 14 then becomes the new `firstGPos`. So 14 bases go, not 12.
/// That is fastp's rule, and it is why polyG trimming can take a base or
/// two of real insert with it.
#[test]
fn trims_a_polyg_tail_back_to_the_reachable_g() {
    let seq = b"ACGTACGTACGTACGTGGGGGGGGGGGG";
    let mut r = ReadRec::from_seq(seq);
    let removed = trim_poly_g(&mut r, 10);
    assert_eq!(removed, 14);
    assert_eq!(r.bases(), b"ACGTACGTACGTAC");
}

#[test]
fn trims_a_clean_polyg_tail_exactly() {
    // No G adjacent to the run, so the cut is exactly the run.
    let seq = b"ACGTACGTACGTACAAGGGGGGGGGGGG";
    let mut r = ReadRec::from_seq(seq);
    assert_eq!(trim_poly_g(&mut r, 10), 12);
    assert_eq!(r.bases(), b"ACGTACGTACGTACAA");
}

#[test]
fn leaves_a_clean_read_alone() {
    let seq = b"ACGTACGTACGTACGTACGTACGTACGT";
    let mut r = ReadRec::from_seq(seq);
    assert_eq!(trim_poly_g(&mut r, 10), 0);
    assert_eq!(r.bases(), seq);
    let mut r2 = ReadRec::from_seq(seq);
    assert!(trim_poly_x(&mut r2, 10).is_none());
    assert_eq!(r2.bases(), seq);
}

#[test]
fn empty_read_is_safe() {
    let mut r = ReadRec::from_seq(b"");
    assert_eq!(trim_poly_g(&mut r, 10), 0);
    assert!(trim_poly_x(&mut r, 10).is_none());
}
