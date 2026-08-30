//! Unit tests for [`filter`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;

fn opts() -> FilterOptions {
    FilterOptions {
        qual: Some(QualFilter::default()),
        length: Some(LengthFilter::default()),
        complexity: Some(ComplexityFilter::default()),
    }
}

#[test]
fn a_good_read_passes() {
    let seq = b"ACGTACGTACGTACGTACGT";
    let qual = b"IIIIIIIIIIIIIIIIIIII";
    let r = ReadRec::new(seq, qual);
    assert_eq!(pass_filter(&r, &opts()), FilterVerdict::Pass);
}

#[test]
fn low_quality_fails() {
    let seq = b"ACGTACGTACGTACGTACGT";
    let qual = b"####################";
    let r = ReadRec::new(seq, qual);
    assert_eq!(pass_filter(&r, &opts()), FilterVerdict::FailQuality);
}

#[test]
fn too_many_ns_fails() {
    let seq = b"ANNNNNNGTACGTACGTACG";
    let qual = b"IIIIIIIIIIIIIIIIIIII";
    let r = ReadRec::new(seq, qual);
    assert_eq!(pass_filter(&r, &opts()), FilterVerdict::FailNBase);
}

#[test]
fn short_reads_fail_on_length() {
    let seq = b"ACGTAC";
    let qual = b"IIIIII";
    let r = ReadRec::new(seq, qual);
    assert_eq!(pass_filter(&r, &opts()), FilterVerdict::FailLength);
}

#[test]
fn a_homopolymer_fails_complexity_when_enabled() {
    let seq = b"AAAAAAAAAAAAAAAAAAAA";
    let qual = b"IIIIIIIIIIIIIIIIIIII";
    let r = ReadRec::new(seq, qual);
    let mut o = opts();
    assert_eq!(pass_filter(&r, &o), FilterVerdict::Pass, "off by default");
    o.complexity = Some(ComplexityFilter { enabled: true, threshold: 0.3 });
    assert_eq!(pass_filter(&r, &o), FilterVerdict::FailComplexity);
}

#[test]
fn fixed_trims_apply_without_quality_cutting() {
    let seq = b"ACGTACGTACGT";
    let qual = b"IIIIIIIIIIII";
    let mut r = ReadRec::new(seq, qual);
    assert!(trim_and_cut(&mut r, 2, 3, &QualityCut::default()));
    assert_eq!(r.bases(), b"GTACGTA");
}

#[test]
fn cut_tail_removes_a_bad_tail() {
    //  20 good bases then 8 terrible ones
    let seq = b"ACGTACGTACGTACGTACGTACGTACGT";
    let qual = b"IIIIIIIIIIIIIIIIIIII########";
    let mut r = ReadRec::new(seq, qual);
    let cut = QualityCut {
        enabled_tail: true,
        window_size_tail: 4,
        quality_tail: 20,
        ..Default::default()
    };
    assert!(trim_and_cut(&mut r, 0, 0, &cut));
    assert!(r.len() <= 21 && r.len() >= 18, "kept {}", r.len());
    assert!(r.quals().iter().all(|&q| q >= b'5'), "kept a bad base");
}

/// A read with no good window anywhere is cut back to a single base rather
/// than dropped: the backward scan runs out of room, `t` lands one window
/// short of the start, and `rlen` comes out as 1 — which clears fastp's
/// `rlen <= 0` guard. The filters are what actually discard it, and the
/// quality check gets there first because fastp orders it before length.
#[test]
fn a_read_with_no_good_window_is_cut_to_one_base_then_filtered_out() {
    let seq = b"ACGTACGTACGTACGTACGT";
    let qual = b"####################";
    let mut r = ReadRec::new(seq, qual);
    let cut = QualityCut {
        enabled_tail: true,
        window_size_tail: 4,
        quality_tail: 30,
        ..Default::default()
    };
    assert!(trim_and_cut(&mut r, 0, 0, &cut));
    assert_eq!(r.len(), 1);
    assert_eq!(pass_filter(&r, &opts()), FilterVerdict::FailQuality);
}

/// `trim_and_cut` does return false when the fixed trims alone eat the read.
#[test]
fn fixed_trims_larger_than_the_read_drop_it() {
    let seq = b"ACGTAC";
    let qual = b"IIIIII";
    let mut r = ReadRec::new(seq, qual);
    assert!(!trim_and_cut(&mut r, 4, 4, &QualityCut::default()));
}

#[test]
fn fasta_input_skips_quality_stages() {
    let mut r = ReadRec::from_seq(b"ACGTACGTACGT");
    let cut = QualityCut { enabled_tail: true, window_size_tail: 4, quality_tail: 30, ..Default::default() };
    assert!(trim_and_cut(&mut r, 1, 1, &cut));
    assert_eq!(r.bases(), b"CGTACGTACG");

    let r2 = ReadRec::from_seq(b"ACGTACGTACGTACGTACGT");
    assert_eq!(pass_filter(&r2, &opts()), FilterVerdict::Pass);
}
