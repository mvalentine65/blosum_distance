//! Unit tests for [`pipeline`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;
use crate::reads::seqops::reverse_complement;

fn truseq_r1() -> Vec<u8> {
    b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec()
}

#[test]
fn a_single_read_loses_its_adapter_and_passes() {
    let insert = b"ACGTTGCAAGGTCCATTGACGATCGGATCCAAGGTTCCAAGG";
    let mut seq = insert.to_vec();
    seq.extend_from_slice(&truseq_r1());
    let qual = vec![b'I'; seq.len()];

    let opts = TrimOptions { adapter_r1: Some(PreparedAdapter::new(&truseq_r1())), ..Default::default() };
    let mut stats = TrimStats::default();
    let kept = process_single(&seq, &qual, &opts, &mut stats, &mut TrimScratch::default()).expect("read should survive");
    assert_eq!(kept.bases(), &insert[..]);
    assert_eq!(stats.adapter_trimmed_reads, 1);
    assert_eq!(stats.reads_passed, 1);
}

#[test]
fn a_low_quality_read_is_dropped_and_counted() {
    let seq = b"ACGTTGCAAGGTCCATTGACGATCGGATCCAAGGTTCCAAGG";
    let qual = vec![b'#'; seq.len()];
    let opts = TrimOptions::default();
    let mut stats = TrimStats::default();
    assert!(process_single(seq, &qual, &opts, &mut stats, &mut TrimScratch::default()).is_none());
    assert_eq!(stats.failed_quality, 1);
    assert_eq!(stats.reads_passed, 0);
}

#[test]
fn a_short_insert_pair_is_trimmed_by_overlap_without_knowing_the_adapter() {
    let insert = b"ACGTTGCAAGGTCCATTGACGATCGGATCCAAGGTTCCAAGGTTCCAAGGTTCCAA";
    let mut read1 = insert.to_vec();
    read1.extend_from_slice(b"AGATCGGAAGAGCACACGTCT");
    let mut read2 = reverse_complement(insert);
    read2.extend_from_slice(b"AGATCGGAAGAGCGTCGTGTA");
    let q1 = vec![b'I'; read1.len()];
    let q2 = vec![b'I'; read2.len()];

    // No adapter list at all — the overlap has to carry it.
    let opts = TrimOptions::default();
    let mut stats = TrimStats::default();
    let mut scratch = TrimScratch::default();
    let out = process_pair(&read1, &q1, &read2, &q2, &opts, &mut stats, &mut scratch);
    assert_eq!(out.r1.expect("r1 kept").bases(), &insert[..]);
    assert_eq!(out.r2.expect("r2 kept").bases(), &reverse_complement(insert)[..]);
    assert_eq!(stats.adapter_trimmed_bases, 42);
}

#[test]
fn an_adapter_dimer_pair_is_dropped() {
    // Nothing but adapter in both mates.
    let a1 = truseq_r1();
    let a2 = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_vec();
    let q1 = vec![b'I'; a1.len()];
    let q2 = vec![b'I'; a2.len()];
    let opts = TrimOptions {
        adapter_r1: Some(PreparedAdapter::new(&a1)),
        adapter_r2: Some(PreparedAdapter::new(&a2)),
        overlap: false,
        ..Default::default()
    };
    let mut stats = TrimStats::default();
    let mut scratch = TrimScratch::default();
    let out = process_pair(&a1, &q1, &a2, &q2, &opts, &mut stats, &mut scratch);
    assert!(out.r1.is_none() && out.r2.is_none());
    assert_eq!(stats.failed_adapter_dimer, 2);
}

#[test]
fn a_clean_pair_survives_untouched() {
    let insert = b"ACGTTGCAAGGTCCATTGACGATCGGATCCAAGGTTCCAAGGTTCCAAGGTTCCAA";
    let read1 = insert.to_vec();
    let read2 = reverse_complement(insert);
    let q1 = vec![b'I'; read1.len()];
    let q2 = vec![b'I'; read2.len()];
    let opts = TrimOptions { adapter_r1: Some(PreparedAdapter::new(&truseq_r1())), ..Default::default() };
    let mut stats = TrimStats::default();
    let mut scratch = TrimScratch::default();
    let out = process_pair(&read1, &q1, &read2, &q2, &opts, &mut stats, &mut scratch);
    assert_eq!(out.r1.expect("kept").bases(), &read1[..]);
    assert_eq!(out.r2.expect("kept").bases(), &read2[..]);
    assert_eq!(stats.adapter_trimmed_bases, 0);
    assert_eq!(stats.reads_passed, 2);
}

#[test]
fn everything_off_is_a_pass_through() {
    let seq = b"ACGTTGCAAGGTCCATTGACGATCGGATCCAAGGTTCCAAGG";
    let qual = vec![b'I'; seq.len()];
    let opts = TrimOptions {
        adapter_trimming: false,
        trim_poly_g: false,
        trim_poly_x: false,
        overlap: false,
        filters: FilterOptions::default(),
        ..Default::default()
    };
    let mut stats = TrimStats::default();
    let kept = process_single(seq, &qual, &opts, &mut stats, &mut TrimScratch::default()).expect("kept");
    assert_eq!(kept.bases(), &seq[..]);
}
