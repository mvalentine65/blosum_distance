//! Unit tests for [`seed`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;

#[test]
fn truseq_is_filterable_and_finds_its_own_offset() {
    let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let prep = PreparedAdapter::new(adapter);
    assert!(prep.can_filter());

    let mut read = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let at = read.len();
    read.extend_from_slice(adapter);

    let mut out = Vec::new();
    prep.candidates(&read, read.len() - adapter.len(), &mut Vec::<u64>::new(), &mut out);
    assert!(out.contains(&at), "offset {at} missing from {out:?}");
}

/// The guarantee that matters: an adapter carrying the maximum tolerated
/// number of mismatches is still proposed as a candidate.
#[test]
fn a_match_at_the_mismatch_limit_survives_the_filter() {
    let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let prep = PreparedAdapter::new(adapter);
    let allowed = adapter.len() / 8; // 4

    let mut damaged = adapter.to_vec();
    // Spread the mismatches so they land in distinct blocks.
    for k in 0..allowed {
        let i = k * SEED_K + 2;
        damaged[i] = if damaged[i] == b'A' { b'C' } else { b'A' };
    }
    let mut read = b"TTTTTTTTTTTTTTTTTTTT".to_vec();
    let at = read.len();
    read.extend_from_slice(&damaged);

    let mut out = Vec::new();
    prep.candidates(&read, read.len() - adapter.len(), &mut Vec::<u64>::new(), &mut out);
    assert!(out.contains(&at), "pigeonhole broken: {at} not in {out:?}");
}

#[test]
fn a_clean_read_yields_few_candidates() {
    let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let prep = PreparedAdapter::new(adapter);
    let read = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    let mut out = Vec::new();
    prep.candidates(read, read.len() - adapter.len(), &mut Vec::<u64>::new(), &mut out);
    assert!(out.is_empty(), "expected no candidates, got {out:?}");
}

#[test]
fn ns_break_seeds_without_breaking_soundness() {
    let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let prep = PreparedAdapter::new(adapter);
    let mut read = b"ACGTACGTACGTACGTACGT".to_vec();
    let at = read.len();
    let mut damaged = adapter.to_vec();
    damaged[1] = b'N'; // kills block 0 only
    read.extend_from_slice(&damaged);
    let mut out = Vec::new();
    prep.candidates(&read, read.len() - adapter.len(), &mut Vec::<u64>::new(), &mut out);
    assert!(out.contains(&at), "later blocks should still anchor it");
}

#[test]
fn a_short_adapter_disables_filtering() {
    let prep = PreparedAdapter::new(b"AGATCGGA");
    assert!(!prep.can_filter(), "8 bp cannot guarantee the pigeonhole");
}

#[test]
fn candidates_come_back_sorted_and_unique() {
    let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let prep = PreparedAdapter::new(adapter);
    let mut read = Vec::new();
    for _ in 0..3 {
        read.extend_from_slice(adapter);
    }
    let mut out = Vec::new();
    prep.candidates(&read, read.len() - adapter.len(), &mut Vec::<u64>::new(), &mut out);
    let mut sorted = out.clone();
    sorted.sort_unstable();
    sorted.dedup();
    assert_eq!(out, sorted);
}
