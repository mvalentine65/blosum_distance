//! Unit tests for [`adapter`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;

/// fastp's own `AdapterTrimmer::test()`, first case.
#[test]
fn fastp_trim_by_sequence_case() {
    let seq = b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG";
    let adapter = b"TTTTCCACGGGGATACTACTG";
    let mut r = ReadRec::from_seq(seq);
    let hit = trim_by_sequence(&mut r, &PreparedAdapter::new(adapter), 4, &mut AdapterScratch::default());
    assert!(hit.is_some());
    assert_eq!(r.bases(), b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAA");
}

/// fastp's own `AdapterTrimmer::test()`, multi-sequence case.
#[test]
fn fastp_trim_by_multi_sequences_case() {
    let seq = b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGGAAATTTCCCGGGAAATTTCCCGGGATCGATCGATCGATCGAATTCC";
    let adapters: Vec<PreparedAdapter> = [
        b"GCTAGCTAGCTAGCTA".as_slice(),
        b"AAATTTCCCGGGAAATTTCCCGGG".as_slice(),
        b"ATCGATCGATCGATCG".as_slice(),
        b"AATTCCGGAATTCCGG".as_slice(),
    ]
    .iter()
    .map(|a| PreparedAdapter::new(a))
    .collect();
    let mut r = ReadRec::from_seq(seq);
    trim_by_multi_sequences(&mut r, &adapters, &mut AdapterScratch::default());
    assert_eq!(
        r.bases(),
        b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG"
    );
}

#[test]
fn trims_a_real_truseq_tail() {
    // 40 bp of insert then the start of TruSeq R1
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let mut r = ReadRec::from_seq(seq);
    let hit = trim_by_sequence(&mut r, &PreparedAdapter::new(adapter), 4, &mut AdapterScratch::default()).expect("adapter should be found");
    assert_eq!(r.bases(), b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
    assert_eq!(hit.removed, 33);
}

#[test]
fn leaves_a_clean_read_alone() {
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let mut r = ReadRec::from_seq(seq);
    assert!(trim_by_sequence(&mut r, &PreparedAdapter::new(adapter), 4, &mut AdapterScratch::default()).is_none());
    assert_eq!(r.bases(), seq);
}

#[test]
fn an_adapter_only_read_is_emptied() {
    let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let mut r = ReadRec::from_seq(adapter);
    let hit = trim_by_sequence(&mut r, &PreparedAdapter::new(adapter), 4, &mut AdapterScratch::default()).expect("should match at zero");
    assert!(r.is_empty(), "left {:?}", r.bases());
    assert!(hit.removed > 0);
}

#[test]
fn match_req_scales_with_the_list() {
    assert_eq!(match_req_for(1), 4);
    assert_eq!(match_req_for(17), 5);
    assert_eq!(match_req_for(257), 6);
}
