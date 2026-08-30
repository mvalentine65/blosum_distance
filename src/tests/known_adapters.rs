//! Unit tests for [`known_adapters`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;

#[test]
fn table_loaded_and_ordered() {
    assert_eq!(KNOWN_ADAPTERS.len(), 234);
    let mut sorted = KNOWN_ADAPTERS.to_vec();
    sorted.sort_by(|a, b| a.0.cmp(b.0));
    assert_eq!(
        sorted.iter().map(|x| x.0).collect::<Vec<_>>(),
        KNOWN_ADAPTERS.iter().map(|x| x.0).collect::<Vec<_>>(),
        "table must be in std::map order to match fastp's tie break"
    );
}

#[test]
fn finds_the_truseq_pair() {
    let (adapter, name) =
        match_known_adapter(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAGGGGG").unwrap();
    assert_eq!(adapter, b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
    assert!(name.contains("TruSeq"), "got {name}");

    let (_, name2) =
        match_known_adapter(b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT").unwrap();
    assert!(name2.contains("TruSeq"), "got {name2}");
}

#[test]
fn a_short_or_unknown_candidate_matches_nothing() {
    assert!(match_known_adapter(b"AGATCGG").is_none());
    assert!(match_known_adapter(b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT").is_none());
}
