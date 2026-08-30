//! Unit tests for [`detect`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;

/// fastp's own `NucleotideTree::test()`.
#[test]
fn fastp_nucleotide_tree_test_case() {
    let mut tree = NucleotideTree::new();
    for _ in 0..100 {
        tree.add_seq(b"AAAATTTT");
        tree.add_seq(b"AAAATTTTGGGG");
        tree.add_seq(b"AAAATTTTGGGGCCCC");
        tree.add_seq(b"AAAATTTTGGGGCCAA");
    }
    tree.add_seq(b"AAAATTTTGGGACCCC");
    let mut reached_leaf = true;
    let path = tree.dominant_path(&mut reached_leaf);
    assert_eq!(path, b"AAAATTTTGGGGCC".to_vec());
}

/// fastp's own `Evaluator::test()`.
#[test]
fn fastp_seq2int_roundtrip() {
    let s = b"ATCGATCGAT";
    let key = seq2int(s, 0, 10, -1);
    assert_eq!(int2seq(key as u32, 10), s.to_vec());
}

#[test]
fn rolling_and_fresh_keys_agree() {
    let s = b"ACGTACGTACGTACGT";
    let fresh = seq2int(s, 1, 10, -1);
    let rolled = seq2int(s, 1, 10, seq2int(s, 0, 10, -1));
    assert_eq!(fresh, rolled);
}

#[test]
fn an_n_poisons_the_key() {
    assert_eq!(seq2int(b"ACGTNCGTAC", 0, 10, -1), -1);
}

/// Synthetic library: 20k reads of varied insert, a third carrying TruSeq.
#[test]
fn detects_truseq_in_a_contaminated_library() {
    let truseq = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let bases = b"ACGT";
    let mut store: Vec<Vec<u8>> = Vec::new();
    let mut seed = 12345u64;
    let mut rand = move || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 33) as usize
    };
    for i in 0..20_000 {
        let insert_len = 60 + rand() % 60;
        let mut read: Vec<u8> = (0..insert_len).map(|_| bases[rand() % 4]).collect();
        if i % 3 == 0 {
            read.extend_from_slice(truseq);
        }
        read.resize(150, b'A');
        store.push(read);
    }
    let reads: Vec<&[u8]> = store.iter().map(|v| v.as_slice()).collect();
    match detect_adapter(&reads, 1) {
        Detected::Known { seq, name } => {
            assert_eq!(seq, truseq.to_vec());
            assert!(name.contains("TruSeq"), "got {name}");
        }
        other => panic!("expected the TruSeq entry, got {other:?}"),
    }
}

#[test]
fn a_clean_library_yields_nothing() {
    let bases = b"ACGT";
    let mut seed = 999u64;
    let mut rand = move || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 33) as usize
    };
    let store: Vec<Vec<u8>> = (0..20_000)
        .map(|_| (0..150).map(|_| bases[rand() % 4]).collect())
        .collect();
    let reads: Vec<&[u8]> = store.iter().map(|v| v.as_slice()).collect();
    assert_eq!(detect_adapter(&reads, 1), Detected::None);
}

#[test]
fn too_small_a_sample_refuses_to_guess() {
    let truseq = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let store: Vec<Vec<u8>> = (0..500)
        .map(|_| {
            let mut r = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
            r.extend_from_slice(truseq);
            r
        })
        .collect();
    let reads: Vec<&[u8]> = store.iter().map(|v| v.as_slice()).collect();
    assert_eq!(detect_adapter(&reads, 1), Detected::None);
}
