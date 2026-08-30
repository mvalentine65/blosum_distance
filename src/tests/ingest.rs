//! Unit tests for [`ingest`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;
use crate::reads::seqops::reverse_complement;

fn p(s: &str) -> PathBuf {
    PathBuf::from(s)
}

#[test]
fn pairs_the_real_library_naming() {
    let paths = vec![
        p("/in/DDM8637_CKDL230042833-1A_HGC33DSX7_L2_1.fq.gz"),
        p("/in/DDM8637_CKDL230042833-1A_HGC33DSX7_L2_2.fq.gz"),
    ];
    assert_eq!(
        group_inputs(&paths),
        vec![InputGroup::Paired(paths[0].clone(), paths[1].clone())]
    );
}

#[test]
fn pairs_r1_r2_and_orders_them() {
    // R2 listed first: the group must still come back R1-then-R2.
    let paths = vec![p("/in/lib_R2.fastq"), p("/in/lib_R1.fastq")];
    assert_eq!(
        group_inputs(&paths),
        vec![InputGroup::Paired(p("/in/lib_R1.fastq"), p("/in/lib_R2.fastq"))]
    );
}

#[test]
fn an_unmatched_file_stays_single() {
    let paths = vec![p("/in/lib_1.fq"), p("/in/other.fa")];
    assert_eq!(
        group_inputs(&paths),
        vec![
            InputGroup::Single(p("/in/lib_1.fq")),
            InputGroup::Single(p("/in/other.fa")),
        ]
    );
}

#[test]
fn same_stem_in_different_directories_does_not_pair() {
    let paths = vec![p("/a/lib_1.fq"), p("/b/lib_2.fq")];
    assert!(matches!(group_inputs(&paths)[0], InputGroup::Single(_)));
    assert_eq!(group_inputs(&paths).len(), 2);
}

#[test]
fn several_libraries_pair_independently() {
    let paths = vec![
        p("/in/a_1.fq.gz"),
        p("/in/b_1.fq.gz"),
        p("/in/a_2.fq.gz"),
        p("/in/b_2.fq.gz"),
    ];
    let groups = group_inputs(&paths);
    assert_eq!(groups.len(), 2);
    assert_eq!(
        groups[0],
        InputGroup::Paired(p("/in/a_1.fq.gz"), p("/in/a_2.fq.gz"))
    );
    assert_eq!(
        groups[1],
        InputGroup::Paired(p("/in/b_1.fq.gz"), p("/in/b_2.fq.gz"))
    );
}

#[test]
fn a_digit_that_is_not_a_mate_marker_is_left_alone() {
    // No underscore before the digit, so this is not a mate suffix.
    assert_eq!(split_mate("sample2.fq"), None);
    assert_eq!(split_mate("lib_L001.fq"), None);
}

/// The staged head must reach the sink too — losing it would silently drop
/// the first quarter-million reads of every library.
#[test]
fn short_input_still_drains_through_finish() {
    let insert = b"ACGTTGCAAGGTCCATTGACGATCGGATCCAAGGTTCCAAGGTTCCAAGGTTCCAA";
    let r2 = reverse_complement(insert);
    let q = vec![b'I'; insert.len()];

    let mut ing = Ingest::new(TrimOptions::default());
    let mut kept: Vec<Vec<u8>> = Vec::new();
    {
        let mut sink = |b: &[u8]| kept.push(b.to_vec());
        for _ in 0..10 {
            ing.push_pair(insert, &q, &r2, &q, &mut sink);
        }
        // Nothing has been emitted yet: it is all staged awaiting detection.
        ing.finish(&mut sink);
    }
    assert_eq!(kept.len(), 20, "10 pairs, both mates surviving");
    assert_eq!(ing.stats().reads_in, 20);
}

#[test]
fn finish_is_idempotent() {
    let mut ing = Ingest::new(TrimOptions::default());
    let mut kept = 0usize;
    {
        let mut sink = |_: &[u8]| kept += 1;
        let q = vec![b'I'; 8];
        ing.push_single(b"ACGTACGT", &q, &mut sink);
        ing.finish(&mut sink);
        ing.finish(&mut sink);
    }
    // The read is below the length filter, so nothing is emitted -- the
    // point is that the second finish does not replay the first.
    assert_eq!(ing.stats().reads_in, 1);
    assert_eq!(kept, 0);
}

/// Adapter-bearing reads must collapse once trimmed, which is the whole
/// reason trimming has to happen before hashing.
#[test]
fn trimming_makes_two_reads_identical() {
    let insert = b"ACGTTGCAAGGTCCATTGACGATCGGATCCAAGGTTCCAAGGTTCCAAGGTTCCAA";
    let mut a = insert.to_vec();
    a.extend_from_slice(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
    let mut b = insert.to_vec();
    b.extend_from_slice(b"AGATCGGAAGAGCAC"); // same adapter, shorter remnant
    let qa = vec![b'I'; a.len()];
    let qb = vec![b'I'; b.len()];

    let opts = TrimOptions {
        adapter_r1: Some(PreparedAdapter::new(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA")),
        ..Default::default()
    };
    let mut ing = Ingest::new(opts);
    let mut kept: Vec<Vec<u8>> = Vec::new();
    {
        let mut sink = |x: &[u8]| kept.push(x.to_vec());
        ing.push_single(&a, &qa, &mut sink);
        ing.push_single(&b, &qb, &mut sink);
        ing.finish(&mut sink);
    }
    assert_eq!(kept.len(), 2);
    assert_eq!(kept[0], kept[1], "untrimmed these would be distinct");
    assert_eq!(kept[0], insert.to_vec());
}
