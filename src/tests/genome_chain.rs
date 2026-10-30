//! Unit tests for [`genome_splice::chain`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;

fn ex(start: i64, end: i64, k1: i64, k2: i64) -> ChainExon {
    ChainExon { start, end, k1, k2, src: SRC_INPUT, codon: 0, bits: 0.0, run: false, owes: false }
}

fn row(gene: &str, ex: Vec<ChainExon>) -> Chain {
    Chain { gene: gene.into(), model: "m".into(), scaffold: "s".into(), strand: b'+', klo: 1, khi: 800, flank_only: false,
            passive: false, imx: false, alt_of: None, rebased: None, ex }
}

// a slot (nodes 613-708) whose sibling variant lies in the gap before row 0's variant
const A: (i64, i64, i64, i64) = (1000, 1159, 559, 612);
const VAR: (i64, i64, i64, i64) = (10000, 10287, 613, 708);
const NEXT: (i64, i64, i64, i64) = (22000, 22122, 709, 749);

fn mask_of(row0: Vec<ChainExon>) -> Option<(i64, i64)> {
    let e = |t: (i64, i64, i64, i64)| ex(t.0, t.1, t.2, t.3);
    let cs = vec![row("r#0", row0.clone()), row("r#1", vec![e(A), e(VAR), e(NEXT)])];
    sibling_mask(&cs, 0, &row0[0], &row0[1], A.1 + 1, row0[1].start - 1, true)
}

#[test]
fn test_sibling_mask_masks_the_slot_before_a_split_alternative() {
    // row 0's variant is split in two (618-661, 669-708): together they are the slot's alternative
    let r0 = vec![ex(A.0, A.1, A.2, A.3), ex(20000, 20131, 618, 661), ex(21000, 21119, 669, 708), ex(NEXT.0, NEXT.1, NEXT.2, NEXT.3)];
    assert_eq!(mask_of(r0), Some((VAR.0, VAR.1)));
}

#[test]
fn test_sibling_mask_leaves_the_slot_open_for_a_lone_partial_exon() {
    // only the first half: not an alternative of the full-length variant
    let r0 = vec![ex(A.0, A.1, A.2, A.3), ex(20000, 20131, 618, 661), ex(NEXT.0, NEXT.1, NEXT.2, NEXT.3)];
    assert_eq!(mask_of(r0), None);
}

#[test]
fn test_sibling_mask_masks_the_slot_before_a_full_length_alternative() {
    let r0 = vec![ex(A.0, A.1, A.2, A.3), ex(20000, 20287, 613, 708), ex(NEXT.0, NEXT.1, NEXT.2, NEXT.3)];
    assert_eq!(mask_of(r0), Some((VAR.0, VAR.1)));
}
