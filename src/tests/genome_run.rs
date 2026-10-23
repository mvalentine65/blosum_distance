//! Unit tests for [`genome_splice::run`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;

fn ex(start: i64, end: i64, k1: i64, k2: i64, src: u8) -> ChainExon {
    ChainExon { start, end, k1, k2, src, codon: 0, bits: 0.0, run: false, owes: false }
}

fn row(gene: &str, ex: Vec<ChainExon>) -> Chain {
    Chain { gene: gene.into(), model: "AA1".into(), scaffold: "NT_033778.4".into(), strand: b'+', klo: 1, khi: 525,
            flank_only: false, passive: false, imx: false, alt_of: None, rebased: None, ex }
}

// Drosophila straw: exon 43037, then X1 (A/F's block) or X2 (D's and G's), then the ends
fn straw() -> Vec<Chain> {
    let e37 = ex(5439570, 5439754, 240, 293, SRC_INPUT);
    let x1 = ex(5440610, 5440805, 297, 345, SRC_ORF);
    let x2 = ex(5453281, 5453476, 297, 345, SRC_SIBFILL);
    vec![
        row("AA1#2", vec![e37, x1, ex(5440858, 5441070, 353, 433, SRC_INPUT), ex(5446535, 5446765, 438, 499, SRC_INPUT)]),
        row("AA1#2_1", vec![e37, x2, ex(5453559, 5454151, 354, 521, SRC_INPUT)]),
        row("AA1#2_2", vec![e37, x1, ex(5454421, 5455013, 357, 521, SRC_INPUT)]),
    ]
}

#[test]
fn test_an_end_takes_the_copy_beside_it() {
    let mut cs = straw();
    let log = nearest_copies(&mut cs, 0.75);
    // 43049's row gives up X1 (13.6 kb away) for X2 (945 nt away): G, not E
    assert_eq!((cs[2].ex[1].start, cs[2].ex[1].end, cs[2].ex[1].src), (5453281, 5453476, SRC_SIBFILL));
    assert_eq!(log.lines().count(), 1);
    // A/F keeps X1 (X2 lies outside its gap) and D keeps X2 (nearer 43048 than X1)
    assert_eq!(cs[0].ex[1].start, 5440610);
    assert_eq!(cs[1].ex[1].start, 5453281);
}

#[test]
fn test_a_start_keeps_its_local_copy() {
    // P450 alternative start: start B with its own exon 2, the sibling's exon 2 further on
    let body = ex(22639189, 22639279, 92, 125, SRC_INPUT);
    let mut cs = vec![
        row("clan4#30", vec![ex(22637013, 22637169, 5, 53, SRC_INPUT), ex(22638408, 22638517, 57, 93, SRC_SPLICE), body]),
        row("clan4#30_1", vec![ex(22635270, 22635426, 8, 52, SRC_INPUT), ex(22636271, 22636380, 57, 93, SRC_SIBFILL), body]),
    ];
    assert_eq!(nearest_copies(&mut cs, 0.75), "");
    assert_eq!(cs[1].ex[1].start, 22636271);
}

// A row's recovered exon e (1400-1600) whose own flank is b; the sibling shares a and offers y
fn swap_case(y: ChainExon) -> Vec<Chain> {
    let a = ex(1000, 1100, 1, 30, SRC_INPUT);
    vec![
        row("g#1", vec![a, y, ex(3000, 3100, 91, 120, SRC_INPUT)]),
        row("g#1_1", vec![a, ex(1400, 1600, 31, 90, SRC_ORF), ex(1800, 1900, 91, 120, SRC_INPUT)]),
    ]
}

#[test]
fn test_a_copy_overlapping_the_exon_is_not_swapped_in() {
    let mut cs = swap_case(ex(1500, 1700, 31, 90, SRC_ORF));
    assert_eq!(nearest_copies(&mut cs, 0.75), "");
    assert_eq!(cs[1].ex[1].start, 1400);
}

#[test]
fn test_a_copy_of_unlike_length_is_not_swapped_in() {
    // 101 nt for 201: under 0.75
    let mut cs = swap_case(ex(1650, 1750, 31, 90, SRC_ORF));
    assert_eq!(nearest_copies(&mut cs, 0.75), "");
    // 151 nt for 201: just over 0.75, and nearer b
    let mut cs = swap_case(ex(1610, 1760, 31, 90, SRC_ORF));
    assert_eq!(nearest_copies(&mut cs, 0.75).lines().count(), 1);
    assert_eq!((cs[1].ex[1].start, cs[1].ex[1].src), (1610, SRC_SIBFILL));
}
