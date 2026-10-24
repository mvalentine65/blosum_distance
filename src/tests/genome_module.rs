//! Unit tests for [`genome_splice::module`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;
use crate::genome_splice::sites::code;

// straw (Drosophila Laccase2): the exon of isoforms A/F/E and its mutually exclusive copy in D/G
const AF: &[u8] = b"VMNPLDAQCNRQRNDAICVSQLKNALEIDRGILAEKPDVKIFLPFRFFVYRAEDLFQPNTYNRFL";
const DG: &[u8] = b"ILNPLDAICDRKRADAVCVSNLKNAKKVDKGVLVERPDVKIFLPFRFFVYEPKALFIPNTYNRFL";

/// One codon per residue, so a peptide reads back as itself in frame 0.
fn back(p: &[u8]) -> Vec<u8> {
    p.iter().flat_map(|&a| match a {
        b'A' => *b"GCT", b'C' => *b"TGT", b'D' => *b"GAT", b'E' => *b"GAA", b'F' => *b"TTT", b'G' => *b"GGT",
        b'H' => *b"CAT", b'I' => *b"ATT", b'K' => *b"AAA", b'L' => *b"CTT", b'M' => *b"ATG", b'N' => *b"AAT",
        b'P' => *b"CCT", b'Q' => *b"CAA", b'R' => *b"CGT", b'S' => *b"TCT", b'T' => *b"ACT", b'V' => *b"GTT",
        b'W' => *b"TGG", _ => *b"TAT",
    }).collect()
}

fn exon(start: i64, end: i64) -> ChainExon {
    ChainExon { start, end, k1: 297, k2: 345, src: 0, codon: 0, bits: 0.0, run: false, owes: false }
}

#[test]
fn test_local_finds_the_query_inside_a_longer_sequence() {
    let t: Vec<u8> = [&b"GGGGGGGGGG"[..], AF, &b"GGGGG"[..]].concat();
    let (score, first, last) = local(AF, &t);
    let own: i32 = AF.iter().map(|&a| b62(a, a)).sum();
    assert_eq!((score, first, last), (own, 10, 10 + AF.len() - 1));
}

#[test]
fn test_local_scores_the_paralog_below_the_query_itself() {
    let own = local(AF, AF).0;
    let other = local(AF, DG).0;
    assert!(other > own / 2 && other < own, "paralog {other} vs self {own}");
}

#[test]
fn test_seq_copy_finds_the_other_copy_and_skips_masked_ones() {
    // the A/F exon at 101.., 600 nt of spacer, then the D/G copy; lo..hi spans both
    let spacer = vec![b'C'; 600];
    let mut sc = vec![b'C'; 100];
    sc.extend(back(AF));
    sc.extend(&spacer);
    let dg_start = sc.len() as i64 + 1;
    sc.extend(back(DG));
    sc.extend(vec![b'C'; 100]);
    let x = exon(101, 100 + 3 * AF.len() as i64);
    let lo = 1;
    let hi = sc.len() as i64;
    let g = code(1).unwrap();
    // the sibling's own exon is masked (taken), so only its copy can be found
    let (e, margin) = seq_copy(g, b'+', &x, &sc, lo, hi, &[(x.start, x.end)], 0.75).expect("copy found");
    assert_eq!((e.start, e.end), (dg_start, dg_start + 3 * DG.len() as i64 - 1));
    assert_eq!((e.k1, e.k2, e.src), (297, 345, SRC_SIBFILL));
    assert!(margin >= 60.0);
    // with the copy masked too, nothing is left to find
    let both = [(x.start, x.end), (e.start, e.end)];
    assert!(seq_copy(g, b'+', &x, &sc, lo, hi, &both, 0.75).is_none());
}

#[test]
fn test_score_gate_drops_a_hit_scoring_under_the_ratio_of_its_exon() {
    // the exon scores 40 bits; a copy needs 0.3 of that
    assert!(!gate_keeps(11.0, 0.3, 40.0, 10.0));
    assert!(gate_keeps(12.5, 0.3, 40.0, 10.0));
    assert!(gate_keeps(1.0, 0.3, 3.0, 10.0)); // a near-zero exon says nothing about its copies
    assert!(gate_keeps(1.0, 0.0, 40.0, 10.0)); // ratio 0: gate off
}
