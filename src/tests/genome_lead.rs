//! Unit tests for [`genome_splice::lead`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;
use crate::genome_splice::sites::{acceptor_default, STANDARD};

const UPSTREAM: &str = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
const INTRON: &str = "CTCTCTCCTTCCTCTCTCTCCCTTCCTCTCTCTCCTTCTCTCTCTCC";
const PPT_AG: &str = "TTTTTTTTTTTTTCAG";
const DONOR: &str = "GTAAGT";

/// upstream + lead + donor + intron + acceptor + first exon; the lead exon's span and
/// the first exon's 5' base, 1-based, on the + strand
fn locus(lead: &str, first: &str) -> (Vec<u8>, (i64, i64), i64) {
    let s = format!("{UPSTREAM}{lead}{DONOR}{INTRON}{PPT_AG}{first}");
    let a = UPSTREAM.len() as i64 + 1;
    let b = a + lead.len() as i64 - 1;
    let five = (UPSTREAM.len() + lead.len() + DONOR.len() + INTRON.len() + PPT_AG.len()) as i64 + 1;
    (s.into_bytes(), (a, b), five)
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| match b { b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', x => x }).collect()
}

#[test]
fn test_lead_exon_in_frame() {
    let (sq, (a, b), five) = locus("ATGGCCGCCGCCCAG", "GCAGCAGCAGCAGCATAA");
    let ld = find_lead(STANDARD, &sq, true, five, 0, 0, &acceptor_default()).expect("a lead exon");
    assert_eq!((ld.start, ld.end, ld.atg), (a, b, a));
    assert_eq!(ld.intron, (DONOR.len() + INTRON.len() + PPT_AG.len()) as i64);
}

#[test]
fn test_lead_exon_completes_a_split_codon() {
    // one lead base and the first exon's two leading bases make GTC
    let (sq, (a, b), five) = locus("ATGGCAG", "TCGCAGCAGCAGCATAA");
    let ld = find_lead(STANDARD, &sq, true, five, 2, 0, &acceptor_default()).expect("a lead exon");
    assert_eq!((ld.start, ld.end), (a, b));
}

#[test]
fn test_no_lead_exon_when_the_split_codon_is_a_stop() {
    // T + AA reads TAA across the intron
    let (sq, _, five) = locus("ATGGCAT", "AAGCAGCAGCAGCATAA");
    assert_eq!(find_lead(STANDARD, &sq, true, five, 2, 0, &acceptor_default()), None);
}

#[test]
fn test_no_lead_exon_before_a_start_codon() {
    let (sq, _, five) = locus("ATGGCCGCCGCCCAG", "ATGCAGCAGCAGCATAA");
    assert_eq!(find_lead(STANDARD, &sq, true, five, 0, 0, &acceptor_default()), None);
}

#[test]
fn test_no_lead_exon_past_another_exon() {
    let (sq, (_, b), five) = locus("ATGGCCGCCGCCCAG", "GCAGCAGCAGCAGCATAA");
    // another row's exon ends inside the intron
    assert_eq!(find_lead(STANDARD, &sq, true, five, 0, b + 10, &acceptor_default()), None);
}

#[test]
fn test_no_lead_exon_without_an_acceptor() {
    let (mut sq, _, five) = locus("ATGGCCGCCGCCCAG", "GCAGCAGCAGCAGCATAA");
    sq[(five - 2) as usize] = b'C';
    assert_eq!(find_lead(STANDARD, &sq, true, five, 0, 0, &acceptor_default()), None);
}

#[test]
fn test_lead_exon_on_the_minus_strand() {
    let (fwd, (a, b), five) = locus("ATGGCCGCCGCCCAG", "GCAGCAGCAGCAGCATAA");
    let sq = revcomp(&fwd);
    let l = sq.len() as i64;
    let flip = |p: i64| l + 1 - p;
    let ld = find_lead(STANDARD, &sq, false, flip(five), 0, l + 1, &acceptor_default()).expect("a lead exon");
    assert_eq!((ld.start, ld.end, ld.atg), (flip(b), flip(a), flip(a)));
}

#[test]
fn test_the_better_donor_wins() {
    // a nearer lead exon behind a weaker donor (GTGAGC), a farther one behind GTAAGT
    let near = "ATGGCCCAGGTGAGC";
    let s = format!("{UPSTREAM}ATGGCCGCCCAG{DONOR}{INTRON}{near}{INTRON}{PPT_AG}GCAGCAGCAGCATAA");
    let sq = s.into_bytes();
    let five = (s_len_before(&[UPSTREAM, "ATGGCCGCCCAG", DONOR, INTRON, near, INTRON, PPT_AG])) as i64 + 1;
    let acc = acceptor_default();
    let ld = find_lead(STANDARD, &sq, true, five, 0, 0, &acc).expect("a lead exon");
    let far_end = (UPSTREAM.len() + "ATGGCCGCCCAG".len()) as i64;
    let near_end = (s_len_before(&[UPSTREAM, "ATGGCCGCCCAG", DONOR, INTRON]) + "ATGGCCCAG".len()) as i64;
    let score = |e: i64| site_scores(&sq, true, e, five, &acc).0;
    assert!(score(far_end) > score(near_end));
    assert_eq!(ld.end, far_end);
}

fn s_len_before(parts: &[&str]) -> usize {
    parts.iter().map(|p| p.len()).sum()
}

#[test]
fn test_a_lead_exon_fits_the_nodes_the_row_lacks() {
    assert!(crate::genome_splice::gate::lead_fits(4, 3));
    assert!(crate::genome_splice::gate::lead_fits(3, 5));
    assert!(!crate::genome_splice::gate::lead_fits(9, 5));
    // a fragment short of a long exon 1 takes no few-codon lead exon
    assert!(!crate::genome_splice::gate::lead_fits(4, 210));
}
