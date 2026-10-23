//! Lead exons: a row whose first exon opens on an acceptor, with no start codon
//! of its own, takes the exon 1 a splice-signal search finds upstream. The lead
//! exon starts at the farthest in-frame ATG before a stop, ends on the best
//! donor within LEAD_WINDOW, completes the codon split across its intron, and
//! is about as long as the model nodes the row lacks (gate::lead_fits).

use super::chain::{left_bound, right_bound, Chain, ChainExon, SRC_LEAD};
use super::gate::{keep_lead, lead_fits, lead_opens, LEAD_MAX_NT, LEAD_MIN_INTRON, LEAD_WINDOW};
use super::refine::Refine;
use super::sites::{site_scores, translate, AccTable, Code};
use std::collections::HashMap;
use std::fmt::Write as _;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Lead {
    /// genomic span of the lead exon
    pub start: i64,
    pub end: i64,
    /// first base of its start codon, coding orientation
    pub atg: i64,
    pub donor: f64,
    pub acceptor: f64,
    pub intron: i64,
}

fn comp(b: u8) -> u8 {
    match b { b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', _ => b'N' }
}

/// The lead exon of a first exon whose 5' base is `five`, its first full codon
/// `f` bases in. `bound`: the nearest base of another exon upstream (0 or
/// len + 1 when none); neither the lead exon nor its intron reaches it.
pub fn find_lead(code: Code, sq: &[u8], plus: bool, five: i64, f: i64, bound: i64, acc: &AccTable) -> Option<Lead> {
    let l = sq.len() as i64;
    // genomic position u bases upstream of five, and the base there on the gene strand
    let pos = |u: i64| if plus { five - u } else { five + u };
    let up = |u: i64| -> u8 {
        let p = pos(u);
        if p < 1 || p > l { return b'N'; }
        let b = sq[(p - 1) as usize].to_ascii_uppercase();
        if plus { b } else { comp(b) }
    };
    let down = |k: i64| -> u8 { up(-k) };
    if up(2) != b'A' || up(1) != b'G' { return None; }
    let acceptor = site_scores(sq, plus, five, five, acc).1;
    if !lead_opens(acceptor) { return None; }
    if down(0) == b'A' && down(1) == b'T' && down(2) == b'G' { return None; }
    let inside = |u: i64| { let p = pos(u); p >= 1 && p <= l && if plus { p > bound } else { p < bound } };
    // bases of the split codon the lead exon holds, then the first exon's
    let r = (3 - f.rem_euclid(3)) % 3;
    let tail: Vec<u8> = (0..(3 - r) % 3).map(down).collect();
    let mut best: Option<Lead> = None;
    // ue: the lead exon's last base, u bases upstream; its intron is ue - 1 long
    for ue in (LEAD_MIN_INTRON + 1)..=LEAD_WINDOW {
        if !inside(ue) { break; }
        if up(ue - 1) != b'G' || !matches!(up(ue - 2), b'T' | b'C') { continue; }
        if r > 0 {
            let mut c = [0u8; 3];
            for (i, x) in ((ue..ue + r).rev().map(up).chain(tail.iter().copied())).enumerate() { c[i] = x; }
            if translate(code, &c) == b'*' { continue; }
        }
        // codons upstream of the split one: the farthest ATG before a stop
        let (mut cs, mut atg) = (ue + r + 2, None);
        while cs - ue < LEAD_MAX_NT && inside(cs) {
            let c = [up(cs), up(cs - 1), up(cs - 2)];
            if translate(code, &c) == b'*' { break; }
            if &c == b"ATG" { atg = Some(cs); }
            cs += 3;
        }
        let Some(us) = atg else { continue };
        let donor = site_scores(sq, plus, pos(ue), five, acc).0;
        if !keep_lead(donor) { continue; }
        // the best donor; the nearer of two that tie
        if best.is_some_and(|b| donor <= b.donor) { continue; }
        let (a, b) = (pos(us), pos(ue));
        best = Some(Lead { start: a.min(b), end: a.max(b), atg: a, donor, acceptor, intron: ue - 1 });
    }
    best
}

/// Adds a lead exon to each row that earns one; rows of the log
/// (gene scaffold strand start end donor acceptor intron).
pub fn add_lead_exons(cs: &mut [Chain], rf: &Refine, genome: &HashMap<String, Vec<u8>>, code: Code, acc: &AccTable) -> String {
    let mut iv: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    for c in cs.iter() {
        let v = iv.entry(c.scaffold.clone()).or_default();
        for e in &c.ex { v.push((e.start, e.end)); }
    }
    for v in iv.values_mut() { v.sort_by_key(|x| x.0); }
    let mut log = String::new();
    for (i, c) in cs.iter_mut().enumerate() {
        if c.passive || c.flank_only || c.imx || c.ex.len() < 2 || c.ex[0].src == SRC_LEAD || c.ex[0].codon == 0 { continue; }
        if rf.ex[i][0].start_g != 0 { continue; }
        let Some(sq) = genome.get(&c.scaffold) else { continue };
        let plus = c.strand == b'+';
        let (a, b) = rf.span(c, i, 0);
        let five = if plus { a } else { b };
        let f = if plus { c.ex[0].codon - five } else { five - c.ex[0].codon };
        let ivs = &iv[&c.scaffold];
        let bound = if plus { left_bound(ivs, a) } else { right_bound(ivs, b, sq.len() as i64 + 1) };
        let Some(ld) = find_lead(code, sq, plus, five, f, bound, acc) else { continue };
        if !lead_fits((ld.end - ld.start + 1) / 3, c.ex[0].k1 - c.klo) { continue; }
        c.ex.insert(0, ChainExon { start: ld.start, end: ld.end, k1: 0, k2: 0, src: SRC_LEAD, codon: ld.atg, bits: 0.0, run: false, owes: false });
        let _ = writeln!(log, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", c.gene, c.scaffold, c.strand as char, ld.start, ld.end, ld.donor, ld.acceptor, ld.intron);
    }
    log
}

#[cfg(test)]
#[path = "../tests/genome_lead.rs"]
mod tests;
