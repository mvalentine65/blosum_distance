//! Refined splice sites of every exon in the chains. Each junction A->B is
//! aligned once; its path gives A's donor and B's acceptor. An exon takes its
//! acceptor from the junction before it and its donor from the junction after
//! it; gene ends keep their input coordinates.

use super::chain::{Chain, ChainExon, SRC_CUT, SRC_INPUT, SRC_NAME};
use super::junction::{JxStatus, Junction};
use super::sites::{base, revcomp, splice_scores};
use super::splice::MIN_INTRON;
use std::collections::HashMap;
use std::fmt::Write as _;

#[derive(Clone, Debug)]
pub struct JxSites {
    pub ok: bool,
    pub don_g: i64, // A's last base before the donor (genomic)
    pub acc_g: i64, // B's first base after the acceptor (genomic)
    pub don_site: [u8; 2],
    pub acc_site: [u8; 2],
    pub phase: i32,
    pub n_inner: i64,
    /// the path runs from A into B as one exon, with no frameshift or stop
    pub joined: bool,
}

impl Default for JxSites {
    fn default() -> Self {
        JxSites { ok: false, don_g: 0, acc_g: 0, don_site: *b"--", acc_site: *b"--", phase: -1, n_inner: 0, joined: false }
    }
}

/// anchors: sites only from a path through both anchors
pub fn junction_sites(jx: &Junction, wlo: i64, anchors: bool) -> JxSites {
    let mut js = JxSites::default();
    let n = jx.res.exons.len();
    // a path that skips either anchor says nothing about its sites
    if jx.status != JxStatus::Ok || n == 0 || (anchors && !(jx.res.exons[0].lo < jx.axe && jx.res.exons[n - 1].hi >= jx.bxs)) { return js; }
    js.joined = n == 1 && jx.res.dis.is_empty();
    if n < 2 { return js; }
    let (a, b) = (&jx.res.exons[0], &jx.res.exons[n - 1]);
    let dna = jx.dna();
    js.don_g = wlo - 1 + jx.pos(a.hi);
    js.acc_g = wlo - 1 + jx.pos(b.lo);
    if a.don >= 0 && a.don + 1 < jx.d { js.don_site = [dna[a.don as usize], dna[a.don as usize + 1]]; }
    if b.acc >= 1 { js.acc_site = [dna[b.acc as usize - 1], dna[b.acc as usize]]; }
    js.phase = b.phase;
    js.n_inner = n as i64 - 2;
    js.ok = true;
    js
}

/// codons a joined gap may hold beyond the model nodes it skips
const JOIN_EXTRA: f64 = 5.0;
/// a gap shorter than this holds no real intron
const JOIN_NO_INTRON: i64 = 50;
/// nt an intron may reach past a joined gap into either exon
const JOIN_SLOP: i64 = 3;

/// A GT..AG intron of at least MIN_INTRON nt fits in scaffold positions lo..=hi
/// (1-based), read on the gene strand.
fn intron_fits(sc: &[u8], strand: u8, lo: i64, hi: i64) -> bool {
    let (lo, hi) = ((lo - JOIN_SLOP).max(1) as usize, (hi + JOIN_SLOP).min(sc.len() as i64) as usize);
    if hi < lo + MIN_INTRON - 1 { return false; }
    let mut s: Vec<u8> = sc[lo - 1..hi].iter().map(|x| x.to_ascii_uppercase()).collect();
    if strand == b'-' { s = revcomp(&s); }
    let gt = (0..s.len() - 1).find(|&i| s[i] == b'G' && s[i + 1] == b'T');
    let ag = (1..s.len()).rev().find(|&j| s[j - 1] == b'A' && s[j] == b'G');
    matches!((gt, ag), (Some(i), Some(j)) if j + 1 >= i + MIN_INTRON)
}

/// donor score that marks a further exon past a chain's last one
const STOP_DONOR: f64 = 6.0;
/// nt inside the last exon where that donor may already lie
const STOP_DONOR_IN: i64 = 15;

/// A GT/GC donor scoring STOP_DONOR or more lies between STOP_DONOR_IN nt before the
/// exon end ce and the stop codon starting at p (1-based, read on the gene strand).
fn donor_ahead(sq: &[u8], plus: bool, ce: i64, p: i64) -> bool {
    let l = sq.len() as i64;
    // three exon bases before the first donor tried, six past the last
    let (lo, hi) = if plus { ((ce - STOP_DONOR_IN - 2).max(1), (p + 8).min(l)) } else { ((p - 8).max(1), (ce + STOP_DONOR_IN + 2).min(l)) };
    let s: Vec<u8> = if plus { sq[(lo - 1) as usize..hi as usize].to_vec() } else { revcomp(&sq[(lo - 1) as usize..hi as usize]) };
    let t: Vec<u8> = s.iter().map(|&x| base(x)).collect();
    let (mut ss5, mut ss3) = (vec![0f64; t.len()], vec![0f64; t.len()]);
    splice_scores(&t, &mut ss5, &mut ss3);
    // first base past the exon, and the stop codon's first base, in t
    let (off, end) = if plus { (ce - lo + 1, p - lo) } else { (hi - ce + 1, hi - p) };
    ((off - STOP_DONOR_IN).max(0)..=end).any(|k| {
        let k = k as usize;
        k + 1 < t.len() && t[k] == 2 && (t[k + 1] == 3 || t[k + 1] == 1) && ss5[k] >= STOP_DONOR
    })
}

#[derive(Clone, Debug)]
pub struct RefExon {
    pub acc_g: i64,
    pub don_g: i64,
    pub acc_site: [u8; 2],
    pub don_site: [u8; 2],
    pub phase_in: i32,
    pub n_inner_before: i64,
    /// gene ends: coding start (first base of the start codon) / coding end
    /// (last base of the stop codon), genomic; 0 when not extended
    pub start_g: i64,
    pub stop_g: i64,
}

impl Default for RefExon {
    fn default() -> Self {
        RefExon { acc_g: 0, don_g: 0, acc_site: *b"--", don_site: *b"--", phase_in: -1, n_inner_before: 0, start_g: 0, stop_g: 0 }
    }
}

pub struct Refine {
    pub ex: Vec<Vec<RefExon>>,
    /// join[i][j]: exon j of chain i and exon j + 1 are one exon
    pub join: Vec<Vec<bool>>,
}

impl Refine {
    pub fn new(cs: &[Chain]) -> Refine {
        Refine {
            ex: cs.iter().map(|c| vec![RefExon::default(); c.ex.len().max(1)]).collect(),
            join: cs.iter().map(|c| vec![false; c.ex.len().max(1)]).collect(),
        }
    }
    pub fn add(&mut self, ci: usize, ia: usize, ib: usize, js: &JxSites) {
        if js.joined && ib == ia + 1 { self.join[ci][ia] = true; }
        if !js.ok { return; }
        let a = &mut self.ex[ci][ia];
        if js.don_g != 0 {
            a.don_g = js.don_g;
            a.don_site = js.don_site;
        }
        let b = &mut self.ex[ci][ib];
        if js.acc_g != 0 {
            b.acc_g = js.acc_g;
            b.acc_site = js.acc_site;
        }
        b.phase_in = js.phase;
        b.n_inner_before = js.n_inner;
    }
    /// Refined start/end of exon j of chain i.
    pub fn span(&self, c: &Chain, i: usize, j: usize) -> (i64, i64) {
        let e = &c.ex[j];
        let r = &self.ex[i][j];
        let (mut ns, mut ne) = (e.start, e.end);
        if c.strand == b'+' {
            if r.acc_g != 0 { ns = r.acc_g; }
            if r.don_g != 0 { ne = r.don_g; }
            if r.start_g != 0 { ns = r.start_g; }
            if r.stop_g != 0 { ne = r.stop_g; }
        } else {
            if r.acc_g != 0 { ne = r.acc_g; }
            if r.don_g != 0 { ns = r.don_g; }
            if r.start_g != 0 { ne = r.start_g; }
            if r.stop_g != 0 { ns = r.stop_g; }
        }
        (ns, ne)
    }
    /// Gene ends: a chain's first exon, when it starts within START_NODES of
    /// the model's first node, extends to the farthest in-frame ATG upstream
    /// before a stop (within START_NT); its last exon, when it ends within
    /// STOP_NODES of the last node, extends through the first in-frame stop
    /// within STOP_NT. Alignments fade a few codons short of both. A last exon
    /// further from the model's end still reads on to a stop within STOP_FAR_NT
    /// unless a donor lies on the way (donor_ahead): a further exon.
    pub fn extend_ends(&mut self, cs: &[Chain], genome: &HashMap<String, Vec<u8>>, m_of: impl Fn(&Chain) -> Option<usize>) {
        const START_NODES: i64 = 30;
        const START_NT: i64 = 300;
        const STOP_NODES: i64 = 20;
        const STOP_NT: i64 = 90;
        const STOP_FAR_NT: i64 = 1500;
        for (i, c) in cs.iter().enumerate() {
            let (Some(m), Some(sq)) = (m_of(c), genome.get(&c.scaffold)) else { continue };
            let n = c.ex.len();
            if n == 0 { continue; }
            let plus = c.strand == b'+';
            // codon whose coding-first base is p (1-based)
            let codon = |p: i64| -> Option<[u8; 3]> {
                let (lo, hi) = if plus { (p, p + 2) } else { (p - 2, p) };
                if lo < 1 || hi > sq.len() as i64 { return None; }
                let t = &sq[(lo - 1) as usize..hi as usize];
                Some(if plus { [t[0], t[1], t[2]] } else { [comp(t[2]), comp(t[1]), comp(t[0])] })
            };
            let stop = |x: &[u8; 3]| matches!(x, b"TAA" | b"TAG" | b"TGA");
            if c.ex[0].k1 <= START_NODES {
                let (a, b) = self.span(c, i, 0);
                let cs0 = if plus { a } else { b };
                let step = if plus { -3 } else { 3 };
                let (mut p, mut atg) = (cs0, None);
                while (p - cs0).abs() <= START_NT {
                    let Some(x) = codon(p) else { break };
                    if stop(&x) { break; }
                    if &x == b"ATG" { atg = Some(p); }
                    p += step;
                }
                if let Some(p) = atg { if p != cs0 { self.ex[i][0].start_g = p; } }
            }
            let near = m as i64 - c.ex[n - 1].k2 <= STOP_NODES;
            let (a, b) = self.span(c, i, n - 1);
            let ce = if plus { b } else { a };
            let step = if plus { 3 } else { -3 };
            let mut p = ce + step / 3;
            while (p - ce).abs() <= STOP_FAR_NT {
                let Some(x) = codon(p) else { break };
                if stop(&x) {
                    if (near && (p - ce).abs() <= STOP_NT) || !donor_ahead(sq, plus, ce, p) {
                        self.ex[i][n - 1].stop_g = if plus { p + 2 } else { p - 2 };
                    }
                    break;
                }
                p += step;
            }
        }
    }

    /// gene scaffold strand idx start end new_start new_end acc_site don_site phase_in inner_before refined source
    /// Merge exons whose junction path is one exon: the merged exon keeps the
    /// first one's acceptor and the last one's donor. Rows of exonfill.joins.tsv
    /// (gene, merged exon, members start-end:source in gene order).
    pub fn apply_joins(&mut self, cs: &mut [Chain], genome: &HashMap<String, Vec<u8>>) -> String {
        let mut rows = String::new();
        for (i, c) in cs.iter_mut().enumerate() {
            if !self.join[i].iter().any(|&x| x) { continue; }
            let n = c.ex.len();
            let (mut ex, mut rf) = (Vec::with_capacity(n), Vec::with_capacity(n));
            let mut j = 0;
            // a path through a real intron can also come back as one exon; a split
            // exon's gap holds at most a few codons beyond the nodes it skips,
            // unless it is too short for an intron or no GT..AG intron fits in it
            let (sc, strand) = (genome.get(&c.scaffold), c.strand);
            let fits = |a: &ChainExon, b: &ChainExon| {
                let (lo, hi) = (a.end.min(b.end) + 1, a.start.max(b.start) - 1);
                (hi - lo + 1) as f64 / 3.0 - (b.k1 - a.k2 - 1) as f64 <= JOIN_EXTRA
                    || hi - lo + 1 < JOIN_NO_INTRON
                    || sc.is_some_and(|s| !intron_fits(s, strand, lo, hi))
            };
            while j < n {
                let mut k = j;
                // the two pieces of a cut exon stay apart
                while k + 1 < n && self.join[i][k] && c.ex[k].src != SRC_CUT && c.ex[k + 1].src != SRC_CUT && fits(&c.ex[k], &c.ex[k + 1]) { k += 1; }
                if k == j {
                    ex.push(c.ex[j]);
                    rf.push(self.ex[i][j].clone());
                } else {
                    let m = &c.ex[j..=k];
                    let src = if m.iter().any(|e| e.src == SRC_INPUT) { SRC_INPUT } else { m[0].src };
                    ex.push(ChainExon { start: m.iter().map(|e| e.start).min().unwrap(), end: m.iter().map(|e| e.end).max().unwrap(),
                                        k1: m[0].k1, k2: m[m.len() - 1].k2, src });
                    let mut r = self.ex[i][j].clone();
                    let last = &self.ex[i][k];
                    r.don_g = last.don_g;
                    r.don_site = last.don_site;
                    r.stop_g = last.stop_g;
                    rf.push(r);
                    let members: Vec<String> = m.iter().map(|e| format!("{}-{}:{}", e.start, e.end, SRC_NAME[e.src as usize])).collect();
                    let _ = writeln!(rows, "{}\t{}\t{}", c.gene, ex.len() - 1, members.join(";"));
                }
                j = k + 1;
            }
            self.join[i] = vec![false; ex.len().max(1)];
            c.ex = ex;
            self.ex[i] = rf;
        }
        rows
    }
    pub fn write(&self, cs: &[Chain]) -> String {
        let mut s = String::new();
        for (i, c) in cs.iter().enumerate() {
            for (j, e) in c.ex.iter().enumerate() {
                let r = &self.ex[i][j];
                let (ns, ne) = self.span(c, i, j);
                let what = match (r.acc_g != 0, r.don_g != 0) { (true, true) => "both", (true, false) => "acc", (false, true) => "don", _ => "none" };
                let _ = writeln!(s, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", c.gene, c.scaffold, c.strand as char, j,
                                 e.start, e.end, ns, ne, String::from_utf8_lossy(&r.acc_site), String::from_utf8_lossy(&r.don_site),
                                 r.phase_in, r.n_inner_before, what, SRC_NAME[e.src as usize]);
            }
        }
        s
    }
}

fn comp(b: u8) -> u8 {
    match b { b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', _ => b'N' }
}
