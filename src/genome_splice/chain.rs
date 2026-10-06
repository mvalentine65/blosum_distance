//! Exon chains in; genome and the gap and flank windows to fill, out.
//!
//! chains.tsv, tab separated, '#' lines ignored:
//!   G  gene  model  scaffold  strand  k_lo  k_hi  [flank_only|passive|imx]
//!   E  gene  start  end  k1  k2  [codon]
//! codon: genomic first base (coding orientation) of one codon of the hit; it
//! fixes the exon's reading frame. Without it the frame is found by score.
//! Per chain, exons in model order. A gap window lies between consecutive
//! exons missing min_gap..max_gap-1 nodes; a flank window lies beyond the
//! first/last exon when k_lo/k_hi leave min_gap..max_end-1 nodes uncovered, reaching up
//! to `flank` nt out and stopping at the nearest other exon on the scaffold.
//! Chains out carry each exon's source and the tags alt_of= and rebased=.

use flate2::read::MultiGzDecoder;
use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

pub const SRC_INPUT: u8 = 0;
pub const SRC_ORF: u8 = 1;
pub const SRC_SPLICE: u8 = 2;
pub const SRC_ALT: u8 = 3;
pub const SRC_TAIL: u8 = 4;
pub const SRC_CUT: u8 = 5;
pub const SRC_NAME: [&str; 6] = ["input", "orf", "splice", "alt", "tail", "cut"];

#[derive(Clone, Copy, Debug)]
pub struct ChainExon {
    pub start: i64,
    pub end: i64,
    pub k1: i64,
    pub k2: i64,
    pub src: u8,
    /// genomic first base (coding orientation) of one codon of the hit; 0: not given
    pub codon: i64,
}

#[derive(Clone, Debug)]
pub struct Chain {
    pub gene: String,
    pub model: String,
    pub scaffold: String,
    pub strand: u8,
    pub klo: i64,
    pub khi: i64,
    pub flank_only: bool,
    pub passive: bool,
    /// an IMX isoform row: gets no alternatives of its own (its base row does)
    pub imx: bool,
    /// an alternative isoform: parent chain and the start-end of the exon it replaces
    pub alt_of: Option<String>,
    /// rebased onto a recovered tail: the chain that owns it
    pub rebased: Option<String>,
    pub ex: Vec<ChainExon>,
}

pub fn sort_exons(ex: &mut [ChainExon]) {
    ex.sort_by(|a, b| (a.k1, a.start).cmp(&(b.k1, b.start)));
}

pub fn read_chains(path: &str) -> Result<Vec<Chain>, String> {
    let f = File::open(path).map_err(|e| format!("{path}: {e}"))?;
    let mut chains: Vec<Chain> = Vec::new();
    let mut idx: HashMap<String, usize> = HashMap::new();
    for line in BufReader::new(f).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.is_empty() || line.starts_with('#') { continue; }
        let f: Vec<&str> = line.split('\t').filter(|s| !s.is_empty()).collect();
        if f[0].starts_with('G') && f.len() >= 7 {
            if idx.contains_key(f[1]) { continue; }
            idx.insert(f[1].to_string(), chains.len());
            chains.push(Chain {
                gene: f[1].into(), model: f[2].into(), scaffold: f[3].into(), strand: f[4].as_bytes()[0],
                klo: f[5].parse().unwrap_or(0), khi: f[6].parse().unwrap_or(0),
                flank_only: f.len() >= 8 && f[7] == "flank_only", passive: f.len() >= 8 && f[7] == "passive",
                imx: f.len() >= 8 && f[7] == "imx", alt_of: None, rebased: None, ex: Vec::new(),
            });
        } else if f[0].starts_with('E') && f.len() >= 6 {
            if let Some(&i) = idx.get(f[1]) {
                chains[i].ex.push(ChainExon {
                    start: f[2].parse().unwrap_or(0), end: f[3].parse().unwrap_or(0),
                    k1: f[4].parse().unwrap_or(0), k2: f[5].parse().unwrap_or(0), src: SRC_INPUT,
                    codon: f.get(6).and_then(|x| x.parse().ok()).unwrap_or(0),
                });
            }
        }
    }
    for c in chains.iter_mut() { sort_exons(&mut c.ex); }
    Ok(chains)
}

/// Scaffolds some chain names, upper-cased. FASTA, plain or gzipped.
pub fn load_genome(path: &str, chains: &[Chain]) -> Result<HashMap<String, Vec<u8>>, String> {
    let need: HashSet<&str> = chains.iter().map(|c| c.scaffold.as_str()).collect();
    let f = File::open(path).map_err(|e| format!("{path}: {e}"))?;
    let mut rd: Box<dyn Read> = if path.ends_with(".gz") { Box::new(MultiGzDecoder::new(BufReader::with_capacity(1 << 20, f))) } else { Box::new(f) };
    let mut text = Vec::new();
    rd.read_to_end(&mut text).map_err(|e| format!("{path}: {e}"))?;
    let mut g: HashMap<String, Vec<u8>> = HashMap::new();
    let mut cur: Option<String> = None;
    let mut buf: Vec<u8> = Vec::new();
    for line in text.split(|&c| c == b'\n') {
        if line.first() == Some(&b'>') {
            if let Some(n) = cur.take() { g.insert(n, std::mem::take(&mut buf)); }
            let name = String::from_utf8_lossy(&line[1..]).split_whitespace().next().unwrap_or("").to_string();
            if need.contains(name.as_str()) && !g.contains_key(&name) { cur = Some(name); }
            buf.clear();
        } else if cur.is_some() {
            buf.extend(line.iter().filter(|c| !c.is_ascii_whitespace()).map(|c| c.to_ascii_uppercase()));
        }
    }
    if let Some(n) = cur { g.insert(n, buf); }
    Ok(g)
}

pub const GAP: u8 = 0;
pub const LEAD: u8 = 1;
pub const TRAIL: u8 = 2;
pub const KIND_NAME: [&str; 3] = ["gap", "lead", "trail"];

#[derive(Clone, Debug)]
pub struct Window {
    pub kind: u8,
    pub ci: usize,
    pub ia: usize,
    pub ib: Option<usize>,
    pub k1: i64,
    pub k2: i64,
    pub gs: i64,
    pub ge: i64,
    /// other chains with this same flank window: they get its exons too
    pub also: Vec<usize>,
    /// gap part left unsearched (read as N): beyond a sibling chain's alternative of either exon
    pub mask: Option<(i64, i64)>,
}

pub struct ChainOpts {
    pub min_gap: i64,
    pub max_gap: i64,
    /// most missing nodes past a chain end for a lead/trail window
    pub max_end: i64,
    pub flank: i64,
    /// chains sharing a flank window all get its exons (off: the first only, as the C version)
    pub share: bool,
    /// gap windows stop at a sibling chain's alternative of either exon
    pub siblings: bool,
}

fn left_bound(iv: &[(i64, i64)], pos: i64) -> i64 {
    let mut i = iv.len() as i64 - 1;
    while i >= 0 && iv[i as usize].0 >= pos { i -= 1; }
    while i >= 0 {
        if iv[i as usize].1 < pos { return iv[i as usize].1; }
        i -= 1;
    }
    0
}

fn right_bound(iv: &[(i64, i64)], pos: i64, l: i64) -> i64 {
    iv.iter().find(|x| x.0 > pos).map(|x| x.0).unwrap_or(l)
}

/// only: build just these (chain, GAP|LEAD|TRAIL) windows; GAP means all of that chain's gaps
pub fn chain_windows(cs: &[Chain], genome: &HashMap<String, Vec<u8>>, o: &ChainOpts, only: Option<&HashSet<(usize, u8)>>) -> Vec<Window> {
    let mut iv: HashMap<&str, Vec<(i64, i64)>> = HashMap::new();
    for c in cs {
        let v = iv.entry(c.scaffold.as_str()).or_default();
        for e in &c.ex { v.push((e.start, e.end)); }
    }
    for v in iv.values_mut() { v.sort_by_key(|x| x.0); }
    let mut w: Vec<Window> = Vec::new();
    for (i, c) in cs.iter().enumerate() {
        let Some(s) = genome.get(&c.scaffold) else { continue };
        if c.ex.is_empty() || c.passive { continue; }
        let l = s.len() as i64;
        let ivs = &iv[c.scaffold.as_str()];
        if !c.flank_only && only.is_none_or(|s| s.contains(&(i, GAP))) {
            for j in 1..c.ex.len() {
                let (a, b) = (&c.ex[j - 1], &c.ex[j]);
                let miss = b.k1 - a.k2 - 1;
                if miss < o.min_gap || miss >= o.max_gap { continue; }
                let gs = if a.start < b.start { a.end + 1 } else { b.end + 1 };
                let ge = if a.start < b.start { b.start - 1 } else { a.start - 1 };
                if ge < gs || gs < 1 || ge > l { continue; }
                let mask = if o.siblings { sibling_mask(cs, i, a, b, gs, ge, false) } else { None };
                w.push(Window { kind: GAP, ci: i, ia: j - 1, ib: Some(j), k1: a.k2 + 1, k2: b.k1 - 1, gs, ge, also: Vec::new(), mask });
            }
        }
        for side in 0..2 {
            let e = if side == 0 { &c.ex[0] } else { &c.ex[c.ex.len() - 1] };
            let miss = if side == 0 { e.k1 - c.klo } else { c.khi - e.k2 };
            if miss < o.min_gap || miss >= o.max_end { continue; }
            let left = (side == 0) == (c.strand == b'+');
            let (mut gs, mut ge) = if left {
                (left_bound(ivs, e.start).max(e.start - 1 - o.flank) + 1, e.start - 1)
            } else {
                (e.end + 1, right_bound(ivs, e.end, l).min(e.end + o.flank))
            };
            gs = gs.max(1);
            ge = ge.min(l);
            if ge < gs { continue; }
            let kind = if side == 0 { LEAD } else { TRAIL };
            if only.is_some_and(|s| !s.contains(&(i, kind))) { continue; }
            let dup = w.iter().position(|x| {
                x.kind == kind && x.gs == gs && x.ge == ge && cs[x.ci].strand == c.strand
                    && cs[x.ci].scaffold == c.scaffold && cs[x.ci].model == c.model
            });
            if let Some(d) = dup {
                if o.share { w[d].also.push(i); }
                continue;
            }
            let (k1, k2) = if side == 0 { (c.klo, e.k1 - 1) } else { (e.k2 + 1, c.khi) };
            w.push(Window { kind, ci: i, ia: if side == 0 { 0 } else { c.ex.len() - 1 }, ib: None, k1, k2, gs, ge, also: Vec::new(), mask: None });
        }
    }
    w
}

/// Part of the gap a..b (gs..ge) that belongs to other isoforms: the block of
/// module copies in it. A module copy is an exon of a sibling chain (same
/// model) lying in the gap that covers at least half of a's nodes or b's, or
/// has a twin (another such exon on the same nodes at another locus). In a
/// module each alternative's own neighbouring exons sit next to it, outside the
/// other alternatives' block; constitutive exons past the block stay searched.
/// alts_only: copies are only alternatives to a or b (half of both node spans), no twins.
pub fn sibling_mask(cs: &[Chain], i: usize, a: &ChainExon, b: &ChainExon, gs: i64, ge: i64, alts_only: bool) -> Option<(i64, i64)> {
    let c = &cs[i];
    let same_nodes = |x: &ChainExon, e: &ChainExon| {
        let ov = x.k2.min(e.k2) - x.k1.max(e.k1) + 1;
        2 * ov >= (x.k2 - x.k1).min(e.k2 - e.k1) + 1
    };
    let mut sib: Vec<ChainExon> = Vec::new();
    for (k, o) in cs.iter().enumerate() {
        if k == i || o.passive || o.model != c.model || o.scaffold != c.scaffold || o.strand != c.strand { continue; }
        for x in &o.ex {
            if x.start < gs || x.end > ge || c.ex.iter().any(|y| y.start <= x.end && x.start <= y.end) { continue; }
            if !sib.iter().any(|y| y.start == x.start && y.end == x.end) { sib.push(*x); }
        }
    }
    let alt = |x: &ChainExon, e: &ChainExon| {
        let ov = x.k2.min(e.k2) - x.k1.max(e.k1) + 1;
        2 * ov >= (x.k2 - x.k1).max(e.k2 - e.k1) + 1
    };
    let copy = |x: &ChainExon| if alts_only { alt(x, a) || alt(x, b) } else {
        same_nodes(x, a) || same_nodes(x, b) || sib.iter().any(|y| (y.end < x.start || x.end < y.start) && same_nodes(x, y))
    };
    let mut span: Option<(i64, i64)> = None;
    for x in sib.iter().filter(|x| copy(x)) {
        span = Some(match span { Some((lo, hi)) => (lo.min(x.start), hi.max(x.end)), None => (x.start, x.end) });
    }
    span
}

/// Block of other modules in gap gs..ge of chain i: exons there of sibling chains that share `share` and pass `test`.
pub fn module_block(cs: &[Chain], i: usize, share: &ChainExon, gs: i64, ge: i64, test: impl Fn(&ChainExon) -> bool) -> Option<(i64, i64)> {
    let c = &cs[i];
    let mut span: Option<(i64, i64)> = None;
    for (k, o) in cs.iter().enumerate() {
        if k == i || o.passive || o.model != c.model || o.scaffold != c.scaffold || o.strand != c.strand { continue; }
        if !o.ex.iter().any(|y| y.start == share.start && y.end == share.end) { continue; }
        for x in &o.ex {
            if x.start < gs || x.end > ge || !test(x) || c.ex.iter().any(|y| y.start <= x.end && x.start <= y.end) { continue; }
            span = Some(match span { Some((lo, hi)) => (lo.min(x.start), hi.max(x.end)), None => (x.start, x.end) });
        }
    }
    span
}

/// Every consecutive exon pair; pairs may share a few nodes, but not one inside the other.
pub fn chain_junctions(cs: &[Chain]) -> Vec<(usize, usize, usize)> {
    let mut j = Vec::new();
    for (i, c) in cs.iter().enumerate() {
        for k in 1..c.ex.len() {
            let (a, b) = (&c.ex[k - 1], &c.ex[k]);
            if b.k1 <= a.k1 || b.k2 <= a.k2 { continue; }
            j.push((i, k - 1, k));
        }
    }
    j
}

/// Add a recovered exon unless it overlaps an exon of the chain on the genome,
/// or over half the shorter node span; keeps model order.
pub fn add_exon(c: &mut Chain, e: &ChainExon) -> bool {
    for x in &c.ex {
        if e.start <= x.end && x.start <= e.end { return false; }
        let ov = e.k2.min(x.k2) - e.k1.max(x.k1) + 1;
        if 2 * ov > (e.k2 - e.k1).min(x.k2 - x.k1) + 1 { return false; }
    }
    c.ex.push(*e);
    sort_exons(&mut c.ex);
    true
}

pub fn write_chains(cs: &[Chain]) -> String {
    let mut s = String::from("# G gene model scaffold strand k_lo k_hi [flag] [alt_of=parent:start-end] [rebased=owner]\n# E gene start end k1 k2 source\n");
    for c in cs {
        let flag = if c.flank_only { "\tflank_only" } else if c.passive { "\tpassive" } else if c.imx { "\timx" } else { "" };
        let mut tags = String::new();
        if let Some(a) = &c.alt_of { let _ = write!(tags, "\talt_of={a}"); }
        if let Some(r) = &c.rebased { let _ = write!(tags, "\trebased={r}"); }
        let pad = if flag.is_empty() && !tags.is_empty() { "\t-" } else { "" };
        let _ = writeln!(s, "G\t{}\t{}\t{}\t{}\t{}\t{}{}{}{}", c.gene, c.model, c.scaffold, c.strand as char, c.klo, c.khi, flag, pad, tags);
        for e in &c.ex {
            let _ = writeln!(s, "E\t{}\t{}\t{}\t{}\t{}\t{}", c.gene, e.start, e.end, e.k1, e.k2, SRC_NAME[e.src as usize]);
        }
    }
    s
}
