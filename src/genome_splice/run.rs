//! Driver: chains in, the fill / refine / pseudo stages, the same output files
//! as the C bathfill (<prefix>.windows.tsv, .segments.tsv, .fills.tsv,
//! .exons.tsv, .chains.tsv, .junctions.tsv, .refined.tsv, .disablements.tsv,
//! .pseudo.tsv, .gff3), genomic 1-based coordinates.

use super::chain::{add_exon, SRC_ALT, SRC_INPUT, SRC_ORF, chain_junctions, chain_windows, load_genome, read_chains, write_chains, Chain, ChainExon, ChainOpts, GAP, KIND_NAME, SRC_NAME, SRC_SPLICE};
use super::hmm::{read_hmms, Hmm};
use super::junction::{junction, JxStatus};
use super::module::alternatives;
use super::orf::{orf_window, Aligner, OrfOpts};
use super::pseudo::Pseudo;
use super::refine::{junction_sites, JxSites, Refine};
use super::gff::write_gff;
use super::score::{score_chain, ExonScore};
use super::splice::{Params, DIS_STOP, MIN_INTRON};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Map f over 0..n on `threads` workers (each with its own aligner); results in index order.
fn par_map<T: Send, F: Fn(usize, &mut Aligner) -> T + Sync>(n: usize, threads: usize, f: F) -> Vec<T> {
    let next = AtomicUsize::new(0);
    let mut parts: Vec<Vec<(usize, T)>> = Vec::new();
    std::thread::scope(|sc| {
        let hs: Vec<_> = (0..threads.max(1).min(n.max(1)))
            .map(|_| sc.spawn(|| {
                let mut al = Aligner::default();
                let mut got = Vec::new();
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= n { break; }
                    got.push((i, f(i, &mut al)));
                }
                got
            }))
            .collect();
        for h in hs { parts.push(h.join().expect("exonfill worker panicked")); }
    });
    let mut all: Vec<(usize, T)> = parts.into_iter().flatten().collect();
    all.sort_by_key(|x| x.0);
    all.into_iter().map(|x| x.1).collect()
}

pub struct Opts {
    pub orf: OrfOpts,
    pub chain: ChainOpts,
    pub prm: Params,
    /// fill: bits a gap exon needs beyond log2 of its search space (gap x missing nodes x 3)
    pub fill_margin: f64,
    pub min_seg: i64,
    pub keep: usize,
    pub min_strict: i64,
    pub splice: bool,
    /// internal recovered exons with rev below this are dropped after refine (-inf: off)
    pub min_rev: f64,
    /// extend gene ends to the start codon and through the stop codon
    pub ends: bool,
    /// score stage: also BATH-style bits and rev for every exon (else frames and nodes only)
    pub score_full: bool,
    /// sites and disablements from the stop-free junction path when it scores within stop_margin and adds no
    /// frameshift; sites only from paths through both anchors (off: the C version)
    pub stop_sites: bool,
    /// half-bits the stop-free path may lose and still be used
    pub stop_margin: f64,
    /// merge two exons when their junction path is one exon with no disablement
    pub join: bool,
    /// further fill rounds for chains that just gained exons: their gaps, or past the end that grew (0: the C version)
    pub flank_rounds: usize,
    /// isoform chains for mutually exclusive alternatives of internal exons (off: the C version)
    pub module: bool,
    /// alternatives: least length ratio and node overlap of the larger span, tested on the
    /// aligned hit, or (alt_refined) on the refined exons after the hit passes a loose gate
    pub alt_size: f64,
    pub alt_nodes_large: f64,
    pub alt_refined: bool,
    /// refine: stop codon score on the junction path (forced_splice's -4)
    pub refine_stop: f64,
    /// refine: frameshift score on the junction path (NaN: fspen)
    pub refine_fs: f64,
    pub threads: usize,
}

impl Default for Opts {
    fn default() -> Self {
        Opts {
            orf: OrfOpts { margin: 5, min_aa: 10, thr: 1.0, near: 100, revthr: 3.0, minm: 12, all: false, decoy: false, overlap: 0, anchored: false, rank: false },
            chain: ChainOpts { min_gap: 10, max_gap: 250, flank: 15000, share: true, siblings: true },
            prm: Params::default(),
            fill_margin: 4.7, min_seg: 30, keep: 30, min_strict: 2, splice: true, min_rev: 0.0, ends: true, score_full: false, stop_sites: true, stop_margin: f64::INFINITY, join: true, flank_rounds: 3, module: true, alt_size: 0.75, alt_nodes_large: 0.75, alt_refined: true, refine_stop: -4.0, refine_fs: f64::NAN,
            threads: std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1),
        }
    }
}

impl Opts {
    fn set(&mut self, k: &str, v: f64) -> Result<(), String> {
        match k {
            "margin" => self.orf.margin = v as i64,
            "min_aa" => self.orf.min_aa = v as usize,
            "bits" => self.orf.thr = v as f32,
            "near" => self.orf.near = v as i64,
            "rev" => self.orf.revthr = v as f32,
            "min_match" => self.orf.minm = v as i32,
            "all" => self.orf.all = v != 0.0,
            "decoy" => self.orf.decoy = v != 0.0,
            "overlap" => self.orf.overlap = v as i64,
            "anchored" => self.orf.anchored = v != 0.0,
            "rank" => self.orf.rank = v != 0.0,
            "min_gap" => self.chain.min_gap = v as i64,
            "max_gap" => self.chain.max_gap = v as i64,
            "flank" => self.chain.flank = v as i64,
            "share_flanks" => self.chain.share = v != 0.0,
            "sibling_mask" => self.chain.siblings = v != 0.0,
            "flank_rounds" => self.flank_rounds = v as usize,
            "module" => self.module = v != 0.0,
            "alt_size" => self.alt_size = v,
            "alt_nodes_large" => self.alt_nodes_large = v,
            "alt_refined" => self.alt_refined = v != 0.0,
            "fill_margin" => self.fill_margin = v,
            "min_seg" => self.min_seg = v as i64,
            "ext" => self.prm.ext = v as i64,
            "w_in" => self.prm.w_in = v as i64,
            "slack" => self.prm.slack = v as i64,
            "stop" => self.prm.stop = v,
            "fspen" => self.prm.fs = v,
            "xsc" => self.prm.xsc = v,
            "keep" => self.keep = v as usize,
            "max_cells" => self.prm.max_cells = v as i64,
            "min_strict" => self.min_strict = v as i64,
            "nosplice" => self.splice = v == 0.0,
            "min_rev" => self.min_rev = v,
            "ends" => self.ends = v != 0.0,
            "score_full" => self.score_full = v != 0.0,
            "stop_sites" => self.stop_sites = v != 0.0,
            "stop_margin" => self.stop_margin = v,
            "join" => self.join = v != 0.0,
            "refine_stop" => self.refine_stop = v,
            "refine_fs" => self.refine_fs = v,
            "skip_open" => self.prm.skip_open = v,
            "cpu" => self.threads = (v as usize).max(1),
            _ => return Err(format!("unknown option '{k}'")),
        }
        Ok(())
    }
}

struct Out {
    files: HashMap<&'static str, String>,
}

impl Out {
    fn buf(&mut self, name: &'static str) -> &mut String { self.files.entry(name).or_default() }
    fn write(&self, prefix: &str) -> Result<(), String> {
        for (name, s) in &self.files {
            let path = if *name == "gff3" { format!("{prefix}.gff3") } else { format!("{prefix}.{name}.tsv") };
            std::fs::write(&path, s).map_err(|e| format!("{path}: {e}"))?;
        }
        Ok(())
    }
}

/// One window's output, assembled in window order afterwards.
#[derive(Default)]
struct WinOut {
    windows: String,
    segments: String,
    fills: String,
    exons: String,
    pend: Vec<(usize, ChainExon, u8)>, // chain, exon, window kind
}

/// seq (starting at genomic position start) with the mask range read as N
fn masked(seq: &[u8], start: i64, mask: Option<(i64, i64)>) -> Vec<u8> {
    let mut s = seq.to_vec();
    if let Some((a, b)) = mask {
        let (i, j) = ((a - start).max(0) as usize, ((b - start + 1).max(0) as usize).min(s.len()));
        if i < j { s[i..j].fill(b'N'); }
    }
    s
}

/// N-terminal (flank_only) rows get no gap windows: copy in a sibling's recovered exons from gaps sharing a flank.
fn copy_to_flank_rows(cs: &mut [Chain]) {
    let same = |x: &ChainExon, y: &ChainExon| x.start == y.start && x.end == y.end;
    let mut add: Vec<(usize, ChainExon)> = Vec::new();
    for (i, c) in cs.iter().enumerate() {
        if !c.flank_only { continue; }
        for s in cs.iter() {
            if s.flank_only || s.passive || s.model != c.model || s.scaffold != c.scaffold || s.strand != c.strand { continue; }
            for w in 1..c.ex.len() {
                let (a, b) = (&c.ex[w - 1], &c.ex[w]);
                if !s.ex.iter().any(|x| same(x, a) || same(x, b)) { continue; }
                let (lo, hi) = (a.end.min(b.end), a.start.max(b.start));
                for x in &s.ex {
                    if (x.src == SRC_SPLICE || x.src == SRC_ORF) && a.k2 < x.k2 && x.k1 < b.k1
                        && x.start - lo - 1 >= MIN_INTRON as i64 && hi - x.end - 1 >= MIN_INTRON as i64 {
                        add.push((i, *x));
                    }
                }
            }
        }
    }
    for (i, x) in add { add_exon(&mut cs[i], &x); }
}

#[allow(clippy::too_many_arguments)]
fn fill_gap(al: &mut Aligner, hid: usize, hmm: &Hmm, o: &Opts, c: &Chain, a: &ChainExon, b: &ChainExon, sc: &[u8],
            mask: Option<(i64, i64)>, jid: &str, ci: usize, wo: &mut WinOut) {
    let lo = a.start.min(b.start);
    let hi = a.end.max(b.end);
    let seq = masked(&sc[(lo - 1) as usize..hi as usize], lo, mask);
    let jx = junction(al, hid, hmm, &o.prm, o.keep, &seq, c.strand,
                      a.k1, a.k2, a.start - lo + 1, a.end - lo + 1, b.k1, b.k2, b.start - lo + 1, b.end - lo + 1);
    let st = c.strand as char;
    if jx.status != JxStatus::Ok {
        let _ = writeln!(wo.fills, "{}\t{}\t{}\t-\t-\t{}\t{}", c.gene, c.scaffold, st, jx.status.name(), jid);
        return;
    }
    let res = &jx.res;
    // evidence in bits against log2 of the search space: gap searched (not masked) x missing nodes x 3 frames
    let e = (res.free - res.null) / 2.0;
    let masked = mask.map_or(0, |(x, y)| (y.min(hi) - x.max(lo) + 1).max(0));
    let space = (jx.gap - masked).max(1) * (b.k1 - a.k2 - 1).max(1) * 3;
    let need = (space as f64).log2() + o.fill_margin;
    let mut segs: Vec<String> = Vec::new();
    for &(s0, s1) in &res.segs {
        if s1 - s0 + 1 < o.min_seg { continue; }
        let (p, q) = (lo - 1 + jx.pos(s0), lo - 1 + jx.pos(s1));
        segs.push(format!("{}-{}", p.min(q), p.max(q)));
    }
    let kept = e >= need && !segs.is_empty();
    let _ = writeln!(wo.fills, "{}\t{}\t{}\t{}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", c.gene, c.scaffold, st,
                     if segs.is_empty() { "-".to_string() } else { segs.join(";") }, e, need, res.free / 2.0, res.null / 2.0, jx.gap,
                     res.nmatch, res.kf, res.kl, res.nfs, res.nstop, res.nnonc, kept as i32, jid);
    let n = res.exons.len();
    let dna = jx.dna();
    for (i, x) in res.exons.iter().enumerate() {
        let role = if i == 0 { "A" } else if i == n - 1 { "B" } else { "gap" };
        let (p, q) = (lo - 1 + jx.pos(x.lo), lo - 1 + jx.pos(x.hi));
        let ac = if x.acc >= 0 { lo - 1 + jx.pos(x.acc) } else { 0 };
        let dn = if x.don >= 0 { lo - 1 + jx.pos(x.don) } else { 0 };
        let asite = if x.acc >= 1 { format!("{}{}", dna[x.acc as usize - 1] as char, dna[x.acc as usize] as char) } else { "--".into() };
        let dsite = if x.don >= 0 && x.don + 1 < jx.d { format!("{}{}", dna[x.don as usize] as char, dna[x.don as usize + 1] as char) } else { "--".into() };
        let _ = writeln!(wo.exons, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", c.gene, c.scaffold, st,
                         jid, i, role, p.min(q), p.max(q), ac, dn, asite, dsite, x.phase, x.nmatch, x.kf, x.kl, kept as i32);
        if kept && i > 0 && i < n - 1 && x.nmatch > 0 && x.hi - x.lo + 1 >= o.min_seg {
            wo.pend.push((ci, ChainExon { start: p.min(q), end: p.max(q), k1: x.kf, k2: x.kl, src: SRC_SPLICE }, GAP));
        }
    }
}

/// One distinct junction's refine/pseudo result.
/// One refined junction, reusable across refine passes (rows are formatted
/// per pass, since exon indices change when exons are dropped).
struct JxOut {
    sites: JxSites,
    pseudo: Pseudo,
    status: &'static str,
    gap: i64,
    dis: Vec<(usize, i64, i64, bool)>, // kind, genomic position, length, tight
    /// permissive path score minus the stop-free path's (NaN: not aligned)
    stopfree: f64,
}

/// Refine every junction of the chains. Junctions already in `cache` are not
/// aligned again. Rows go to `out` when given.
#[allow(clippy::too_many_arguments)]
fn refine_all(models: &[Hmm], mid: &HashMap<String, usize>, o: &Opts, cs: &[Chain], genome: &HashMap<String, Vec<u8>>,
              _do_refine: bool, do_pseudo: bool, cache: &mut HashMap<String, JxOut>, out: Option<&mut Out>) -> (Refine, Vec<Pseudo>) {
    let mut prm = o.prm;
    prm.stop = o.refine_stop; if !o.refine_fs.is_nan() { prm.fs = o.refine_fs; } prm.xsc = 0.0; prm.null_run = false; prm.min_gap_nt = 0; // annotate, as forced_splice
    // stop-free path for splice sites (inserts may not hold stops either)
    let mut bar = prm;
    bar.stop = -1000.0;
    const DIS: [&str; 3] = ["fs", "stop", "noncanonical"];
    // distinct junctions, in order of first appearance
    let pairs = chain_junctions(cs);
    let mut seen: HashMap<String, usize> = HashMap::new();
    let mut uniq: Vec<(usize, usize, usize, String)> = Vec::new();
    let mut which: Vec<Option<usize>> = Vec::with_capacity(pairs.len());
    for &(ci, ia, ib) in &pairs {
        let c = &cs[ci];
        let (a, b) = (&c.ex[ia], &c.ex[ib]);
        if !genome.contains_key(&c.scaffold) || !mid.contains_key(&c.model) { which.push(None); continue; }
        let key = format!("{}|{}|{}|{}|{}|{}|{}", c.model, c.scaffold, c.strand as char, a.start, a.end, b.start, b.end);
        let n = uniq.len();
        let idx = *seen.entry(key.clone()).or_insert(n);
        if idx == n { uniq.push((ci, ia, ib, key)); }
        which.push(Some(idx));
    }
    let todo: Vec<usize> = (0..uniq.len()).filter(|&i| !cache.contains_key(&uniq[i].3)).collect();
    let res: Vec<JxOut> = par_map(todo.len(), o.threads, |t, al| {
        let (ci, ia, ib, _) = &uniq[todo[t]];
        let c = &cs[*ci];
        let (a, b) = (&c.ex[*ia], &c.ex[*ib]);
        let sc = &genome[&c.scaffold];
        let hid = mid[&c.model];
        let lo = a.start.min(b.start);
        let hi = a.end.max(b.end);
        let seq = &sc[(lo - 1) as usize..hi as usize];
        let run = |al: &mut Aligner, p: &Params| junction(al, hid, &models[hid], p, o.keep, seq, c.strand,
                                                          a.k1, a.k2, a.start - lo + 1, a.end - lo + 1, b.k1, b.k2, b.start - lo + 1, b.end - lo + 1);
        let jx = run(al, &prm);
        // a permissive path with no stop is also the best stop-free one
        let has_stop = jx.status == JxStatus::Ok && jx.res.dis.iter().any(|d| d.kind == DIS_STOP);
        let jb = (o.stop_sites && has_stop).then(|| run(al, &bar)).filter(|b| b.status == JxStatus::Ok);
        let stopfree = jb.as_ref().map_or(f64::NAN, |b| jx.res.free - b.res.free);
        // the stop-free path gives sites and disablements unless avoiding the
        // stops costs at least stop_margin; then the stops are real
        // nor when it trades the stops for more frameshifts
        let pj = match &jb { Some(b) if stopfree < o.stop_margin && b.res.nfs <= jx.res.nfs => b, _ => &jx };
        let sites = junction_sites(pj, lo, o.stop_sites);
        let (mut pseudo, mut dis) = (Pseudo::default(), Vec::new());
        if pj.status == JxStatus::Ok {
            pseudo = Pseudo::count(&pj.res);
            dis = pj.res.dis.iter().map(|d| (d.kind as usize, lo - 1 + pj.pos(d.pos), d.n as i64, d.tight)).collect();
        }
        JxOut { sites, pseudo, status: jx.status.name(), gap: jx.gap, dis, stopfree }
    });
    for (t, r) in todo.iter().zip(res) { cache.insert(uniq[*t].3.clone(), r); }
    let mut rf = Refine::new(cs);
    let mut ps = vec![Pseudo::default(); cs.len()];
    for (n, &(ci, ia, ib)) in pairs.iter().enumerate() {
        let Some(idx) = which[n] else { continue };
        let r = &cache[&uniq[idx].3];
        rf.add(ci, ia, ib, &r.sites);
        ps[ci].merge(&r.pseudo);
    }
    if o.ends {
        rf.extend_ends(cs, genome, |c| mid.get(&c.model).map(|&h| models[h].m));
    }
    let Some(out) = out else { return (rf, ps) };
    let mut written = vec![false; uniq.len()];
    for (n, &(ci, ia, ib)) in pairs.iter().enumerate() {
        let Some(idx) = which[n] else { continue };
        if written[idx] { continue; }
        written[idx] = true;
        let (c, r) = (&cs[ci], &cache[&uniq[idx].3]);
        let _ = write!(out.buf("junctions"), "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tx{}\n", c.gene, c.scaffold, c.strand as char, ia, ib,
                       r.status, r.gap, r.sites.n_inner, String::from_utf8_lossy(&r.sites.acc_site),
                       String::from_utf8_lossy(&r.sites.don_site), r.sites.phase, idx);
        if o.stop_sites && !r.stopfree.is_nan() {
            let _ = writeln!(out.buf("stopfree"), "{}\t{}\t{}\tx{}\t{:.2}", c.gene, c.scaffold, c.strand as char, idx, r.stopfree);
        }
        if do_pseudo {
            for &(kind, pos, len, tight) in &r.dis {
                let _ = writeln!(out.buf("disablements"), "{}\t{}\t{}\tx{}\t{}\t{}\t{}\t{}", c.gene, c.scaffold, c.strand as char, idx,
                                 DIS[kind], pos, len, tight as i32);
            }
        }
    }
    if do_pseudo {
        for (i, c) in cs.iter().enumerate() {
            if c.ex.len() > 1 { let s = ps[i].row(c, o.min_strict); out.buf("pseudo").push_str(&s); }
        }
    }
    (rf, ps)
}

/// One isoform chain per mutually exclusive alternative of an internal exon of
/// a base chain (not an IMX row, which shares its base's alternatives): the
/// chain with that exon swapped. Each exon trio is searched once; an
/// alternative already on some chain is not one.
/// alt_refined: length ratio a hit needs to be refined at all
const ALT_LOOSE_SIZE: f64 = 0.5;

/// alt_refined: an alternative chain stays when its refined alternative exon and the
/// refined exon it replaces meet alt_size and alt_nodes_large (residues and node
/// ranges from the score stage). Returns a keep flag per chain.
fn alt_keep(cs: &[Chain], es: &[Vec<ExonScore>], o: &Opts) -> Vec<bool> {
    let idx: HashMap<&str, usize> = cs.iter().enumerate().map(|(i, c)| (c.gene.as_str(), i)).collect();
    let nodes = |x: &ExonScore| {
        let v: Vec<u32> = x.nodes.iter().copied().filter(|&n| n > 0).collect();
        (!v.is_empty()).then(|| (*v.iter().min().unwrap() as i64, *v.iter().max().unwrap() as i64))
    };
    cs.iter().enumerate().map(|(i, c)| {
        let Some(tag) = &c.alt_of else { return true };
        let Some((par, span)) = tag.rsplit_once(':') else { return true };
        let Some((xs, xe)) = span.split_once('-').and_then(|(a, b)| Some((a.parse::<i64>().ok()?, b.parse::<i64>().ok()?))) else { return true };
        let (Some(&pi), Some(j)) = (idx.get(par), c.ex.iter().position(|e| e.src == SRC_ALT)) else { return true };
        let Some(k) = cs[pi].ex.iter().position(|e| e.start <= xe && xs <= e.end) else { return true };
        let (a, x) = (&es[i][j], &es[pi][k]);
        let (Some((a1, a2)), Some((x1, x2))) = (nodes(a), nodes(x)) else { return false };
        let ratio = a.aa.min(x.aa) as f64 / a.aa.max(x.aa).max(1) as f64;
        let ov = (a2.min(x2) - a1.max(x1) + 1).max(0) as f64;
        let larger = ((a2 - a1).max(x2 - x1) + 1) as f64;
        ratio >= o.alt_size && ov >= o.alt_nodes_large * larger
    }).collect()
}

fn add_alternatives(cs: &mut Vec<Chain>, models: &[Hmm], mid: &HashMap<String, usize>, genome: &HashMap<String, Vec<u8>>, o: &Opts) {
    let mut seen: HashMap<String, usize> = HashMap::new();
    let mut todo: Vec<(usize, usize)> = Vec::new();
    for (ci, c) in cs.iter().enumerate() {
        if c.passive || c.imx || c.alt_of.is_some() || !mid.contains_key(&c.model) || !genome.contains_key(&c.scaffold) { continue; }
        for j in 1..c.ex.len().saturating_sub(1) {
            let (p, x, n) = (c.ex[j - 1], c.ex[j], c.ex[j + 1]);
            let key = format!("{}|{}|{}|{}|{}|{}|{}|{}|{}|{}", c.model, c.scaffold, c.strand as char, p.start, p.end, x.start, x.end, n.start, n.end, ci);
            if seen.insert(key, todo.len()).is_none() { todo.push((ci, j)); }
        }
    }
    let found: Vec<Vec<ChainExon>> = par_map(todo.len(), o.threads, |t, al| {
        let (ci, j) = todo[t];
        let c = &cs[ci];
        let hid = mid[&c.model];
        // with the refined test the hit only passes a loose gate here
        let (size, large) = if o.alt_refined { (ALT_LOOSE_SIZE, 0.0) } else { (o.alt_size, o.alt_nodes_large) };
        alternatives(al, hid, &models[hid], c, j, &genome[&c.scaffold], size, large)
    });
    let mut taken: HashMap<(String, u8), Vec<(i64, i64)>> = HashMap::new();
    for c in cs.iter() {
        let v = taken.entry((c.scaffold.clone(), c.strand)).or_default();
        v.extend(c.ex.iter().map(|e| (e.start, e.end)));
    }
    let mut new: Vec<Chain> = Vec::new();
    let mut count: HashMap<usize, usize> = HashMap::new();
    for ((ci, j), alts) in todo.iter().zip(found) {
        for e in alts {
            let v = taken.entry((cs[*ci].scaffold.clone(), cs[*ci].strand)).or_default();
            if v.iter().any(|&(a, b)| e.start <= b && a <= e.end) { continue; }
            v.push((e.start, e.end));
            let c = &cs[*ci];
            let n = count.entry(*ci).or_default();
            *n += 1;
            let mut alt = c.clone();
            alt.gene = format!("{}~m{}", c.gene, n);
            alt.alt_of = Some(format!("{}:{}-{}", c.gene, c.ex[*j].start, c.ex[*j].end));
            alt.ex[*j] = e;
            // An N-terminal row's recovered start pieces belong to the exon being
            // replaced, not to its alternative, which starts the row itself.
            if c.flank_only && c.ex[..*j].iter().all(|x| x.src != SRC_INPUT) {
                alt.ex.drain(..*j);
            }
            new.push(alt);
        }
    }
    cs.extend(new);
}

/// Where the genome comes from: a FASTA path (plain or .gz), or sequences by name.
pub enum GenomeSrc {
    Path(String),
    Seqs(HashMap<String, Vec<u8>>),
}

pub fn run(models_path: &str, chains_path: &str, genome_src: GenomeSrc, prefix: &str, o: &Opts, stages: &str) -> Result<(), String> {
    let mut cs = read_chains(chains_path)?;
    let (models, genome) = std::thread::scope(|sc| {
        let m = sc.spawn(|| read_hmms(models_path, o.threads));
        let g = match genome_src {
            GenomeSrc::Path(p) => load_genome(&p, &cs),
            GenomeSrc::Seqs(g) => Ok(g),
        };
        (m.join().expect("model reader panicked"), g)
    });
    let (models, genome) = (models?, genome?);
    let mut mid: HashMap<String, usize> = HashMap::new();
    for (i, h) in models.iter().enumerate() { mid.entry(h.name.clone()).or_insert(i); }
    let (do_fill, do_refine, do_pseudo) = (stages.contains("fill"), stages.contains("refine"), stages.contains("pseudo"));
    let do_score = stages.contains("score");
    let do_splice = do_fill && o.splice;
    let mut out = Out { files: HashMap::new() };
    if do_fill { out.buf("windows"); out.buf("segments"); }
    if do_splice { out.buf("fills"); out.buf("exons"); }
    if do_fill { out.buf("dropped"); }
    if do_refine || do_pseudo { out.buf("junctions"); }
    if do_refine { out.buf("refined"); }
    if do_pseudo { out.buf("disablements"); out.buf("pseudo"); }
    if do_score { out.buf("exonscores"); }

    if do_fill {
        // a chain that gains exons is searched again: its gaps, or past the end that grew
        let (mut base, mut only): (usize, Option<HashSet<(usize, u8)>>) = (0, None);
        for _round in 0..=o.flank_rounds {
            let wins = chain_windows(&cs, &genome, &o.chain, only.as_ref());
            if wins.is_empty() { break; }
            let outs: Vec<WinOut> = par_map(wins.len(), o.threads, |i, al| {
                let mut wo = WinOut::default();
                let w = &wins[i];
                let c = &cs[w.ci];
                let Some(&hid) = mid.get(&c.model) else { return wo };
                let hmm = &models[hid];
                if w.ge - w.gs < 30 { return wo; }
                let sc = &genome[&c.scaffold];
                let wl = w.ge - w.gs + 1;
                let mut k1 = 1.max(w.k1);
                let mut k2 = (hmm.m as i64).min(w.k2);
                if k2 < k1 { k1 = 1.max((hmm.m as i64).min(k1)); k2 = k1; }
                let id = format!("w{}", base + i);
                let lead = format!("{}\t{}\t{}", c.gene, c.scaffold, c.strand as char);
                let _ = writeln!(wo.windows, "{id}\t{}\t{}\t{}\t{}\t{}\t{k1}\t{k2}\t{}\t{}", KIND_NAME[w.kind as usize], c.gene, c.model,
                                 c.scaffold, c.strand as char, w.gs, w.ge);
                let e = &c.ex[w.ia];
                let (alo, ahi) = if w.kind == GAP { (1, wl) } else if w.ge < e.start { (0, wl) } else { (1, 0) };
                // a gap part that belongs to other isoforms reads as N
                let seg: Vec<u8> = masked(&sc[(w.gs - 1) as usize..w.ge as usize], w.gs, w.mask);
                let mut kept: Vec<ChainExon> = Vec::new();
                orf_window(al, hid, hmm, &o.orf, &id, &lead, w.gs - 1, c.strand, k1, k2, alo, ahi, &seg, &mut wo.segments, &mut kept);
                for k in kept {
                    wo.pend.push((w.ci, k, w.kind));
                    for &c2 in &w.also { wo.pend.push((c2, k, w.kind)); }
                }
                if do_splice && w.kind == GAP {
                    let jid = format!("j{}", base + i);
                    let (a, b) = (c.ex[w.ia], c.ex[w.ib.unwrap()]);
                    fill_gap(al, hid, hmm, o, c, &a, &b, sc, w.mask, &jid, w.ci, &mut wo);
                }
                wo
            });
            let mut pend: Vec<(usize, ChainExon, u8)> = Vec::new();
            for wo in outs {
                out.buf("windows").push_str(&wo.windows);
                out.buf("segments").push_str(&wo.segments);
                if do_splice { out.buf("fills").push_str(&wo.fills); out.buf("exons").push_str(&wo.exons); }
                pend.extend(wo.pend);
            }
            let mut grew: HashSet<(usize, u8)> = HashSet::new();
            for pass in 0..2 {
                for (ci, e, kind) in &pend {
                    if (e.src == SRC_SPLICE) == (pass == 0) && add_exon(&mut cs[*ci], e) { grew.insert((*ci, *kind)); }
                }
            }
            base += wins.len();
            if grew.is_empty() { break; }
            only = Some(grew);
        }
        copy_to_flank_rows(&mut cs);
    }
    let mut cache: HashMap<String, JxOut> = HashMap::new();
    if do_fill && do_refine && o.min_rev.is_finite() {
        // drop internal recovered exons that fit their chain, at refined
        // boundaries, no better than reversed; then refine the new junctions
        let (rf0, _) = refine_all(&models, &mid, o, &cs, &genome, true, false, &mut cache, None);
        let gated = par_map(cs.len(), o.threads, |i, al| {
            let c = &cs[i];
            let n = c.ex.len();
            let test = |j: usize| j > 0 && j + 1 < n && c.ex[j].src != SRC_INPUT;
            if !(0..n).any(test) { return Vec::new(); }
            let (Some(&hid), Some(sc)) = (mid.get(&c.model), genome.get(&c.scaffold)) else { return Vec::new() };
            let spans: Vec<(i64, i64)> = (0..n).map(|j| rf0.span(c, i, j)).collect();
            score_chain(al, hid, &models[hid], c, &spans, sc, false, false, test)
        });
        let d = out.buf("dropped");
        for (c, v) in cs.iter_mut().zip(&gated) {
            if v.is_empty() { continue; }
            let mut j = 0;
            c.ex.retain(|e| {
                let r = v[j].rev;
                j += 1;
                let drop = r < o.min_rev; // NaN (untested) keeps
                if drop {
                    let _ = writeln!(d, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}", c.gene, c.scaffold, c.strand as char, e.start, e.end,
                                     e.k1, e.k2, SRC_NAME[e.src as usize], r);
                }
                !drop
            });
        }
    }
    if do_fill && o.module {
        add_alternatives(&mut cs, &models, &mid, &genome, o);
    }
    let (mut rf, ps) = if do_refine || do_pseudo {
        let (r, p) = refine_all(&models, &mid, o, &cs, &genome, do_refine, do_pseudo, &mut cache, Some(&mut out));
        (Some(r), Some(p))
    } else { (None, None) };
    if let (true, true, Some(r)) = (do_refine, o.join, rf.as_mut()) {
        let s = r.apply_joins(&mut cs);
        out.buf("joins").push_str(&s);
    }
    let mut ps = ps;
    let mut es: Option<Vec<Vec<ExonScore>>> = do_score.then(|| {
        let rf = if do_refine { rf.as_ref() } else { None };
        par_map(cs.len(), o.threads, |i, al| {
            let c = &cs[i];
            let spans: Vec<(i64, i64)> = (0..c.ex.len())
                .map(|j| match rf { Some(r) => r.span(c, i, j), None => (c.ex[j].start, c.ex[j].end) })
                .collect();
            match (mid.get(&c.model), genome.get(&c.scaffold)) {
                (Some(&hid), Some(sc)) if !c.ex.is_empty() => score_chain(al, hid, &models[hid], c, &spans, sc, true, o.score_full, |_| o.score_full),
                _ => vec![ExonScore::default(); c.ex.len()],
            }
        })
    });
    // alternatives judged on their refined exons: failing chains go before anything is written
    if let (true, Some(v)) = (o.module && o.alt_refined, es.as_ref()) {
        let keep = alt_keep(&cs, v, o);
        if keep.iter().any(|&k| !k) {
            let mut it = keep.iter();
            cs.retain(|_| *it.next().unwrap());
            if let Some(r) = rf.as_mut() {
                let mut it = keep.iter();
                r.ex.retain(|_| *it.next().unwrap());
                let mut it = keep.iter();
                r.join.retain(|_| *it.next().unwrap());
            }
            if let Some(p) = ps.as_mut() { let mut it = keep.iter(); p.retain(|_| *it.next().unwrap()); }
            if let Some(v) = es.as_mut() { let mut it = keep.iter(); v.retain(|_| *it.next().unwrap()); }
        }
    }
    if do_fill { out.buf("chains").push_str(&write_chains(&cs)); }
    if let (true, Some(r)) = (do_refine, rf.as_ref()) { let s = r.write(&cs); out.buf("refined").push_str(&s); }
    let rf = if do_refine { rf.as_ref() } else { None };
    if let Some(v) = es.as_ref() {
        let b = out.buf("exonscores");
        b.push_str("# gene exon scaffold strand start end k1 k2 source frame aa bits rev nodes stops\n");
        for (c, v) in cs.iter().zip(v.iter()) {
            for (j, (e, x)) in c.ex.iter().zip(v).enumerate() {
                let _ = writeln!(b, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}\t{}", c.gene, j + 1, c.scaffold, c.strand as char,
                                 e.start, e.end, e.k1, e.k2, SRC_NAME[e.src as usize], x.frame, x.aa, x.bits, x.rev,
                                 x.nodes.iter().map(|n| n.to_string()).collect::<Vec<_>>().join(","),
                                 x.stops.iter().map(|n| n.to_string()).collect::<Vec<_>>().join(","));
            }
        }
    }
    let g = write_gff(&cs, rf, if do_pseudo { ps.as_deref() } else { None }, es.as_deref(), o.min_strict);
    out.buf("gff3").push_str(&g);
    out.write(prefix)
}

/// exonfill_model_cols(models, threads): {model name: 0-based MSA column (MAP) of each match node}.
#[pyfunction]
#[pyo3(signature = (models, threads=1))]
pub fn exonfill_model_cols(py: Python<'_>, models: &str, threads: usize) -> PyResult<HashMap<String, Vec<i64>>> {
    let hs = py.detach(|| read_hmms(models, threads.max(1))).map_err(PyRuntimeError::new_err)?;
    Ok(hs.into_iter().map(|h| (h.name, h.map[1..].iter().map(|&c| c as i64 - 1).collect())).collect())
}

/// exonfill_run(models, chains, genome, prefix, opts=None, stages="fill,refine,pseudo,score");
/// genome is a FASTA path or a dict of scaffold sequences.
#[pyfunction]
#[pyo3(signature = (models, chains, genome, prefix, opts=None, stages="fill,refine,pseudo,score"))]
pub fn exonfill_run(py: Python<'_>, models: &str, chains: &str, genome: &Bound<'_, PyAny>, prefix: &str,
                    opts: Option<HashMap<String, f64>>, stages: &str) -> PyResult<()> {
    let mut o = Opts::default();
    for (k, v) in opts.unwrap_or_default() { o.set(&k, v).map_err(PyRuntimeError::new_err)?; }
    // genome: FASTA path, {scaffold: sequence}, or {scaffold: (length, [(start, sequence), ...])}
    // for scaffolds only partly needed: the rest reads as N (zero bytes, never touched)
    let src = match genome.extract::<String>() {
        Ok(p) => GenomeSrc::Path(p),
        Err(_) => {
            let d = genome.cast::<pyo3::types::PyDict>()?;
            let mut m: HashMap<String, Vec<u8>> = HashMap::new();
            for (k, v) in d.iter() {
                let name: String = k.extract()?;
                let seq = if let Ok(full) = v.extract::<String>() {
                    full.bytes().filter(|c| !c.is_ascii_whitespace()).map(|c| c.to_ascii_uppercase()).collect()
                } else {
                    let (len, segs): (usize, Vec<(usize, String)>) = v.extract()?;
                    let mut buf = vec![0u8; len];
                    for (start, piece) in segs {
                        let at = start.saturating_sub(1).min(len);
                        let n = piece.len().min(len - at);
                        for (b, c) in buf[at..at + n].iter_mut().zip(piece.bytes()) { *b = c.to_ascii_uppercase(); }
                    }
                    buf
                };
                m.insert(name, seq);
            }
            GenomeSrc::Seqs(m)
        }
    };
    py.detach(|| run(models, chains, src, prefix, &o, stages)).map_err(PyRuntimeError::new_err)
}
