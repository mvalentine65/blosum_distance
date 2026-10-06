//! ORF stage: align stop-to-stop ORFs to the missing nodes of a window, judge
//! and dedup them; anchor frame detection and trimming.

use super::align::{residue_code, AlnHit, Dp, Profile};
use super::chain::{ChainExon, SRC_ORF};
use super::hmm::Hmm;
use super::sites::{revcomp, translate, Code};
use std::fmt::Write as _;

#[cfg(target_arch = "x86_64")]
pub use super::avx::XRow as ChainRow;
#[cfg(not(target_arch = "x86_64"))]
pub struct ChainRow;

/// Bits by which the fast aligner's scores of a chain may differ from the exact Forward.
const CHAIN_TOL: f32 = 0.005;

pub struct OrfOpts {
    pub margin: i64,
    pub min_aa: usize,
    pub thr: f32,
    pub near: i64,
    pub revthr: f32,
    pub minm: i32,
    pub all: bool,
    pub decoy: bool,
    pub overlap: i64,
    pub anchored: bool,
    /// flank windows: charge log2(1 + ORFs nearer the anchor) instead of log2(ORFs)
    pub rank: bool,
}

/// A profile slice kept while consecutive calls use the same model and nodes.
/// With AVX2 the striped profile runs BATH's SIMD path; otherwise the generic one.
#[derive(Default)]
pub struct Aligner {
    key: Option<(usize, usize, usize)>,
    prof: Option<Profile>,
    #[cfg(target_arch = "x86_64")]
    oprof: Option<super::avx::OProfile>,
    #[cfg(target_arch = "x86_64")]
    avxdp: super::avx::AvxDp,
    pub dp: Dp,
    /// spliced-aligner buffers, reused across junctions
    pub sw: super::splice::SpliceWork,
}

impl Aligner {
    pub fn set(&mut self, hmm_id: usize, hmm: &Hmm, a: usize, b: usize) {
        if self.key == Some((hmm_id, a, b)) {
            return;
        }
        let p = Profile::new(hmm, a, b);
        #[cfg(target_arch = "x86_64")]
        { self.oprof = if super::avx::available() { Some(super::avx::OProfile::new(&p)) } else { None }; }
        self.prof = Some(p);
        self.key = Some((hmm_id, a, b));
    }
    fn run(&mut self, pep: &[u8]) -> (AlnHit, Vec<(usize, usize)>) {
        let seq: Vec<u8> = pep.iter().map(|&c| residue_code(c)).collect();
        let p = self.prof.as_ref().unwrap();
        #[cfg(target_arch = "x86_64")]
        if let Some(om) = self.oprof.as_ref() {
            let r = super::avx::align(p, om, &seq, &mut self.avxdp);
            self.avxdp.trim();
            if let Some(r) = r { return r; }
            // rare: a worker does not keep the generic matrices
            let r = super::align::align_both(p, &seq, &mut self.dp);
            self.dp = Dp::default();
            return r;
        }
        super::align::align_both(p, &seq, &mut self.dp)
    }
    /// ln C per residue of `pep` against the current slice (see avx::fwd_lnc).
    pub fn fwd_lnc(&mut self, pep: &[u8]) -> Vec<f64> {
        let seq: Vec<u8> = pep.iter().map(|&c| residue_code(c)).collect();
        #[cfg(target_arch = "x86_64")]
        if let Some(om) = self.oprof.as_ref() { return super::avx::fwd_lnc(om, &seq, &mut self.avxdp); }
        super::align::fwd_lnc(self.prof.as_ref().unwrap(), &seq, &mut self.dp)
    }
    /// Forward bits of a whole chain `pep`, by a Forward that cannot lose a path.
    /// Also the state after each row of `marks` (ascending), for `fwd_bits_chain_from`.
    pub fn fwd_bits_chain(&self, pep: &[u8], marks: &[usize]) -> (f32, Vec<ChainRow>) {
        let seq: Vec<u8> = pep.iter().map(|&c| residue_code(c)).collect();
        #[cfg(target_arch = "x86_64")]
        if let Some(om) = self.oprof.as_ref() {
            let mut snaps = Vec::with_capacity(marks.len());
            let b = super::avx::fwd_bits_x(om, &seq, 0, None, marks, &mut snaps);
            return (b, snaps);
        }
        let _ = marks;
        (super::align::fwd_bits_exact(self.prof.as_ref().unwrap(), &seq), Vec::new())
    }
    /// The same for a chain that has the first `row` residues of the one `snap`
    /// was taken from (and its length): only the rows after `row` are computed.
    pub fn fwd_bits_chain_from(&self, pep: &[u8], snap: Option<&ChainRow>, row: usize) -> f32 {
        let seq: Vec<u8> = pep.iter().map(|&c| residue_code(c)).collect();
        #[cfg(target_arch = "x86_64")]
        if let (Some(om), Some(s)) = (self.oprof.as_ref(), snap) {
            return super::avx::fwd_bits_x(om, &seq, row, Some(s), &[], &mut Vec::new());
        }
        let _ = (snap, row);
        super::align::fwd_bits_exact(self.prof.as_ref().unwrap(), &seq)
    }
    /// Forward bits of `pep` against the current slice.
    pub fn fwd_bits(&mut self, pep: &[u8]) -> f32 {
        #[cfg(target_arch = "x86_64")]
        if let Some(om) = self.oprof.as_ref() {
            let seq: Vec<u8> = pep.iter().map(|&c| residue_code(c)).collect();
            return super::avx::fwd_bits(om, &seq, &mut self.avxdp);
        }
        self.one(pep).bits
    }
    pub fn one(&mut self, pep: &[u8]) -> AlnHit {
        self.run(pep).0
    }
}

pub fn translate_frame(code: Code, nt: &[u8], f: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(nt.len() / 3 + 1);
    let mut i = f;
    while i + 2 < nt.len() {
        out.push(translate(code, &nt[i..i + 3]));
        i += 3;
    }
    out
}

/// Frame (0-2 from the anchor's first base) whose translation best aligns to
/// nodes k1..k2; stops read as X. With keep > 0, also the kept span next to
/// the junction (side 0: A keeps its tail; side 1: B keeps its head) and the
/// node the cut falls at. Returns (frame relative to lo, lo, hi, k).
pub fn anchor_frame(al: &mut Aligner, hmm_id: usize, hmm: &Hmm, code: Code, nt: &[u8], k1: usize, k2: usize,
                    side: i32, keep: usize, pad: i64) -> (i64, i64, i64, usize) {
    let n = nt.len() as i64;
    let (mut lo, mut hi) = (0i64, n);
    let mut kcut = if side == 0 { k1 } else { k2 };
    al.set(hmm_id, hmm, k1, k2);
    let (mut best, mut bsc) = (0usize, f32::NEG_INFINITY);
    let frames: Vec<Vec<u8>> = (0..3)
        .map(|f| translate_frame(code, nt, f).into_iter().map(|c| if c == b'*' { b'X' } else { c }).collect())
        .collect();
    // frame by Forward alone (the bits a full run reports); only the best one is decoded
    for f in 0..3 {
        if frames[f].is_empty() {
            continue;
        }
        let b = al.fwd_bits(&frames[f]);
        if b > bsc {
            bsc = b;
            best = f;
        }
    }
    if keep > 0 && k2 + 1 - k1 > keep {
        // frame 0 when none scored
        let tr = al.one_trace(&frames[best]);
        let cut = if side == 0 { k2 - keep + 1 } else { k1 + keep - 1 };
        let mut res = 0usize;
        for &(node, i) in &tr {
            if side == 0 && node >= cut { res = i; break; }
            if side == 1 && node <= cut { res = i; }
        }
        if res > 0 {
            if side == 0 { lo = (best as i64 + 3 * (res as i64 - 1) - pad).max(0); } else { hi = (best as i64 + 3 * res as i64 + pad).min(n); }
            kcut = cut;
        }
    }
    ((((best as i64 - lo) % 3) + 3) % 3, lo, hi, kcut)
}

impl Aligner {
    /// Matched (node, residue) pairs of the optimal-accuracy path.
    pub fn one_trace(&mut self, pep: &[u8]) -> Vec<(usize, usize)> {
        self.run(pep).1
    }
    /// The same for a whole chain. The fast aligner's path stands only if its Forward
    /// and Backward scores both match the exact Forward: then neither lost a path.
    pub fn chain_trace(&mut self, pep: &[u8]) -> Vec<(usize, usize)> {
        let seq: Vec<u8> = pep.iter().map(|&c| residue_code(c)).collect();
        let p = self.prof.as_ref().unwrap();
        #[cfg(target_arch = "x86_64")]
        if let Some(om) = self.oprof.as_ref() {
            let r = super::avx::align(p, om, &seq, &mut self.avxdp);
            let bck = self.avxdp.bck_bits;
            self.avxdp.trim();
            if let Some((h, tr)) = r {
                let exact = super::avx::fwd_bits_x(om, &seq, 0, None, &[], &mut Vec::new());
                if (h.bits - exact).abs() <= CHAIN_TOL && (bck - exact).abs() <= CHAIN_TOL { return tr; }
            }
            let r = super::align::align_both(p, &seq, &mut self.dp);
            self.dp = Dp::default();
            return r.1;
        }
        super::align::align_both(p, &seq, &mut self.dp).1
    }
}

/// A hit this much inside short-period runs (period 1-3, >= PERIOD_RUN residues) may be a repeat,
const MAX_PERIODIC: f32 = 0.3;
const PERIOD_RUN: usize = 10;
/// unless the model's consensus over the hit's nodes is at least this much the runs' residues.
const REF_UNIT: f32 = 0.35;

/// Share of residues inside runs repeating the residue p back (p = 1..3, best p), and the runs' residues.
fn periodic(aa: &[u8]) -> (f32, [bool; 21]) {
    let n = aa.len();
    let (mut best, mut unit) = (0usize, [false; 21]);
    for p in 1..=3usize {
        if n <= p { continue; }
        let (mut cov, mut u, mut i) = (0usize, [false; 21], 0usize);
        while i + p < n {
            if aa[i] != aa[i + p] { i += 1; continue; }
            let mut j = i;
            while j + p < n && aa[j] == aa[j + p] { j += 1; }
            if j - i + p >= PERIOD_RUN {
                cov += j - i + p;
                for &c in &aa[i..j + p] { u[residue_code(c) as usize] = true; }
            }
            i = j;
        }
        if cov > best { best = cov; unit = u; }
    }
    (best as f32 / n.max(1) as f32, unit)
}

/// A short-period repeat (e.g. (TA)n read as IYIY) that the model's consensus at nodes k1..k2 does not share.
pub fn is_repeat(aa: &[u8], hmm: &Hmm, k1: i64, k2: i64) -> bool {
    let (share, unit) = periodic(aa);
    if share < MAX_PERIODIC { return false; }
    let (a, b) = (k1.max(1) as usize, (k2.max(0) as usize).min(hmm.m));
    if b < a { return true; }
    let cons = |k: usize| (0..hmm.mat[k].len()).max_by(|&x, &y| hmm.mat[k][x].total_cmp(&hmm.mat[k][y])).unwrap_or(0);
    let n = (a..=b).filter(|&k| unit[cons(k)]).count();
    (n as f32) / ((b - a + 1) as f32) < REF_UNIT
}

/// is_repeat for an exon (oriented nt) in its frame with the fewest stops.
pub fn exon_is_repeat(code: Code, nt: &[u8], hmm: &Hmm, k1: i64, k2: i64) -> bool {
    (0..3).map(|f| translate_frame(code, nt, f)).filter(|t| !t.is_empty())
        .min_by_key(|t| t.iter().filter(|&&c| c == b'*').count())
        .is_some_and(|t| is_repeat(&t, hmm, k1, k2))
}

#[derive(Clone, Copy)]
struct Cand {
    slo: i64,
    shi: i64,
    nmatch: i32,
    kf: i32,
    kl: i32,
    dist: i64,
    kept: bool,
    bits: f32,
    margin: f32,
    rep: bool,
}

/// One gap/flank window (+ strand sequence). Rows go to `out`; kept exons to `kept`.
#[allow(clippy::too_many_arguments)]
pub fn orf_window(al: &mut Aligner, hmm_id: usize, hmm: &Hmm, code: Code, o: &OrfOpts, id: &str, lead: &str, goff: i64,
                  strand: u8, a0: i64, b0: i64, mut alo: i64, mut ahi: i64, wseq: &[u8],
                  out: &mut String, kept: &mut Vec<ChainExon>) {
    let w = wseq.len() as i64;
    let s: Vec<u8> = if strand == b'-' {
        let t1 = if alo != 0 { w - alo + 1 } else { 0 };
        let t2 = if ahi != 0 { w - ahi + 1 } else { 0 };
        alo = t1;
        ahi = t2;
        revcomp(wseq)
    } else {
        wseq.to_vec()
    };
    let up = (alo > 0 && alo <= w / 2) || (ahi > 0 && ahi <= w / 2);
    let dn = alo > w / 2 || ahi > w / 2;
    let a = 1.max(a0 - o.margin - if up { o.overlap } else { 0 });
    let b = (hmm.m as i64).min(b0 + o.margin + if dn { o.overlap } else { 0 });
    if b < a {
        return;
    }
    let one_sided = o.anchored && up != dn;
    let mut cands: Vec<Cand> = Vec::new();
    let mut score = |aa: &[u8], start: i64, cands: &mut Vec<Cand>| {
        let n = aa.len() as i64;
        let ra: Vec<u8> = aa.iter().rev().copied().collect();
        let (mut ca, mut cb) = (a, b);
        if one_sided {
            if up { cb = b.min(a + n + 10) } else { ca = a.max(b - n - 10) }
        }
        al.set(hmm_id, hmm, ca as usize, cb as usize);
        let (fh, rh) = if o.decoy { (al.one(&ra), al.one(aa)) } else { (al.one(aa), al.one(&ra)) };
        let olo = start + 1;
        let ohi = start + 3 * n;
        let slo = if fh.nmatch > 0 { start + 3 * (fh.rf as i64 - 1) + 1 } else { olo };
        let shi = if fh.nmatch > 0 { start + 3 * fh.rl as i64 } else { ohi };
        let mut d = 1i64 << 30;
        if alo > 0 { d = d.min(if slo > alo { slo - alo } else if alo > shi { alo - shi } else { 0 }); }
        if ahi > 0 { d = d.min(if ahi > shi { ahi - shi } else if slo > ahi { slo - ahi } else { 0 }); }
        let fp: &[u8] = if o.decoy { &ra } else { aa };
        let rep = fh.nmatch > 0 && is_repeat(&fp[(fh.rf - 1) as usize..fh.rl as usize], hmm, fh.kf as i64, fh.kl as i64);
        cands.push(Cand { slo, shi, nmatch: fh.nmatch, kf: fh.kf, kl: fh.kl, dist: d, kept: false, bits: fh.bits, margin: fh.bits - rh.bits, rep });
    };
    for f in 0..3usize {
        let mut aa: Vec<u8> = Vec::new();
        let mut start = f as i64;
        let mut i = f;
        while i + 2 < s.len() {
            let c = translate(code, &s[i..i + 3]);
            if c != b'*' {
                aa.push(c);
            } else {
                if aa.len() >= o.min_aa { score(&aa, start, &mut cands); }
                aa.clear();
                start = i as i64 + 3;
            }
            i += 3;
        }
        if aa.len() >= o.min_aa { score(&aa, start, &mut cands); }
    }
    let adj = if cands.len() > 1 { (cands.len() as f32).log2() } else { 0.0 };
    let flank = up != dn;
    for n in 0..cands.len() {
        let a = if o.rank && flank {
            let d = cands[n].dist;
            (1.0 + cands.iter().filter(|x| x.dist < d).count() as f32).log2()
        } else { adj };
        let c = &mut cands[n];
        c.kept = !c.rep && c.margin >= o.revthr && c.nmatch >= o.minm && (c.bits - a >= o.thr || c.dist <= o.near);
    }
    cands.sort_by(|x, y| y.bits.partial_cmp(&x.bits).unwrap_or(std::cmp::Ordering::Equal));
    for k in 0..cands.len() {
        if !cands[k].kept { continue; }
        for kk in 0..k {
            if !cands[kk].kept { continue; }
            let ov = cands[k].kl.min(cands[kk].kl) - cands[k].kf.max(cands[kk].kf) + 1;
            let sh = (cands[k].kl - cands[k].kf).min(cands[kk].kl - cands[kk].kf) + 1;
            if 2 * ov > sh { cands[k].kept = false; break; }
        }
    }
    for c in &cands {
        if !c.kept && !o.all { continue; }
        let (mut slo, mut shi) = (c.slo, c.shi);
        if strand == b'-' { slo = w - c.shi + 1; shi = w - c.slo + 1; }
        let _ = writeln!(out, "{lead}\t{}\t{}\t{:.2}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{id}", goff + slo, goff + shi,
                         c.bits, c.margin, c.nmatch, c.kf, c.kl, c.dist, c.kept as i32);
        if c.kept {
            kept.push(ChainExon { start: goff + slo, end: goff + shi, k1: c.kf as i64, k2: c.kl as i64, src: SRC_ORF, codon: 0 });
        }
    }
}
