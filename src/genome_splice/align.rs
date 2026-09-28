//! Unilocal alignment of a peptide to a node slice of a model: HMMER's generic
//! profile configuration (p7_ProfileConfig, unilocal), Forward, Backward,
//! posterior decoding and optimal-accuracy traceback, as BATH uses them.
//!
//! Scores are computed in f64 log space with an exact log-sum; profile
//! parameters are rounded to f32 as HMMER stores them.

use super::hmm::{Hmm, DD, DM, II, IM, K, MD, MI, MM};

/// HMMER's default amino acid background (p7_AminoFrequencies).
pub const BG: [f32; K] = [
    0.0787945, 0.0151600, 0.0535222, 0.0668298, 0.0397062, 0.0695071, 0.0229198, 0.0590092, 0.0594422, 0.0963728,
    0.0237718, 0.0414386, 0.0482904, 0.0395639, 0.0540978, 0.0683364, 0.0540687, 0.0673417, 0.0114135, 0.0304133,
];

/// Residue code for anything that is not one of the 20 (scored as X).
pub const RX: u8 = 20;
const NEG: f64 = f64::NEG_INFINITY;
const FLT_MIN: f64 = f32::MIN_POSITIVE as f64;

pub fn residue_code(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' => 0, b'C' => 1, b'D' => 2, b'E' => 3, b'F' => 4, b'G' => 5, b'H' => 6, b'I' => 7, b'K' => 8, b'L' => 9,
        b'M' => 10, b'N' => 11, b'P' => 12, b'Q' => 13, b'R' => 14, b'S' => 15, b'T' => 16, b'V' => 17, b'W' => 18,
        b'Y' => 19, _ => RX,
    }
}

/// A unilocal profile of model nodes a..=b (1-based), as bathfill slices it.
pub struct Profile {
    pub a: usize,
    pub m: usize,
    // transition scores out of node k (0..m); bm[k] is entry into M_{k+1}
    pub(super) mm: Vec<f64>, pub(super) im: Vec<f64>, pub(super) dm: Vec<f64>, pub(super) bm: Vec<f64>,
    pub(super) md: Vec<f64>, pub(super) mi: Vec<f64>, pub(super) ii: Vec<f64>, pub(super) dd: Vec<f64>,
    /// match scores [node][residue 0..=20]
    msc: Vec<[f64; 21]>,
    // the same in probability space: transitions and match odds ratios
    pmm: Vec<f64>, pim: Vec<f64>, pdm: Vec<f64>, pbm: Vec<f64>,
    pmd: Vec<f64>, pmi: Vec<f64>, pii: Vec<f64>, pdd: Vec<f64>,
    odds: Vec<[f64; 21]>,
}

fn lnf(x: f64) -> f64 { (x.ln() as f32) as f64 }

impl Profile {
    pub fn msc_at(&self, k: usize, x: usize) -> f64 { self.msc[k][x] }

    pub fn new(hmm: &Hmm, a: usize, b: usize) -> Profile {
        let m = b - a + 1;
        // slice: node 0 takes node a-1's transitions; nodes 1..=m are a..=b
        let st = |k: usize| &hmm.t[a - 1 + k];
        let smat = |k: usize| &hmm.mat[a - 1 + k];

        // local entry from match occupancy
        let mut occ = vec![0f32; m + 1];
        occ[1] = st(0)[MI] + st(0)[MM];
        for k in 2..=m {
            occ[k] = occ[k - 1] * (st(k - 1)[MM] + st(k - 1)[MI]) + (1.0 - occ[k - 1]) * st(k - 1)[DM];
        }
        let mut z = 0f32;
        for k in 1..=m {
            z += occ[k] * (m - k + 1) as f32;
        }
        let mut p = Profile {
            a, m,
            mm: vec![NEG; m + 1], im: vec![NEG; m + 1], dm: vec![NEG; m + 1], bm: vec![NEG; m + 1],
            md: vec![NEG; m + 1], mi: vec![NEG; m + 1], ii: vec![NEG; m + 1], dd: vec![NEG; m + 1],
            msc: vec![[NEG; 21]; m + 1],
            pmm: Vec::new(), pim: Vec::new(), pdm: Vec::new(), pbm: Vec::new(),
            pmd: Vec::new(), pmi: Vec::new(), pii: Vec::new(), pdd: Vec::new(),
            odds: Vec::new(),
        };
        for k in 1..=m {
            p.bm[k - 1] = lnf((occ[k] / z) as f64);
        }
        for k in 1..m {
            let t = st(k);
            p.mm[k] = lnf(t[MM] as f64); p.mi[k] = lnf(t[MI] as f64); p.md[k] = lnf(t[MD] as f64);
            p.im[k] = lnf(t[IM] as f64); p.ii[k] = lnf(t[II] as f64);
            p.dm[k] = lnf(t[DM] as f64); p.dd[k] = lnf(t[DD] as f64);
        }
        for k in 1..=m {
            let mut sc = [0f32; K];
            for x in 0..K {
                sc[x] = ((smat(k)[x] as f64) / (BG[x] as f64)).ln() as f32;
            }
            // X: background-weighted expectation (esl_abc_FExpectScore)
            let (mut r, mut d) = (0f32, 0f32);
            for x in 0..K { r += sc[x] * BG[x]; d += BG[x]; }
            for x in 0..K { p.msc[k][x] = sc[x] as f64; }
            p.msc[k][RX as usize] = (r / d) as f64;
        }
        let ex = |v: &Vec<f64>| v.iter().map(|x| x.exp()).collect::<Vec<f64>>();
        p.pmm = ex(&p.mm); p.pim = ex(&p.im); p.pdm = ex(&p.dm); p.pbm = ex(&p.bm);
        p.pmd = ex(&p.md); p.pmi = ex(&p.mi); p.pii = ex(&p.ii); p.pdd = ex(&p.dd);
        p.odds = p.msc.iter().map(|row| { let mut o = [0f64; 21]; for x in 0..21 { o[x] = row[x].exp(); } o }).collect();
        p
    }
}

/// Alignment summary: bits over the null, match count, node span (model
/// coordinates), first and last matched residue (1-based).
#[derive(Clone, Copy, Default, Debug)]
pub struct AlnHit {
    pub bits: f32,
    pub nmatch: i32,
    pub kf: i32,
    pub kl: i32,
    pub rf: i32,
    pub rl: i32,
}

/// Reusable DP storage (probability space, rows scaled).
#[derive(Default)]
pub struct Dp {
    fm: Vec<f64>, fi: Vec<f64>, fd: Vec<f64>, fx: Vec<f64>,
    bm: Vec<f64>, bi: Vec<f64>, bd: Vec<f64>, bx: Vec<f64>,
    scale: Vec<f64>,
    pm: Vec<f64>, pi: Vec<f64>, px: Vec<f64>,
    om: Vec<f64>, oi: Vec<f64>, od: Vec<f64>, ox: Vec<f64>,
}

const XN: usize = 0; const XB: usize = 1; const XE: usize = 2; const XC: usize = 3; const XJ: usize = 4;

fn grow(v: &mut Vec<f64>, n: usize) { if v.len() < n { v.resize(n, 0.0); } }

/// Unilocal Forward/Backward/decoding/optimal-accuracy of `seq` (residue codes):
/// summary and matched (model node, residue) pairs.
pub fn align_both(p: &Profile, seq: &[u8], dp: &mut Dp) -> (AlnHit, Vec<(usize, usize)>) {
    run(p, seq, dp)
}

/// ln C(i) per residue (see avx::fwd_lnc), generic path.
pub fn fwd_lnc(p: &Profile, seq: &[u8], dp: &mut Dp) -> Vec<f64> {
    let l = seq.len();
    let mut lnc = vec![0.0f64; l + 1];
    if l == 0 { return lnc; }
    run(p, seq, dp);
    let mut cum = 0.0f64;
    for i in 1..=l {
        cum += dp.scale[i].ln();
        lnc[i] = cum + dp.fx[i * 5 + XC].ln();
    }
    lnc
}

/// Forward and Backward run in probability space with each row rescaled (the
/// forward scale factors are reused by Backward, as HMMER's SIMD code does);
/// the result equals the log-space algorithm without a log per cell.
fn run(p: &Profile, seq: &[u8], dp: &mut Dp) -> (AlnHit, Vec<(usize, usize)>) {
    let l = seq.len();
    let m = p.m;
    if l == 0 { return (AlnHit::default(), Vec::new()); }
    let w = m + 1;
    let cells = (l + 1) * w;
    for v in [&mut dp.fm, &mut dp.fi, &mut dp.fd, &mut dp.bm, &mut dp.bi, &mut dp.bd, &mut dp.pm, &mut dp.pi,
              &mut dp.om, &mut dp.oi, &mut dp.od] { grow(v, cells); }
    for v in [&mut dp.fx, &mut dp.bx, &mut dp.px, &mut dp.ox] { grow(v, (l + 1) * 5); }
    grow(&mut dp.scale, l + 1);

    // length model: unilocal, nj = 0 (J is never entered: E->J is impossible)
    let pmove = 2.0f32 / (l as f32 + 2.0);
    let ploop = 1.0f32 - pmove;
    let (nloop, nmove) = (lnf(ploop as f64), lnf(pmove as f64));
    let (cloop, cmove) = (nloop, nmove);
    let (pnl, pnm, pcl, pcm) = (nloop.exp(), nmove.exp(), cloop.exp(), cmove.exp());
    let p1 = l as f32 / (l as f32 + 1.0);
    let nullsc = ((l as f32) as f64 * (p1 as f64).ln() + (1.0 - p1 as f64).ln()) as f32 as f64;
    let at = |i: usize, k: usize| i * w + k;
    let x = |i: usize, s: usize| i * 5 + s;

    // ---- Forward
    {
        let (fm, fi, fd, fx, sc) = (&mut dp.fm, &mut dp.fi, &mut dp.fd, &mut dp.fx, &mut dp.scale);
        for k in 0..w { fm[k] = 0.0; fi[k] = 0.0; fd[k] = 0.0; }
        fx[x(0, XN)] = 1.0; fx[x(0, XB)] = pnm; fx[x(0, XE)] = 0.0; fx[x(0, XC)] = 0.0; fx[x(0, XJ)] = 0.0;
        sc[0] = 1.0;
        for i in 1..=l {
            let o = &p.odds;
            let r = seq[i - 1] as usize;
            let (pr, cr) = ((i - 1) * w, i * w);
            fm[cr] = 0.0; fi[cr] = 0.0; fd[cr] = 0.0;
            let bprev = fx[x(i - 1, XB)];
            let mut e = 0.0;
            let mut mx = 0.0f64;
            for k in 1..w {
                let mv = (fm[pr + k - 1] * p.pmm[k - 1] + fi[pr + k - 1] * p.pim[k - 1]
                        + bprev * p.pbm[k - 1] + fd[pr + k - 1] * p.pdm[k - 1]) * o[k][r];
                let iv = if k < m { fm[pr + k] * p.pmi[k] + fi[pr + k] * p.pii[k] } else { 0.0 };
                let dv = fm[cr + k - 1] * p.pmd[k - 1] + fd[cr + k - 1] * p.pdd[k - 1];
                fm[cr + k] = mv; fi[cr + k] = iv; fd[cr + k] = dv;
                e += mv + dv;
                mx = mx.max(mv).max(iv).max(dv);
            }
            let cval = fx[x(i - 1, XC)] * pcl + e;
            let nval = fx[x(i - 1, XN)] * pnl;
            let bval = nval * pnm;
            mx = mx.max(e).max(cval).max(nval).max(bval);
            let s = if mx > 0.0 { mx } else { 1.0 };
            let inv = 1.0 / s;
            for k in 1..w { fm[cr + k] *= inv; fi[cr + k] *= inv; fd[cr + k] *= inv; }
            fx[x(i, XE)] = e * inv; fx[x(i, XC)] = cval * inv; fx[x(i, XN)] = nval * inv; fx[x(i, XB)] = bval * inv; fx[x(i, XJ)] = 0.0;
            sc[i] = s;
        }
    }
    let lnscale: f64 = dp.scale[1..=l].iter().map(|s| s.ln()).sum();
    let total = dp.fx[x(l, XC)] * pcm; // scaled
    let fwd = total.ln() + lnscale;

    // ---- Backward, rows scaled by the forward factors: row i holds B_i / prod_{j>i} s_j
    {
        let (fx, bm, bi, bd, bx, sc) = (&dp.fx, &mut dp.bm, &mut dp.bi, &mut dp.bd, &mut dp.bx, &dp.scale);
        let _ = fx;
        let lr = l * w;
        bx[x(l, XC)] = pcm; bx[x(l, XE)] = pcm; bx[x(l, XB)] = 0.0; bx[x(l, XN)] = 0.0; bx[x(l, XJ)] = 0.0;
        bm[lr + m] = pcm; bd[lr + m] = pcm; bi[lr + m] = 0.0;
        for k in (1..m).rev() {
            bm[lr + k] = pcm + bd[lr + k + 1] * p.pmd[k];
            bd[lr + k] = pcm + bd[lr + k + 1] * p.pdd[k];
            bi[lr + k] = 0.0;
        }
        for i in (0..l).rev() {
            let r = seq[i] as usize; // x_{i+1}
            let o = &p.odds;
            let (cr, nr) = (i * w, (i + 1) * w);
            let mut b = 0.0;
            for k in 1..w { b += bm[nr + k] * p.pbm[k - 1] * o[k][r]; }
            let inv = 1.0 / sc[i + 1];
            if i == 0 {
                bx[x(0, XB)] = b * inv;
                bx[x(0, XN)] = (bx[x(1, XN)] * pnl + b * pnm) * inv;
                break;
            }
            let cval = bx[x(i + 1, XC)] * pcl;
            let e = cval;
            let nval = bx[x(i + 1, XN)] * pnl + b * pnm;
            bm[cr + m] = e; bd[cr + m] = e; bi[cr + m] = 0.0;
            for k in (1..m).rev() {
                let nm = bm[nr + k + 1] * o[k + 1][r];
                bm[cr + k] = nm * p.pmm[k] + bi[nr + k] * p.pmi[k] + e + bd[cr + k + 1] * p.pmd[k];
                bi[cr + k] = nm * p.pim[k] + bi[nr + k] * p.pii[k];
                bd[cr + k] = nm * p.pdm[k] + bd[cr + k + 1] * p.pdd[k] + e;
            }
            for k in 1..w { bm[cr + k] *= inv; bi[cr + k] *= inv; bd[cr + k] *= inv; }
            bx[x(i, XB)] = b * inv; bx[x(i, XC)] = cval * inv; bx[x(i, XE)] = e * inv; bx[x(i, XN)] = nval * inv; bx[x(i, XJ)] = 0.0;
        }
    }

    // ---- Posterior decoding (p7_GDecoding); rows normalised
    {
        let (fm, fi, fx, bm, bi, bx, sc) = (&dp.fm, &dp.fi, &dp.fx, &dp.bm, &dp.bi, &dp.bx, &dp.scale);
        let (pm, pi, px) = (&mut dp.pm, &mut dp.pi, &mut dp.px);
        let inv_total = 1.0 / total;
        for i in 1..=l {
            let cr = i * w;
            pm[cr] = 0.0; pi[cr] = 0.0;
            let mut denom = 0.0;
            for k in 1..w {
                pm[cr + k] = fm[cr + k] * bm[cr + k] * inv_total; denom += pm[cr + k];
                pi[cr + k] = if k < m { fi[cr + k] * bi[cr + k] * inv_total } else { 0.0 }; denom += pi[cr + k];
            }
            // specials: F_{i-1} B_i misses the factor s_i
            let f = inv_total / sc[i];
            px[x(i, XN)] = fx[x(i - 1, XN)] * bx[x(i, XN)] * pnl * f;
            px[x(i, XJ)] = 0.0;
            px[x(i, XC)] = fx[x(i - 1, XC)] * bx[x(i, XC)] * pcl * f;
            denom += px[x(i, XN)] + px[x(i, XC)];
            let d = 1.0 / denom;
            for k in 1..w { pm[cr + k] *= d; pi[cr + k] *= d; }
            px[x(i, XN)] *= d; px[x(i, XC)] *= d;
        }
    }

    // ---- Optimal accuracy fill (p7_GOptimalAccuracy)
    let td = |v: f64| if v == NEG { FLT_MIN } else { 1.0 };
    let (pm, pi, px) = (&dp.pm, &dp.pi, &dp.px);
    let (om, oi, od, ox) = (&mut dp.om, &mut dp.oi, &mut dp.od, &mut dp.ox);
    for k in 0..w { om[k] = NEG; oi[k] = NEG; od[k] = NEG; }
    ox[x(0, XN)] = 0.0; ox[x(0, XB)] = 0.0; ox[x(0, XE)] = NEG; ox[x(0, XC)] = NEG; ox[x(0, XJ)] = NEG;
    let (t_jl, t_el) = (td(nloop), FLT_MIN); // J loop finite; E->J impossible
    let (t_cl, t_em) = (td(cloop), 1.0);
    let (t_nl, t_nm, t_jm) = (td(nloop), td(nmove), td(nmove));
    for i in 1..=l {
        let (pr, cr) = ((i - 1) * w, i * w);
        om[cr] = NEG; oi[cr] = NEG; od[cr] = NEG;
        let mut e = NEG;
        for k in 1..m {
            let pp = pm[cr + k];
            let v = f64::max(
                f64::max(td(p.mm[k - 1]) * (om[pr + k - 1] + pp), td(p.im[k - 1]) * (oi[pr + k - 1] + pp)),
                f64::max(td(p.dm[k - 1]) * (od[pr + k - 1] + pp), td(p.bm[k - 1]) * (ox[x(i - 1, XB)] + pp)));
            om[cr + k] = v;
            e = f64::max(e, v);
            oi[cr + k] = f64::max(td(p.mi[k]) * (om[pr + k] + pi[cr + k]), td(p.ii[k]) * (oi[pr + k] + pi[cr + k]));
            od[cr + k] = f64::max(td(p.md[k - 1]) * om[cr + k - 1], td(p.dd[k - 1]) * od[cr + k - 1]);
        }
        let pp = pm[cr + m];
        om[cr + m] = f64::max(
            f64::max(td(p.mm[m - 1]) * (om[pr + m - 1] + pp), td(p.im[m - 1]) * (oi[pr + m - 1] + pp)),
            f64::max(td(p.dm[m - 1]) * (od[pr + m - 1] + pp), td(p.bm[m - 1]) * (ox[x(i - 1, XB)] + pp)));
        oi[cr + m] = NEG;
        od[cr + m] = f64::max(td(p.md[m - 1]) * om[cr + m - 1], td(p.dd[m - 1]) * od[cr + m - 1]);
        e = f64::max(e, f64::max(om[cr + m], od[cr + m]));
        ox[x(i, XE)] = e;
        ox[x(i, XJ)] = f64::max(t_jl * (ox[x(i - 1, XJ)] + px[x(i, XJ)]), t_el * e);
        ox[x(i, XC)] = f64::max(t_cl * (ox[x(i - 1, XC)] + px[x(i, XC)]), t_em * e);
        ox[x(i, XN)] = t_nl * (ox[x(i - 1, XN)] + px[x(i, XN)]);
        ox[x(i, XB)] = f64::max(t_nm * ox[x(i, XN)], t_jm * ox[x(i, XJ)]);
    }

    // ---- Traceback (p7_GOATrace): collect match states
    #[derive(Clone, Copy, PartialEq)]
    enum S { M, D, I, N, C, J, E, B, S }
    let (mut i, mut k) = (l, 0usize);
    let mut matches: Vec<(usize, usize)> = Vec::new();
    let mut sprv = S::C;
    let mut guard = 0usize;
    while sprv != S::S {
        guard += 1;
        if guard > 4 * (l + 2) * (m + 2) { break; }
        // a degenerate matrix (all paths tied at the floor) can walk off the edge
        if (matches!(sprv, S::M | S::D) && k == 0) || (matches!(sprv, S::M | S::I) && i == 0) { break; }
        let scur = match sprv {
            S::M => {
                let path = [td(p.mm[k - 1]) * om[at(i - 1, k - 1)], td(p.im[k - 1]) * oi[at(i - 1, k - 1)],
                            td(p.dm[k - 1]) * od[at(i - 1, k - 1)], td(p.bm[k - 1]) * ox[x(i - 1, XB)]];
                let mut a = 0;
                for j in 1..4 { if path[j] > path[a] { a = j; } }
                if k == 1 { a = 3; } // node 1 is entered only from B
                k -= 1; i -= 1;
                [S::M, S::I, S::D, S::B][a]
            }
            S::D => {
                let s = if td(p.md[k - 1]) * om[at(i, k - 1)] >= td(p.dd[k - 1]) * od[at(i, k - 1)] { S::M } else { S::D };
                k -= 1;
                s
            }
            S::I => {
                let s = if td(p.mi[k]) * om[at(i - 1, k)] >= td(p.ii[k]) * oi[at(i - 1, k)] { S::M } else { S::I };
                i -= 1;
                s
            }
            S::N => if i == 0 { S::S } else { S::N },
            S::C => if t_cl * (ox[x(i - 1, XC)] + px[x(i, XC)]) > t_em * ox[x(i, XE)] { S::C } else { S::E },
            S::J => if t_jl * (ox[x(i - 1, XJ)] + px[x(i, XJ)]) > t_el * ox[x(i, XE)] { S::J } else { S::E },
            S::E => {
                let (mut best, mut s, mut kb) = (NEG, S::S, 0usize);
                for kk in 1..=m {
                    if om[at(i, kk)] >= best { best = om[at(i, kk)]; s = S::M; kb = kk; }
                    if od[at(i, kk)] > best { best = od[at(i, kk)]; s = S::D; kb = kk; }
                }
                k = kb;
                s
            }
            S::B => if t_nm * ox[x(i, XN)] > t_jm * ox[x(i, XJ)] { S::N } else { S::J },
            S::S => S::S,
        };
        if scur == S::M { matches.push((k, i)); }
        if (scur == S::N || scur == S::J || scur == S::C) && scur == sprv { i -= 1; }
        sprv = scur;
    }
    matches.reverse();
    let mut h = AlnHit { bits: ((fwd - nullsc) / std::f64::consts::LN_2) as f32, ..Default::default() };
    let trace: Vec<(usize, usize)> = matches.iter().map(|&(kk, ii)| (kk + p.a - 1, ii)).collect();
    for (n, &(node, ii)) in trace.iter().enumerate() {
        if n == 0 { h.kf = node as i32; h.rf = ii as i32; }
        h.kl = node as i32; h.rl = ii as i32;
        h.nmatch += 1;
    }
    (h, trace)
}
