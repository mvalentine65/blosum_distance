//! Unilocal alignment of a peptide to a node slice of a model: HMMER's generic
//! profile configuration (p7_ProfileConfig, unilocal), Forward, Backward,
//! posterior decoding and optimal-accuracy traceback, as BATH uses them.
//! HMMER and BATH are BSD 3-clause; their notices are in LICENSE-THIRD-PARTY.
//!
//! Scores are computed in f32 log space with an exact log-sum; profile
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
const LNEG: f32 = f32::NEG_INFINITY;
const FLT_MIN: f32 = f32::MIN_POSITIVE;

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

/// Reusable DP storage (log space): one matrix, used in turn for Forward,
/// Forward plus Backward, the posteriors and the optimal-accuracy fill.
#[derive(Default)]
pub struct Dp {
    fm: Vec<f32>, fi: Vec<f32>, fd: Vec<f32>, fx: Vec<f32>,
    bm: Vec<f32>, bi: Vec<f32>, bd: Vec<f32>, bx: Vec<f32>,
    px: Vec<f32>, ox: Vec<f32>,
}

const XN: usize = 0; const XB: usize = 1; const XE: usize = 2; const XC: usize = 3; const XJ: usize = 4;

fn grow(v: &mut Vec<f32>, n: usize) { if v.len() < n { v.resize(n, 0.0); } }

/// Unilocal Forward/Backward/decoding/optimal-accuracy of `seq` (residue codes):
/// summary and matched (model node, residue) pairs.
pub fn align_both(p: &Profile, seq: &[u8], dp: &mut Dp) -> (AlnHit, Vec<(usize, usize)>) {
    run(p, seq, dp)
}

/// ln(e^a + e^b) in f64; a term 40 nats under the other is below its resolution.
#[inline(always)]
fn lsum64(a: f64, b: f64) -> f64 {
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    let d = lo - hi;
    if !(d >= -40.0) { return hi; } // -inf terms too
    hi + d.exp().ln_1p()
}

/// Forward bits of `seq` in f64 log space, two rows: no path can be lost,
/// however far behind the others it falls.
pub fn fwd_bits_exact(p: &Profile, seq: &[u8]) -> f32 {
    let l = seq.len();
    if l == 0 { return 0.0; }
    let m = p.m;
    let w = m + 1;
    let pmove = 2.0f32 / (l as f32 + 2.0);
    let ploop = 1.0f32 - pmove;
    let (nloop, nmove) = (lnf(ploop as f64), lnf(pmove as f64));
    let (cloop, cmove) = (nloop, nmove);
    let p1 = l as f32 / (l as f32 + 1.0);
    let nullsc = ((l as f32) as f64 * (p1 as f64).ln() + (1.0 - p1 as f64).ln()) as f32 as f64;
    let (mut fm, mut fi, mut fd) = (vec![NEG; 2 * w], vec![NEG; 2 * w], vec![NEG; 2 * w]);
    let (mut xn, mut xb, mut xc) = (0.0f64, nmove, NEG);
    for i in 1..=l {
        let r = seq[i - 1] as usize;
        let (pr, cr) = (((i - 1) & 1) * w, (i & 1) * w);
        fm[cr] = NEG; fi[cr] = NEG; fd[cr] = NEG;
        let mut e = NEG;
        for k in 1..w {
            let mv = lsum64(lsum64(fm[pr + k - 1] + p.mm[k - 1], fi[pr + k - 1] + p.im[k - 1]),
                            lsum64(xb + p.bm[k - 1], fd[pr + k - 1] + p.dm[k - 1])) + p.msc[k][r];
            let iv = if k < m { lsum64(fm[pr + k] + p.mi[k], fi[pr + k] + p.ii[k]) } else { NEG };
            let dv = lsum64(fm[cr + k - 1] + p.md[k - 1], fd[cr + k - 1] + p.dd[k - 1]);
            fm[cr + k] = mv; fi[cr + k] = iv; fd[cr + k] = dv;
            e = lsum64(lsum64(e, mv), dv);
        }
        xn += nloop;
        xc = lsum64(xc + cloop, e);
        xb = xn + nmove;
    }
    ((xc + cmove - nullsc) / std::f64::consts::LN_2) as f32
}

/// ln C(i) per residue (see avx::fwd_lnc), generic path.
pub fn fwd_lnc(p: &Profile, seq: &[u8], dp: &mut Dp) -> Vec<f64> {
    let l = seq.len();
    let mut lnc = vec![0.0f64; l + 1];
    if l == 0 { return lnc; }
    run(p, seq, dp);
    for i in 1..=l { lnc[i] = dp.fx[i * 5 + XC] as f64; }
    lnc
}

/// ln(e^a + e^b). A term 18 nats under the other is below f32 resolution.
#[inline(always)]
fn lsum(a: f32, b: f32) -> f32 {
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    let d = lo - hi;
    if !(d >= -18.0) { return hi; } // -inf terms too
    hi + d.exp().ln_1p()
}

/// Forward and Backward in log space (p7_GForward, p7_GBackward): no scaling,
/// so no limit on how far apart the paths of a row may score.
fn run(p: &Profile, seq: &[u8], dp: &mut Dp) -> (AlnHit, Vec<(usize, usize)>) {
    let l = seq.len();
    let m = p.m;
    if l == 0 { return (AlnHit::default(), Vec::new()); }
    let w = m + 1;
    let cells = (l + 1) * w;
    for v in [&mut dp.fm, &mut dp.fi, &mut dp.fd] { grow(v, cells); }
    for v in [&mut dp.bm, &mut dp.bi, &mut dp.bd] { grow(v, 2 * w); }
    for v in [&mut dp.fx, &mut dp.bx, &mut dp.px, &mut dp.ox] { grow(v, (l + 1) * 5); }

    // length model: unilocal, nj = 0 (J is never entered: E->J is impossible)
    let pmove = 2.0f32 / (l as f32 + 2.0);
    let ploop = 1.0f32 - pmove;
    let (nloop, nmove) = (lnf(ploop as f64) as f32, lnf(pmove as f64) as f32);
    let (cloop, cmove) = (nloop, nmove);
    let p1 = l as f32 / (l as f32 + 1.0);
    let nullsc = ((l as f32) as f64 * (p1 as f64).ln() + (1.0 - p1 as f64).ln()) as f32 as f64;
    let at = |i: usize, k: usize| i * w + k;
    let x = |i: usize, s: usize| i * 5 + s;
    let f = |v: &Vec<f64>| v.iter().map(|&t| t as f32).collect::<Vec<f32>>();
    let (tmm, tim, tdm, tbm) = (f(&p.mm), f(&p.im), f(&p.dm), f(&p.bm));
    let (tmd, tmi, tii, tdd) = (f(&p.md), f(&p.mi), f(&p.ii), f(&p.dd));

    // ---- Forward
    {
        let (fm, fi, fd, fx) = (&mut dp.fm, &mut dp.fi, &mut dp.fd, &mut dp.fx);
        for k in 0..w { fm[k] = LNEG; fi[k] = LNEG; fd[k] = LNEG; }
        fx[x(0, XN)] = 0.0; fx[x(0, XB)] = nmove; fx[x(0, XE)] = LNEG; fx[x(0, XC)] = LNEG; fx[x(0, XJ)] = LNEG;
        for i in 1..=l {
            let r = seq[i - 1] as usize;
            let (pr, cr) = ((i - 1) * w, i * w);
            fm[cr] = LNEG; fi[cr] = LNEG; fd[cr] = LNEG;
            let bprev = fx[x(i - 1, XB)];
            let mut e = LNEG;
            for k in 1..w {
                let mv = lsum(lsum(fm[pr + k - 1] + tmm[k - 1], fi[pr + k - 1] + tim[k - 1]),
                              lsum(bprev + tbm[k - 1], fd[pr + k - 1] + tdm[k - 1])) + p.msc[k][r] as f32;
                let iv = if k < m { lsum(fm[pr + k] + tmi[k], fi[pr + k] + tii[k]) } else { LNEG };
                let dv = lsum(fm[cr + k - 1] + tmd[k - 1], fd[cr + k - 1] + tdd[k - 1]);
                fm[cr + k] = mv; fi[cr + k] = iv; fd[cr + k] = dv;
                e = lsum(lsum(e, mv), dv);
            }
            let nval = fx[x(i - 1, XN)] + nloop;
            fx[x(i, XE)] = e; fx[x(i, XC)] = lsum(fx[x(i - 1, XC)] + cloop, e); fx[x(i, XN)] = nval; fx[x(i, XB)] = nval + nmove; fx[x(i, XJ)] = LNEG;
        }
    }
    let total = dp.fx[x(l, XC)] + cmove;

    // ---- Backward, row i in slot i & 1; each finished row is added to Forward's (M and I)
    {
        let (fm, fi) = (&mut dp.fm, &mut dp.fi);
        let (bm, bi, bd, bx) = (&mut dp.bm, &mut dp.bi, &mut dp.bd, &mut dp.bx);
        let lr = (l & 1) * w;
        bx[x(l, XC)] = cmove; bx[x(l, XE)] = cmove; bx[x(l, XB)] = LNEG; bx[x(l, XN)] = LNEG; bx[x(l, XJ)] = LNEG;
        bm[lr + m] = cmove; bd[lr + m] = cmove; bi[lr + m] = LNEG;
        for k in (1..m).rev() {
            bm[lr + k] = lsum(cmove, bd[lr + k + 1] + tmd[k]);
            bd[lr + k] = lsum(cmove, bd[lr + k + 1] + tdd[k]);
            bi[lr + k] = LNEG;
        }
        for k in 1..w { fm[l * w + k] += bm[lr + k]; fi[l * w + k] += bi[lr + k]; }
        for i in (0..l).rev() {
            let r = seq[i] as usize; // x_{i+1}
            let (cr, nr) = ((i & 1) * w, ((i + 1) & 1) * w);
            let mut b = LNEG;
            for k in 1..w { b = lsum(b, bm[nr + k] + tbm[k - 1] + p.msc[k][r] as f32); }
            let nval = lsum(bx[x(i + 1, XN)] + nloop, b + nmove);
            if i == 0 {
                bx[x(0, XB)] = b;
                bx[x(0, XN)] = nval;
                break;
            }
            let e = bx[x(i + 1, XC)] + cloop;
            bm[cr + m] = e; bd[cr + m] = e; bi[cr + m] = LNEG;
            for k in (1..m).rev() {
                let nm = bm[nr + k + 1] + p.msc[k + 1][r] as f32;
                bm[cr + k] = lsum(lsum(nm + tmm[k], bi[nr + k] + tmi[k]), lsum(e, bd[cr + k + 1] + tmd[k]));
                bi[cr + k] = lsum(nm + tim[k], bi[nr + k] + tii[k]);
                bd[cr + k] = lsum(lsum(nm + tdm[k], bd[cr + k + 1] + tdd[k]), e);
            }
            for k in 1..w { fm[i * w + k] += bm[cr + k]; fi[i * w + k] += bi[cr + k]; }
            bx[x(i, XB)] = b; bx[x(i, XC)] = e; bx[x(i, XE)] = e; bx[x(i, XN)] = nval; bx[x(i, XJ)] = LNEG;
        }
    }

    // ---- Posterior decoding (p7_GDecoding), in place; rows normalised
    {
        let (fx, bx) = (&dp.fx, &dp.bx);
        let (pm, pi, px) = (&mut dp.fm, &mut dp.fi, &mut dp.px);
        for i in 1..=l {
            let cr = i * w;
            pm[cr] = 0.0; pi[cr] = 0.0;
            let mut denom = 0.0f32;
            for k in 1..w {
                pm[cr + k] = (pm[cr + k] - total).exp(); denom += pm[cr + k];
                pi[cr + k] = if k < m { (pi[cr + k] - total).exp() } else { 0.0 }; denom += pi[cr + k];
            }
            px[x(i, XN)] = (fx[x(i - 1, XN)] + bx[x(i, XN)] + nloop - total).exp();
            px[x(i, XJ)] = 0.0;
            px[x(i, XC)] = (fx[x(i - 1, XC)] + bx[x(i, XC)] + cloop - total).exp();
            denom += px[x(i, XN)] + px[x(i, XC)];
            let d = 1.0 / denom;
            for k in 1..w { pm[cr + k] *= d; pi[cr + k] *= d; }
            px[x(i, XN)] *= d; px[x(i, XC)] *= d;
        }
    }

    // ---- Optimal accuracy fill (p7_GOptimalAccuracy), in place
    let td = |v: f64| if v == NEG { FLT_MIN } else { 1.0 };
    let px = &dp.px;
    let (om, oi, od, ox) = (&mut dp.fm, &mut dp.fi, &mut dp.fd, &mut dp.ox);
    for k in 0..w { om[k] = LNEG; oi[k] = LNEG; od[k] = LNEG; }
    ox[x(0, XN)] = 0.0; ox[x(0, XB)] = 0.0; ox[x(0, XE)] = LNEG; ox[x(0, XC)] = LNEG; ox[x(0, XJ)] = LNEG;
    let (t_jl, t_el) = (td(nloop as f64), FLT_MIN); // J loop finite; E->J impossible
    let (t_cl, t_em) = (td(cloop as f64), 1.0);
    let (t_nl, t_nm, t_jm) = (td(nloop as f64), td(nmove as f64), td(nmove as f64));
    for i in 1..=l {
        let (pr, cr) = ((i - 1) * w, i * w);
        om[cr] = LNEG; oi[cr] = LNEG; od[cr] = LNEG;
        let mut e = LNEG;
        for k in 1..m {
            let (pp, ppi) = (om[cr + k], oi[cr + k]);
            let v = f32::max(
                f32::max(td(p.mm[k - 1]) * (om[pr + k - 1] + pp), td(p.im[k - 1]) * (oi[pr + k - 1] + pp)),
                f32::max(td(p.dm[k - 1]) * (od[pr + k - 1] + pp), td(p.bm[k - 1]) * (ox[x(i - 1, XB)] + pp)));
            om[cr + k] = v;
            e = f32::max(e, v);
            oi[cr + k] = f32::max(td(p.mi[k]) * (om[pr + k] + ppi), td(p.ii[k]) * (oi[pr + k] + ppi));
            od[cr + k] = f32::max(td(p.md[k - 1]) * om[cr + k - 1], td(p.dd[k - 1]) * od[cr + k - 1]);
        }
        let pp = om[cr + m];
        om[cr + m] = f32::max(
            f32::max(td(p.mm[m - 1]) * (om[pr + m - 1] + pp), td(p.im[m - 1]) * (oi[pr + m - 1] + pp)),
            f32::max(td(p.dm[m - 1]) * (od[pr + m - 1] + pp), td(p.bm[m - 1]) * (ox[x(i - 1, XB)] + pp)));
        oi[cr + m] = LNEG;
        od[cr + m] = f32::max(td(p.md[m - 1]) * om[cr + m - 1], td(p.dd[m - 1]) * od[cr + m - 1]);
        e = f32::max(e, f32::max(om[cr + m], od[cr + m]));
        ox[x(i, XE)] = e;
        ox[x(i, XJ)] = f32::max(t_jl * (ox[x(i - 1, XJ)] + px[x(i, XJ)]), t_el * e);
        ox[x(i, XC)] = f32::max(t_cl * (ox[x(i - 1, XC)] + px[x(i, XC)]), t_em * e);
        ox[x(i, XN)] = t_nl * (ox[x(i - 1, XN)] + px[x(i, XN)]);
        ox[x(i, XB)] = f32::max(t_nm * ox[x(i, XN)], t_jm * ox[x(i, XJ)]);
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
                let (mut best, mut s, mut kb) = (LNEG, S::S, 0usize);
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
    let mut h = AlnHit { bits: ((total as f64 - nullsc) / std::f64::consts::LN_2) as f32, ..Default::default() };
    let trace: Vec<(usize, usize)> = matches.iter().map(|&(kk, ii)| (kk + p.a - 1, ii)).collect();
    for (n, &(node, ii)) in trace.iter().enumerate() {
        if n == 0 { h.kf = node as i32; h.rf = ii as i32; }
        h.kl = node as i32; h.rl = ii as i32;
        h.nmatch += 1;
    }
    (h, trace)
}
