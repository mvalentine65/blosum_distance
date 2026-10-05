//! Forced spliced alignment of a junction between two anchor exons.
//!
//! Nodes k1..k2 (anchor A's first to anchor B's last) are aligned through the
//! locus A..B. Codons inside an anchor keep its frame; no intron may enter an
//! anchor deeper than w_in. Introns score the splice PSSMs of sites, intron open
//! and non-canonical penalties; phase 1/2 introns score the codon assembled
//! across them, node-skipping ones included; 1-2 nt frameshifts and node-skipping introns are penalised. Scores are
//! in half-bits. Intron open, minimum intron and frameshift costs take the
//! values exonerate 2.4.0 uses. The fill stage runs it twice:
//! free, and with no codon starting inside the gap (beyond ext nt of either
//! anchor); the difference is the evidence for new exons. For detection, N
//! codons score xsc and, when stops are barred, inserts may not hold them
//! either, and anchor nodes more than slack from the gap may not place codons
//! in it.

use super::hmm::{Hmm, DD, DM, MD, MM};
use super::sites::{acceptor_default, base, cod64, splice_scores, translate, AccTable, RES_STOP, RES_X};

pub const NEG: f64 = -1e18; // anything above HALF is a real score
pub const HALF: f64 = -1e17;
const HALFBITS: f64 = 2.0 / std::f64::consts::LN_2;
const INTRON_OPEN: f64 = -30.0;
const NONCANON: f64 = -30.0;
pub const MIN_INTRON: usize = 30;
const SKIP_INTRON: usize = 60;
// kernel tiling: ring columns kept per row (> SKIP_INTRON + TILE is not needed:
// rows above are read at most 3 back, own row at most SKIP_INTRON - 1 back)
const RING: usize = 64;
const RM: usize = RING - 1;
const TILE: usize = 32;
// cut: spacer columns for a cut-out gap block; no splice site within WALL_SITES of them (past both PSSMs)
const SPACER: usize = 60;
const WALL_SITES: i64 = 30;

#[derive(Clone, Copy)]
pub struct Params {
    pub stop: f64,
    pub fs: f64,
    /// once per intron that skips model nodes, however many
    pub skip_open: f64,
    pub w_in: i64,
    pub ext: i64,
    pub max_cells: i64,
    pub xsc: f64,
    pub slack: i64,
    pub null_run: bool,
    pub min_gap_nt: i64,
    /// acceptor log-odds
    pub acc: AccTable,
}

impl Default for Params {
    fn default() -> Self {
        Params { stop: -1000.0, fs: -28.0, skip_open: -20.0, w_in: 60, ext: 45, max_cells: 40_000_000, xsc: -4.0, slack: 15, null_run: true, min_gap_nt: 30, acc: acceptor_default() }
    }
}

/// Oriented locus: anchor A start .. anchor B end.
pub struct Locus<'a> {
    pub dna: &'a [u8],
    pub axs: i64, pub axe: i64, // anchor extents, 0-based half-open
    pub bxs: i64, pub bxe: i64,
    pub afo: i64, pub bfo: i64, // frame offsets from the anchor starts
    pub k1: usize, pub k2: usize, // model nodes, A's first .. B's last
    pub ak2: usize, pub bk1: usize, // A's last node, B's first
    /// a block of the gap holding no codon, aligned as (part of) one intron: other isoforms' copies
    pub cut: Option<(i64, i64)>,
}

#[derive(Clone, Copy, Default, Debug)]
pub struct PathExon {
    pub lo: i64, pub hi: i64, // locus coords, inclusive
    pub acc: i64, pub don: i64, // G of AG before, G of GT after; -1 at a path end
    pub phase: i32,
    pub nmatch: i32, pub kf: i64, pub kl: i64,
}

pub const DIS_FS: u8 = 0;
pub const DIS_STOP: u8 = 1;
pub const DIS_NONCANON: u8 = 2;

#[derive(Clone, Copy, Debug)]
pub struct Disable {
    pub kind: u8,
    pub pos: i64,
    pub n: i64,
    pub tight: bool,
}

#[derive(Default, Debug)]
pub struct SpliceResult {
    pub free: f64,
    pub null: f64,
    pub segs: Vec<(i64, i64)>, // gap exon segments, locus coords inclusive
    pub nmatch: i32, pub kf: i64, pub kl: i64,
    pub nfs: i32, pub nstop: i32, pub nnonc: i32,
    pub exons: Vec<PathExon>,
    pub dis: Vec<Disable>,
}

struct Tables {
    lp: usize,
    d: usize,
    em: Vec<[f64; 22]>,
    t: [Vec<f64>; 7],
    codaa: Vec<i32>,
    tb: Vec<u8>,
    don: Vec<f64>,
    acc: Vec<f64>,
    zstart: Vec<i64>,
    cod64: [usize; 64],
}

fn lnhb(p: f32) -> f64 {
    if p > 0.0 { (p as f64).ln() * HALFBITS } else { NEG }
}

fn tables(hmm: &Hmm, loc: &Locus, prm: &Params, wall: Option<(i64, i64)>) -> Tables {
    let lp = loc.k2 - loc.k1 + 1;
    let d = loc.dna.len();
    let mut em = vec![[0f64; 22]; lp];
    for r in 0..lp {
        let node = loc.k1 + r;
        for a in 0..20 {
            let (mv, iv) = (hmm.mat[node][a], hmm.ins(node)[a]);
            em[r][a] = if mv > 0.0 && iv > 0.0 { ((mv / iv) as f64).ln() * HALFBITS } else { NEG };
        }
        em[r][RES_X] = prm.xsc;
        em[r][RES_STOP] = prm.stop;
    }
    let mut t: [Vec<f64>; 7] = Default::default();
    for j in 0..7 { t[j] = vec![NEG; lp + 1]; }
    t[MM][0] = 0.0;
    t[MD][0] = lnhb(hmm.t[loc.k1 - 1][MD]).max((1e-3f64).ln() * HALFBITS);
    t[DD][0] = 0.0;
    for r in 1..lp {
        for j in 0..7 { t[j][r] = lnhb(hmm.t[loc.k1 + r - 1][j]); }
    }
    t[MM][lp] = 0.0;
    t[DM][lp] = 0.0;
    let c64 = cod64();
    let mut tb: Vec<u8> = loc.dna.iter().map(|&c| base(c)).collect();
    tb.push(4);
    let mut codaa = vec![-1i32; d + 1];
    for k in 3..=d {
        let (b1, b2, b3) = (tb[k - 3], tb[k - 2], tb[k - 1]);
        codaa[k] = if b1 < 4 && b2 < 4 && b3 < 4 { c64[(b1 as usize) * 16 + (b2 as usize) * 4 + b3 as usize] as i32 } else { RES_X as i32 };
    }
    let mut ace = loc.axe - prm.w_in;
    let mut bcs = loc.bxs + prm.w_in;
    if ace <= loc.axs { ace = (loc.axs + loc.axe) / 2 + 1; }
    if bcs >= loc.bxe { bcs = (loc.bxs + loc.bxe) / 2; }
    let mut ss5 = vec![0f64; d];
    let mut ss3 = vec![0f64; d];
    splice_scores(&tb[..d], &prm.acc, &mut ss5, &mut ss3);
    let mut don = vec![NEG; d];
    let mut acc = vec![NEG; d];
    let mut zstart = vec![-1i64; d];
    let mut last_core: i64 = -1;
    for k in 0..d {
        let core = (k as i64) < ace || (k as i64) >= bcs;
        let nxt = tb[k + 1];
        let prv = if k > 0 { tb[k - 1] } else { 4 };
        let gt_gc = tb[k] == 2 && (nxt == 3 || nxt == 1);
        let ag = prv == 0 && tb[k] == 2;
        if !core {
            don[k] = ss5[k] + INTRON_OPEN + if gt_gc { 0.0 } else { NONCANON };
            acc[k] = ss3[k] + if ag { 0.0 } else { NONCANON };
        }
        if core { last_core = k as i64; }
        zstart[k] = if core { -1 } else { last_core + 1 };
    }
    if let Some((wlo, whi)) = wall {
        for k in (wlo - WALL_SITES).max(0)..(whi + WALL_SITES).min(d as i64) {
            don[k as usize] = NEG;
            acc[k as usize] = NEG;
        }
    }
    Tables { lp, d, em, t, codaa, tb, don, acc, zstart, cod64: c64 }
}

fn frame_ok(loc: &Locus, flo: i64, fhi: i64) -> Vec<bool> {
    let d = loc.dna.len() as i64;
    (0..=d)
        .map(|k| {
            let x = k - 3;
            if x < 0 { return true; }
            if x >= loc.axs && x < loc.axe { ((x - loc.axs - loc.afo) % 3 + 3) % 3 == 0 }
            else if x >= loc.bxs && x < loc.bxe { ((x - loc.bxs - loc.bfo) % 3 + 3) % 3 == 0 }
            else { !(x >= flo && x < fhi) }
        })
        .collect()
}

// Traceback pointers, packed into one u32 per cell (0 = unset for every field).
const B_BIN: u32 = 0; // 2 bits: bin + 1
const B_ML: u32 = 2; // 3 bits: -1,0,3,4,9 -> 0..4
const B_M: u32 = 5; // 2 bits: -1,0,5,6 -> 0..3
const B_X: u32 = 7; // 2 bits: x + 1
const B_Y: u32 = 9; // 2 bits: y + 1
const B_F: u32 = 11; // 2 bits: f + 1
const B_BP: u32 = 13; // 2 bits: bp + 1
const B_E: u32 = 15; // 2 bits: e (always set)
const B_I: u32 = 17; // 2 bits: i + 1
const B_IL: u32 = 19; // 2 bits: il + 1
const B_A: u32 = 21; // 2 bits: a + 1
const B_IL1: u32 = 23; // 2 bits: il1 + 1
const B_IL2: u32 = 25; // 2 bits: il2 + 1

#[inline(always)]
fn field(w: u32, at: u32, bits: u32) -> u32 { (w >> at) & ((1 << bits) - 1) }
#[inline(always)]
fn dec(w: u32, at: u32) -> i8 { field(w, at, 2) as i8 - 1 }
#[inline(always)]
fn dec_ml(w: u32) -> i8 { [-1, 0, 3, 4, 9][field(w, B_ML, 3) as usize] }
#[inline(always)]
fn dec_m(w: u32) -> i8 { [-1, 0, 5, 6][field(w, B_M, 2) as usize] }

/// Per-worker buffers for the spliced aligner, reused across junctions.
#[derive(Default)]
pub struct SpliceWork {
    ptr: Vec<u32>,
    md: Vec<i32>,
    last: Vec<f64>,
    w: usize,
}

impl SpliceWork {
    fn prepare(&mut self, cells: usize, d: usize) {
        if self.ptr.len() < cells { self.ptr.resize(cells, 0); }
        if self.md.len() < cells { self.md.resize(cells, -1); }
        self.last.clear();
        self.last.resize(d + 1, NEG);
        self.w = d + 1;
    }
}

macro_rules! g { ($v:expr, $i:expr) => { unsafe { *$v.get_unchecked($i) } } }
macro_rules! s { ($v:expr, $i:expr, $x:expr) => { unsafe { *$v.get_unchecked_mut($i) = $x; } } }

/// Best open split-codon entry at an acceptor: the first c maximizing
/// o[c] + a + em[c] over live entries, as the sequential scan picks it.
/// Branchless so it vectorizes.
#[inline(always)]
fn fold<const N: usize>(o: &[f64; N], od: &[i32; N], a: f64, em: &[f64]) -> (f64, i32) {
    let em: &[f64; N] = em[..N].try_into().unwrap();
    let mut v = [f64::NEG_INFINITY; N];
    for c in 0..N { let x = o[c] + a + em[c]; v[c] = if o[c] > HALF { x } else { f64::NEG_INFINITY }; }
    let mut m = f64::NEG_INFINITY;
    for c in 0..N { m = if v[c] > m { v[c] } else { m }; }
    if m <= NEG { return (NEG, -1); }
    let i = (0..N).find(|&c| v[c] == m).unwrap();
    (v[i], od[i])
}

#[inline(always)]
fn fold_any<const N: usize>(avx: bool, o: &[f64; N], od: &[i32; N], a: f64, em: &[f64]) -> (f64, i32) {
    #[cfg(target_arch = "x86_64")]
    if avx { return unsafe { fold_avx(o, od, a, em) }; }
    let _ = avx;
    fold(o, od, a, em)
}

/// fold() in AVX2 lanes; N is 4 or 16.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn fold_avx<const N: usize>(o: &[f64; N], od: &[i32; N], a: f64, em: &[f64]) -> (f64, i32) {
    use std::arch::x86_64::*;
    let em = &em[..N];
    let (ha, av, ninf) = (_mm256_set1_pd(HALF), _mm256_set1_pd(a), _mm256_set1_pd(f64::NEG_INFINITY));
    let mut v = [ninf; 4];
    for j in 0..N / 4 {
        let ov = _mm256_loadu_pd(o.as_ptr().add(4 * j));
        let x = _mm256_add_pd(_mm256_add_pd(ov, av), _mm256_loadu_pd(em.as_ptr().add(4 * j)));
        v[j] = _mm256_blendv_pd(ninf, x, _mm256_cmp_pd::<_CMP_GT_OQ>(ov, ha));
    }
    let mut mv = v[0];
    for j in 1..N / 4 { mv = _mm256_max_pd(mv, v[j]); }
    let t = _mm256_max_pd(mv, _mm256_permute2f128_pd::<1>(mv, mv));
    let t = _mm256_max_pd(t, _mm256_permute_pd::<5>(t));
    let m = _mm256_cvtsd_f64(t);
    if m <= NEG { return (NEG, -1); }
    let mb = _mm256_set1_pd(m);
    let mut bits = 0u32;
    for j in 0..N / 4 { bits |= (_mm256_movemask_pd(_mm256_cmp_pd::<_CMP_EQ_OQ>(v[j], mb)) as u32) << (4 * j); }
    let i = bits.trailing_zeros() as usize;
    // the lane value itself, so a signed zero comes out as the scan would give it
    (o[i] + a + em[i], od[i])
}

/// The Viterbi fill. Row r has consumed r nodes; Bin(r,k) is ready to enter
/// node r+1 with DNA used up to k. A split-codon intron keeps, per zone, the
/// best open entry for each donor-side base (phase 1) or base pair (phase 2);
/// at each acceptor the row folds those into h1/h2, scored with the next node.
/// Every row buffer and pointer word is written once per cell.
#[allow(clippy::too_many_arguments)]
fn kernel(t: &Tables, ok: &[bool], in_gap: &[bool], wl: &[bool], ra: i64, rb: i64, fs: f64, skip_open: f64, ins_stop: bool, p: &mut SpliceWork) {
    #[cfg(target_arch = "x86_64")]
    if super::avx::available() {
        // same code, compiled for AVX2; plain adds and compares, so results match
        return unsafe { kernel_avx2(t, ok, in_gap, wl, ra, rb, fs, skip_open, ins_stop, p) };
    }
    kernel_body::<false>(t, ok, in_gap, wl, ra, rb, fs, skip_open, ins_stop, p)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
unsafe fn kernel_avx2(t: &Tables, ok: &[bool], in_gap: &[bool], wl: &[bool], ra: i64, rb: i64, fs: f64, skip_open: f64, ins_stop: bool, p: &mut SpliceWork) {
    kernel_body::<true>(t, ok, in_gap, wl, ra, rb, fs, skip_open, ins_stop, p)
}

/// wl: bases of the cut wall, which only an intron may consume
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn kernel_body<const AVX: bool>(t: &Tables, ok: &[bool], in_gap: &[bool], wl: &[bool], ra: i64, rb: i64, fs: f64, skip_open: f64, ins_stop: bool, p: &mut SpliceWork) {
    let (lp, d) = (t.lp, t.d);
    let d1 = d + 1;
    let (l_, ls) = (MIN_INTRON as i64, SKIP_INTRON as i64);
    let (neg, half) = (NEG, HALF);
    assert!(ok.len() == d1 && in_gap.len() == d1 && wl.len() == d1 && t.codaa.len() == d1 && t.tb.len() == d1 && t.don.len() == d && t.acc.len() == d);
    assert!(p.ptr.len() >= (lp + 1) * d1 && p.md.len() >= (lp + 1) * d1 && p.last.len() >= d1);
    // Row state lives in per-row rings of RING columns; the fill runs in
    // column tiles, all rows per tile, so long loci stay in cache. Every read
    // reaches at most 59 columns back in its own row or 3 in the row above.
    let nr = (lp + 1) * RING;
    let (mut bin, mut ml, mut xr, mut il) = (vec![neg; nr], vec![neg; nr], vec![neg; nr], vec![neg; nr]);
    let (mut h1, mut h2, mut h1dr, mut h2dr) = (vec![neg; nr], vec![neg; nr], vec![-1i32; nr], vec![-1i32; nr]);
    let (mut yr, mut bp, mut ir) = (vec![neg; nr], vec![neg; nr], vec![neg; nr]);
    let (mut e0, mut e1, mut e2) = (vec![neg; nr], vec![neg; nr], vec![neg; nr]);
    // node-skipping introns after 1 or 2 bases of a split codon
    let (mut il1, mut il2) = (vec![neg; nr], vec![neg; nr]);
    // open split-codon entries per row: o1, o1d, o2, o2d
    let mut ost = vec![([neg; 4], [-1i32; 4], [neg; 16], [-1i32; 16]); lp + 1];
    // per row, emissions regrouped by the acceptor-side bases: [x*4 + b] then 64 + [y*16 + c]
    let mut em12 = vec![[0f64; 128]; lp + 1];
    for r in 0..lp {
        let mut emr = [0f64; 64];
        for i in 0..64 { emr[i] = t.em[r][t.cod64[i]]; }
        for b in 0..4 { for x in 0..16 { em12[r][x * 4 + b] = emr[b * 16 + x]; } }
        for c in 0..16 { for y in 0..4 { em12[r][64 + y * 16 + c] = emr[(c >> 2) * 16 + (c & 3) * 4 + y]; } }
    }
    let tb = &t.tb;
    let (codaa, don, acc, zst) = (&t.codaa, &t.don, &t.acc, &t.zstart);
    let (ptr, mdv, last) = (&mut p.ptr, &mut p.md, &mut p.last);

    for k0 in (0..d1).step_by(TILE) {
        let k1 = (k0 + TILE).min(d1);
        for r in 0..=lp {
            let ro = r * d1;
            let (cr, pr) = (r * RING, r.saturating_sub(1) * RING);
            let (t_mm, t_mi, t_ym, t_yy, t_xm) = (t.t[0][r], t.t[1][r], t.t[3][r], t.t[4][r], t.t[5][r]);
            let emj = if r > 0 { &t.em[r - 1] } else { &t.em[0] };
            let anchor_row = (r as i64) <= ra || (r as i64) >= rb;
            let (t_mx, t_xx) = if r > 0 { (t.t[2][r - 1], t.t[6][r - 1]) } else { (neg, neg) };
            let (mut o1, mut o1d, mut o2, mut o2d) = ost[r];

            for k in k0..k1 {
                let q = ro + k;
                let mut w: u32 = 0;
                let ck = g!(codaa, k);
                let cod_ok = k >= 3 && ck >= 0 && g!(ok, k) && !(anchor_row && g!(in_gap, k));
                let dk = if k < d { g!(don, k) } else { neg };
                let (mut mk, mut xk, mut yk, mut fk, mut ak) = (neg, neg, neg, neg, neg);
                let z = if k < d { g!(zst, k) } else { -1 };
                let ki = k as i64;

                let mut bpv = neg;
                let mut xv = neg;
                if r > 0 {
                    let (mut best, mut ptrv, mut dn) = (neg, 0u32, -1i32);
                    if cod_ok && g!(bin, pr + ((k - 3) & RM)) > half { best = g!(bin, pr + ((k - 3) & RM)) + emj[ck as usize]; ptrv = 1; }
                    if k >= 3 && g!(h1, pr + ((k - 3) & RM)) > best { best = g!(h1, pr + ((k - 3) & RM)); ptrv = 2; dn = g!(h1dr, pr + ((k - 3) & RM)); }
                    if k >= 2 && g!(h2, pr + ((k - 2) & RM)) > best { best = g!(h2, pr + ((k - 2) & RM)); ptrv = 3; dn = g!(h2dr, pr + ((k - 2) & RM)); }
                    if best > half { mk = best; w |= ptrv << B_M; if ptrv > 1 { s!(mdv, q, dn); } }
                    let (vo, ve) = (g!(ml, pr + ((k) & RM)) + t_mx, g!(xr, pr + ((k) & RM)) + t_xx);
                    if vo > half || ve > half {
                        if ve > vo { xk = ve; w |= 2 << B_X; } else { xk = vo; w |= 1 << B_X; }
                        xv = xk;
                    }
                    if mk >= xk { bpv = mk; w |= 1 << B_BP; } else { bpv = xk; w |= 2 << B_BP; }
                }
                s!(xr, cr + ((k) & RM), xv);
                s!(bp, cr + ((k) & RM), bpv);
                let mut yv = neg;
                if cod_ok && (ins_stop || ck as usize != RES_STOP) {
                    // an inserted stop costs what a matched one does
                    let yst = if ck as usize == RES_STOP { emj[RES_STOP] } else { 0.0 };
                    let (vo, ve) = (g!(ml, cr + ((k - 3) & RM)) + t_mi + yst, g!(yr, cr + ((k - 3) & RM)) + t_yy + yst);
                    if vo > half || ve > half {
                        if ve > vo { yk = ve; w |= 2 << B_Y; } else { yk = vo; w |= 1 << B_Y; }
                        yv = yk;
                    }
                }
                s!(yr, cr + ((k) & RM), yv);
                if k >= 1 && !g!(wl, k - 1) {
                    let mut fp = 0u32;
                    if g!(bp, cr + ((k - 1) & RM)) > half { fk = g!(bp, cr + ((k - 1) & RM)) + fs; fp = 2; }
                    if k >= 2 && !g!(wl, k - 2) && g!(bp, cr + ((k - 2) & RM)) > half && g!(bp, cr + ((k - 2) & RM)) + fs > fk { fk = g!(bp, cr + ((k - 2) & RM)) + fs; fp = 3; }
                    w |= fp << B_F;
                }
                let (mut e, mut pe) = (mk, 0u32);
                if xk > e { e = xk; pe = 1; }
                if yk > e { e = yk; pe = 2; }
                if fk > e { e = fk; pe = 3; }
                w |= pe << B_E;
                s!(e0, cr + ((k) & RM), if dk > half && e > half { e + dk } else { neg });
                let (mut icv, mut ilv) = (neg, neg);
                if z >= 0 {
                    let (mut v, mut pp) = (neg, 0u32);
                    if ki - 1 >= z && g!(ir, cr + ((k - 1) & RM)) > half { v = g!(ir, cr + ((k - 1) & RM)); pp = 1; }
                    let dd = ki - l_ + 1;
                    if dd >= z && g!(e0, cr + ((dd as usize) & RM)) > v { v = g!(e0, cr + ((dd as usize) & RM)); pp = 2; }
                    if pp > 0 { icv = v; w |= pp << B_I; }
                    let (mut v, mut pp) = (neg, 0u32);
                    if ki - 1 >= z && g!(il, cr + ((k - 1) & RM)) > half { v = g!(il, cr + ((k - 1) & RM)); pp = 1; }
                    let dd = ki - ls + 1;
                    // a node-skipping intron pays skip_open once, however many nodes it skips
                    if dd >= z && g!(e0, cr + ((dd as usize) & RM)) + skip_open > v { v = g!(e0, cr + ((dd as usize) & RM)) + skip_open; pp = 2; }
                    if r > 0 && g!(il, pr + ((k) & RM)) > half && g!(il, pr + ((k) & RM)) > v { v = g!(il, pr + ((k) & RM)); pp = 3; }
                    if pp > 0 { ilv = v; w |= pp << B_IL; }
                }
                s!(ir, cr + ((k) & RM), icv);
                s!(il, cr + ((k) & RM), ilv);
                // phase 1/2 node-skipping introns open from the e1/e2 donors
                let (mut il1v, mut il2v) = (neg, neg);
                if z >= 0 {
                    let dd = ki - ls + 1;
                    for (ph, st, ev, bit) in [(1usize, &mut il1, &e1, B_IL1), (2usize, &mut il2, &e2, B_IL2)] {
                        let (mut v, mut pp) = (neg, 0u32);
                        if ki - 1 >= z && g!(st, cr + ((k - 1) & RM)) > half { v = g!(st, cr + ((k - 1) & RM)); pp = 1; }
                        if dd >= z && g!(ev, cr + ((dd as usize) & RM)) > half && g!(ev, cr + ((dd as usize) & RM)) + skip_open > v { v = g!(ev, cr + ((dd as usize) & RM)) + skip_open; pp = 2; }
                        if r > 0 && g!(st, pr + ((k) & RM)) > half && g!(st, pr + ((k) & RM)) > v { v = g!(st, pr + ((k) & RM)); pp = 3; }
                        if pp > 0 { if ph == 1 { il1v = v; } else { il2v = v; } w |= pp << bit; }
                    }
                }
                s!(il1, cr + ((k) & RM), il1v);
                s!(il2, cr + ((k) & RM), il2v);
                if k >= 1 && g!(acc, k - 1) > half {
                    let ak1 = g!(acc, k - 1);
                    let mut pa = 0u32;
                    if g!(ir, cr + ((k - 1) & RM)) > half { ak = g!(ir, cr + ((k - 1) & RM)) + ak1; pa = 1; }
                    if g!(il, cr + ((k - 1) & RM)) > half && g!(il, cr + ((k - 1) & RM)) + ak1 > ak { ak = g!(il, cr + ((k - 1) & RM)) + ak1; pa = 2; }
                    w |= pa << B_A;
                }
                let mlk;
                if r == 0 { s!(ml, cr + ((k) & RM), 0.0); w |= 4 << B_ML; mlk = 0.0; } else {
                    let (mut m_, mut pml) = (mk, 1u32);
                    if fk > m_ { m_ = fk; pml = 2; }
                    if ak > m_ { m_ = ak; pml = 3; }
                    s!(ml, cr + ((k) & RM), m_); w |= pml << B_ML; mlk = m_;
                }
                let (mut bi, mut pbi) = (mlk + t_mm, 1u32);
                if yk + t_ym > bi { bi = yk + t_ym; pbi = 2; }
                if xk + t_xm > bi { bi = xk + t_xm; pbi = 3; }
                if bi > half { s!(bin, cr + ((k) & RM), bi); w |= pbi << B_BIN; } else { s!(bin, cr + ((k) & RM), neg); }
                if r == lp { s!(last, k, if bi > half { bi } else { neg }); }
                let (mut e1v, mut e2v) = (neg, neg);
                if dk > half {
                    if k >= 1 && g!(tb, k - 1) < 4 && g!(bin, cr + ((k - 1) & RM)) > half { e1v = g!(bin, cr + ((k - 1) & RM)) + dk; }
                    if k >= 2 && g!(tb, k - 2) < 4 && g!(tb, k - 1) < 4 && g!(bin, cr + ((k - 2) & RM)) > half { e2v = g!(bin, cr + ((k - 2) & RM)) + dk; }
                }
                s!(e1, cr + ((k) & RM), e1v);
                s!(e2, cr + ((k) & RM), e2v);
                let (mut h1v, mut h1d, mut h2v, mut h2d) = (neg, -1i32, neg, -1i32);
                if z >= 0 {
                    if ki - 1 < z || k == 0 {
                        o1 = [neg; 4]; o1d = [-1; 4]; o2 = [neg; 16]; o2d = [-1; 16];
                    }
                    let dd = ki - l_ + 1;
                    if dd >= z {
                        let du = dd as usize;
                        let v1 = g!(e1, cr + ((du) & RM));
                        if v1 > half { let b = g!(tb, du - 1) as usize; if v1 > o1[b] { o1[b] = v1; o1d[b] = dd as i32; } }
                        let v2 = g!(e2, cr + ((du) & RM));
                        if v2 > half { let c = g!(tb, du - 2) as usize * 4 + g!(tb, du - 1) as usize; if v2 > o2[c] { o2[c] = v2; o2d[c] = dd as i32; } }
                    }
                    let ak0 = g!(acc, k);
                    if r < lp && ak0 > half && g!(tb, k + 1) < 4 {
                        if k + 2 <= d && g!(tb, k + 2) < 4 {
                            let x = g!(tb, k + 1) as usize * 4 + g!(tb, k + 2) as usize;
                            (h1v, h1d) = fold_any(AVX, &o1, &o1d, ak0, &em12[r][x * 4..64]);
                        }
                        let y = g!(tb, k + 1) as usize;
                        (h2v, h2d) = fold_any(AVX, &o2, &o2d, ak0, &em12[r][64 + y * 16..]);
                        // a node-skipping split intron; its codon is unscored, -2 marks it for the traceback
                        if k + 2 <= d && g!(tb, k + 2) < 4 && il1v > half && il1v + ak0 > h1v { h1v = il1v + ak0; h1d = -2; }
                        if il2v > half && il2v + ak0 > h2v { h2v = il2v + ak0; h2d = -2; }
                    }
                }
                s!(h1, cr + ((k) & RM), h1v); s!(h1dr, cr + ((k) & RM), h1d);
                s!(h2, cr + ((k) & RM), h2v); s!(h2dr, cr + ((k) & RM), h2d);
                s!(ptr, q, w);
            }
            ost[r] = (o1, o1d, o2, o2d);
        }
    }
}

#[derive(Clone, Copy)]
enum Ev { M { r: usize, k: usize }, S { d: usize, a: usize, ph: usize }, F { k: usize, n: usize }, I { d: usize, a: usize, ph: usize } }

/// ins: where each insert codon ends (k), for the null-run check
struct Trace { ev: Vec<Ev>, ins: Vec<usize>, kstart: usize, kend: usize }

fn traceback(t: &Tables, p: &SpliceWork) -> (f64, Trace) {
    #[derive(Clone, Copy, PartialEq)]
    enum St { Bin, Ml, Bp, E, M, X, Y, F, A, I, Il, Is }
    let w = p.w;
    let (mut r, mut k) = (t.lp, 0usize);
    let mut best = NEG;
    for j in 0..w { if p.last[j] > best { best = p.last[j]; k = j; } }
    let mut tr = Trace { ev: Vec::new(), ins: Vec::new(), kstart: 0, kend: k };
    if best <= HALF { return (best, tr); }
    let mut st = St::Bin;
    let mut a_end: i64 = -1;
    let (mut sa, mut sph) = (0usize, 0usize);
    let limit = 4 * (t.lp + 2) * w;
    for _ in 0..limit {
        let q = r * w + k;
        match st {
            St::Bin => {
                if r == 0 { tr.kstart = k; return (best, tr); }
                st = match dec(p.ptr[q], B_BIN) { 0 => St::Ml, 1 => St::Y, _ => St::X };
            }
            St::Ml => {
                let v = dec_ml(p.ptr[q]);
                if v == 9 { tr.kstart = k; return (best, tr); }
                st = match v { 0 => St::M, 3 => St::F, _ => St::A };
            }
            St::Bp => st = if dec(p.ptr[q], B_BP) == 0 { St::M } else { St::X },
            St::E => st = match field(p.ptr[q], B_E, 2) { 0 => St::M, 1 => St::X, 2 => St::Y, _ => St::F },
            St::M => {
                let v = dec_m(p.ptr[q]);
                if v == 0 {
                    tr.ev.push(Ev::M { r, k });
                    if r == 0 || k < 3 { return (NEG, tr); }
                    r -= 1; k -= 3; st = St::Bin;
                } else {
                    let ph = if v == 5 { 1 } else { 2 };
                    let a = if ph == 1 { k - 3 } else { k - 2 };
                    let d = p.md[q];
                    if d == -2 { sa = a; sph = ph; r -= 1; k = a; st = St::Is; continue; }
                    if d < 0 { return (NEG, tr); }
                    let d = d as usize;
                    tr.ev.push(Ev::S { d, a, ph });
                    tr.ev.push(Ev::I { d, a, ph });
                    r -= 1; k = d - ph; st = St::Bin;
                }
            }
            St::X => { let v = dec(p.ptr[q], B_X); if r == 0 { return (NEG, tr); } r -= 1; st = if v == 1 { St::X } else { St::Ml }; }
            St::Y => { let v = dec(p.ptr[q], B_Y); if k < 3 { return (NEG, tr); } tr.ins.push(k); k -= 3; st = if v == 1 { St::Y } else { St::Ml }; }
            St::F => { let n = dec(p.ptr[q], B_F) as usize; tr.ev.push(Ev::F { k, n }); k -= n; st = St::Bp; }
            St::A => { a_end = k as i64 - 1; st = if dec(p.ptr[q], B_A) == 0 { St::I } else { St::Il }; k -= 1; }
            St::I => {
                match dec(p.ptr[q], B_I) {
                    0 => { loop { if k == 0 { return (NEG, tr); } k -= 1; if dec(p.ptr[r * w + k], B_I) != 0 { break; } } }
                    1 => { let d = k + 1 - MIN_INTRON; tr.ev.push(Ev::I { d, a: a_end as usize, ph: 0 }); k = d; st = St::E; }
                    _ => return (NEG, tr),
                }
            }
            St::Il => {
                match dec(p.ptr[q], B_IL) {
                    0 => { loop { if k == 0 { return (NEG, tr); } k -= 1; if dec(p.ptr[r * w + k], B_IL) != 0 { break; } } }
                    1 => { let d = k + 1 - SKIP_INTRON; tr.ev.push(Ev::I { d, a: a_end as usize, ph: 0 }); k = d; st = St::E; }
                    2 => { if r == 0 { return (NEG, tr); } r -= 1; }
                    _ => return (NEG, tr),
                }
            }
            St::Is => {
                let bit = if sph == 1 { B_IL1 } else { B_IL2 };
                match dec(p.ptr[q], bit) {
                    0 => { loop { if k == 0 { return (NEG, tr); } k -= 1; if dec(p.ptr[r * w + k], bit) != 0 { break; } } }
                    1 => {
                        let d = k + 1 - SKIP_INTRON;
                        tr.ev.push(Ev::S { d, a: sa, ph: sph });
                        tr.ev.push(Ev::I { d, a: sa, ph: sph });
                        k = d - sph; st = St::Bin;
                    }
                    2 => { if r == 0 { return (NEG, tr); } r -= 1; }
                    _ => return (NEG, tr),
                }
            }
        }
    }
    (NEG, tr)
}

fn tight_at(loc: &Locus, x: i64) -> bool {
    if x < loc.axe || x >= loc.bxs { return true; }
    loc.bxs - loc.axe < 30
}

fn summarize(t: &Tables, loc: &Locus, tr: &Trace, res: &mut SpliceResult) {
    let mut pos: Vec<i64> = Vec::new();
    let open = |res: &mut SpliceResult, lo: i64, acc: i64, phase: i32| {
        res.exons.push(PathExon { lo, hi: -1, acc, don: -1, phase, nmatch: 0, kf: 0, kl: 0 });
    };
    open(res, tr.kstart as i64, -1, -1);
    let dna = loc.dna;
    let add = |res: &mut SpliceResult, kind: u8, p: i64, n: i64| {
        let tight = if kind == DIS_NONCANON { true } else { tight_at(loc, p) };
        res.dis.push(Disable { kind, pos: p, n, tight });
    };
    for e in tr.ev.iter().rev() {
        match *e {
            Ev::M { r, k } => {
                let c0 = k as i64 - 3;
                let node = (loc.k1 + r - 1) as i64;
                let x = res.exons.last_mut().unwrap();
                if x.nmatch == 0 { x.kf = node; }
                x.kl = node; x.nmatch += 1;
                if t.codaa[k] == RES_STOP as i32 { res.nstop += 1; add(res, DIS_STOP, c0, 0); }
                if c0 >= loc.axe && c0 < loc.bxs {
                    if res.nmatch == 0 { res.kf = node; }
                    res.kl = node; res.nmatch += 1;
                    for q in c0..k as i64 { pos.push(q); }
                }
            }
            Ev::S { d, a, ph, .. } => {
                let mut c = [0u8; 3];
                let mut j = 0;
                for q in d - ph..d { c[j] = dna[q]; j += 1; }
                let mut q = a + 1;
                while j < 3 { c[j] = dna[q]; j += 1; q += 1; }
                if translate(&c) == b'*' { res.nstop += 1; add(res, DIS_STOP, a as i64 + 1, 0); }
            }
            Ev::F { k, n } => { res.nfs += 1; add(res, DIS_FS, (k - n) as i64, n as i64); }
            Ev::I { d, a, ph } => {
                let gt = dna[d] == b'G' && (dna[d + 1] == b'T' || dna[d + 1] == b'C');
                let ag = a >= 1 && dna[a - 1] == b'A' && dna[a] == b'G';
                if !(gt && ag) { res.nnonc += 1; add(res, DIS_NONCANON, d as i64, 0); }
                let x = res.exons.last_mut().unwrap();
                x.hi = d as i64 - 1; x.don = d as i64;
                open(res, a as i64 + 1, a as i64, ph as i32);
            }
        }
    }
    // stops in inserted codons are stops too
    for &k in &tr.ins {
        if t.codaa[k] == RES_STOP as i32 { res.nstop += 1; add(res, DIS_STOP, k as i64 - 3, 0); }
    }
    res.dis.sort_by_key(|d| d.pos);
    res.exons.last_mut().unwrap().hi = tr.kend as i64 - 1;
    for &q in &pos {
        if let Some(last) = res.segs.last_mut() {
            if q - last.1 <= 6 { if q > last.1 { last.1 = q; } continue; }
        }
        res.segs.push((q, q));
    }
}

pub enum SpliceStatus { Ok, Big, NoPath }

pub fn splice(hmm: &Hmm, loc: &Locus, prm: &Params, res: &mut SpliceResult, work: &mut SpliceWork) -> SpliceStatus {
    // cut WALL_SITES inside the block and ext + WALL_SITES clear of the anchors, so nearby sites survive
    let keep = prm.ext + WALL_SITES;
    let wall = loc.cut.map(|(lo, hi)| ((lo + WALL_SITES).max(loc.axe + keep), (hi - WALL_SITES).min(loc.bxs - keep)));
    let Some((wlo, whi)) = wall.filter(|&(lo, hi)| hi - lo > SPACER as i64) else {
        return splice_on(hmm, loc, prm, res, work, None);
    };
    // the block is cut out; a spacer wall keeps any path across it one intron
    let (l, r) = (wlo as usize, whi as usize);
    let mut dna = Vec::with_capacity(l + SPACER + loc.dna.len() - r);
    dna.extend_from_slice(&loc.dna[..l]);
    dna.resize(l + SPACER, b'N');
    dna.extend_from_slice(&loc.dna[r..]);
    let sh = whi - wlo - SPACER as i64;
    let short = Locus { dna: &dna, bxs: loc.bxs - sh, bxe: loc.bxe - sh, cut: None, ..*loc };
    let st = splice_on(hmm, &short, prm, res, work, Some((wlo, wlo + SPACER as i64)));
    // back to full-locus coords: everything right of the spacer moves by sh
    let at = wlo + SPACER as i64;
    let m = |x: i64| if x >= at { x + sh } else { x };
    for s in &mut res.segs { *s = (m(s.0), m(s.1)); }
    for e in &mut res.exons { e.lo = m(e.lo); e.hi = m(e.hi); e.acc = m(e.acc); e.don = m(e.don); }
    for x in &mut res.dis { x.pos = m(x.pos); }
    st
}

/// wall: locus bases [lo, hi) only an intron may cross, with no splice site near
fn splice_on(hmm: &Hmm, loc: &Locus, prm: &Params, res: &mut SpliceResult, work: &mut SpliceWork, wall: Option<(i64, i64)>) -> SpliceStatus {
    *res = SpliceResult::default();
    let d = loc.dna.len();
    let cells = (loc.k2 - loc.k1 + 2) * (d + 1);
    if cells as i64 > prm.max_cells { return SpliceStatus::Big; }
    let t = tables(hmm, loc, prm, wall);
    work.prepare(cells, d);
    let ra = (loc.ak2 as i64 - loc.k1 as i64 + 1) - prm.slack;
    let rb = (loc.bk1 as i64 - loc.k1 as i64 + 1) + prm.slack;
    let in_gap: Vec<bool> = (0..=d as i64).map(|k| k >= 3 && k - 3 >= loc.axe && k - 3 < loc.bxs).collect();
    let (wlo, whi) = wall.unwrap_or((0, 0));
    let wl: Vec<bool> = (0..=d as i64).map(|k| k >= wlo && k < whi).collect();
    // codon ending at k covers bases k-3..k
    let walled = |ok: &mut Vec<bool>| if whi > wlo { for k in (wlo + 1).max(0)..=(whi + 2).min(d as i64) { ok[k as usize] = false; } };
    let ins_stop = prm.stop > -100.0;
    let flo = loc.axe + prm.ext;
    let fhi = loc.bxs - prm.ext;
    let mut ok = frame_ok(loc, 0, 0);
    walled(&mut ok);
    kernel(&t, &ok, &in_gap, &wl, ra, rb, prm.fs, prm.skip_open, ins_stop, work);
    let (sc, tr) = traceback(&t, work);
    res.free = sc;
    if prm.null_run {
        let fhi = if fhi > flo { fhi } else { flo };
        // a free path with no codon the null run bars is also the null run's best
        let barred = |k: usize| {
            let x = k as i64 - 3;
            x >= flo && x < fhi && !(x >= loc.axs && x < loc.axe) && !(x >= loc.bxs && x < loc.bxe)
        };
        let uses = tr.ev.iter().any(|e| matches!(*e, Ev::M { k, .. } if barred(k))) || tr.ins.iter().any(|&k| barred(k));
        res.null = if fhi == flo || (sc > HALF && !uses) { sc } else {
            let mut ok = frame_ok(loc, flo, fhi);
            walled(&mut ok);
            kernel(&t, &ok, &in_gap, &wl, ra, rb, prm.fs, prm.skip_open, ins_stop, work);
            traceback(&t, work).0
        };
    }
    if res.free > HALF {
        summarize(&t, loc, &tr, res);
        SpliceStatus::Ok
    } else {
        SpliceStatus::NoPath
    }
}
