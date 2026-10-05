//! AVX2 unilocal alignment: BATH's striped (Farrar) Forward, Backward,
//! posterior decoding and optimal-accuracy traceback (impl_avx/fwdback_avx.c,
//! decoding_avx.c, optacc_avx.c, p7_oprofile_avx.c), ported operation for
//! operation so the float results match BATH's. Used when the CPU has AVX2.

#![allow(clippy::missing_safety_doc)]
use super::align::{AlnHit, Profile};
use std::arch::x86_64::*;

const MOVE: usize = 0;
const LOOP: usize = 1;

// cell layout per q: M, D, I
const XM: usize = 0;
const XD: usize = 1;
const XI: usize = 2;
// special cells per row
const SE: usize = 0;
const SN: usize = 1;
const SJ: usize = 2;
const SB: usize = 3;
const SC: usize = 4;
const SSCALE: usize = 5;
const NX: usize = 6;

pub fn available() -> bool {
    is_x86_feature_detected!("avx2")
}

#[inline(always)]
unsafe fn lanes(v: __m256) -> [f32; 8] { std::mem::transmute(v) }

#[target_feature(enable = "avx2")]
unsafe fn rightshiftz(v: __m256) -> __m256 {
    let vi = _mm256_castps_si256(v);
    _mm256_castsi256_ps(_mm256_alignr_epi8::<12>(vi, _mm256_permute2x128_si256::<0x08>(vi, vi)))
}

#[target_feature(enable = "avx2")]
unsafe fn leftshiftz(v: __m256) -> __m256 {
    let vi = _mm256_castps_si256(v);
    _mm256_castsi256_ps(_mm256_alignr_epi8::<4>(_mm256_permute2x128_si256::<0x81>(vi, vi), vi))
}

#[target_feature(enable = "avx2")]
unsafe fn hsum(a: __m256) -> f32 {
    let ai = _mm256_castps_si256(a);
    let t1 = _mm256_castsi256_ps(_mm256_permute2x128_si256::<0x01>(ai, ai));
    let mut t2 = _mm256_add_ps(a, t1);
    let t1 = _mm256_castsi256_ps(_mm256_shuffle_epi32::<0x4e>(_mm256_castps_si256(t2)));
    t2 = _mm256_add_ps(t1, t2);
    let t1 = _mm256_castsi256_ps(_mm256_shuffle_epi32::<0xb1>(_mm256_castps_si256(t2)));
    t2 = _mm256_add_ps(t1, t2);
    f32::from_bits(_mm256_extract_epi32::<0>(_mm256_castps_si256(t2)) as u32)
}

#[target_feature(enable = "avx2")]
unsafe fn hmax(v: __m256) -> f32 {
    let hi = _mm256_extractf128_ps::<1>(v);
    let lo = _mm256_castps256_ps128(v);
    let mut m4 = _mm_max_ps(lo, hi);
    m4 = _mm_max_ps(m4, _mm_shuffle_ps::<0b10_11_00_01>(m4, m4));
    m4 = _mm_max_ps(m4, _mm_shuffle_ps::<0b01_00_11_10>(m4, m4));
    _mm_cvtss_f32(m4)
}

/// esl_sse_expf: Cephes polynomial exp, 4 lanes.
#[target_feature(enable = "avx2")]
unsafe fn sse_expf(mut x: __m128) -> __m128 {
    const P: [f32; 6] = [1.9875691500E-4, 1.3981999507E-3, 8.3334519073E-3, 4.1665795894E-2, 1.6666665459E-1, 5.0000001201E-1];
    const CC: [f32; 2] = [0.693359375, -2.12194440e-4];
    const MAXLOGF: f32 = 88.3762626647949;
    const MINLOGF: f32 = -88.3762626647949;
    let maxmask = _mm_cmpgt_ps(x, _mm_set1_ps(MAXLOGF));
    let minmask = _mm_cmple_ps(x, _mm_set1_ps(MINLOGF));
    let mut fx = _mm_mul_ps(x, _mm_set1_ps(1.44269504088896341f64 as f32));
    fx = _mm_add_ps(fx, _mm_set1_ps(0.5));
    let mut k = _mm_cvttps_epi32(fx);
    let mut tmp = _mm_cvtepi32_ps(k);
    let mut mask = _mm_cmpgt_ps(tmp, fx);
    mask = _mm_and_ps(mask, _mm_set1_ps(1.0));
    fx = _mm_sub_ps(tmp, mask);
    k = _mm_cvttps_epi32(fx);
    tmp = _mm_mul_ps(fx, _mm_set1_ps(CC[0]));
    let mut z = _mm_mul_ps(fx, _mm_set1_ps(CC[1]));
    x = _mm_sub_ps(x, tmp);
    x = _mm_sub_ps(x, z);
    z = _mm_mul_ps(x, x);
    let mut y = _mm_set1_ps(P[0]);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, _mm_set1_ps(P[1])); y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, _mm_set1_ps(P[2])); y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, _mm_set1_ps(P[3])); y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, _mm_set1_ps(P[4])); y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, _mm_set1_ps(P[5])); y = _mm_mul_ps(y, z);
    y = _mm_add_ps(y, x);
    y = _mm_add_ps(y, _mm_set1_ps(1.0));
    k = _mm_add_epi32(k, _mm_set1_epi32(127));
    k = _mm_slli_epi32::<23>(k);
    fx = _mm_castsi128_ps(k);
    y = _mm_mul_ps(y, fx);
    let sel = |a: __m128, b: __m128, m: __m128| _mm_or_ps(_mm_andnot_ps(m, a), _mm_and_ps(m, b));
    y = sel(y, _mm_set1_ps(f32::INFINITY), maxmask);
    y = sel(y, _mm_set1_ps(0.0), minmask);
    y
}

#[target_feature(enable = "avx2")]
unsafe fn avx_expf(v: [f32; 8]) -> __m256 {
    let v = _mm256_loadu_ps(v.as_ptr());
    let lo = sse_expf(_mm256_castps256_ps128(v));
    let hi = sse_expf(_mm256_extractf128_ps::<1>(v));
    _mm256_set_m128(hi, lo)
}

/// Striped profile (p7_oprofile fb_conversion) for one node slice.
pub struct OProfile {
    m: usize,
    q: usize,
    rfv: Vec<__m256>, // 21 residues x q
    tfv: Vec<__m256>, // 7 per q, then q of DD
    e_loop: f32,
    e_move: f32,
}

impl OProfile {
    pub fn new(p: &Profile) -> OProfile {
        unsafe { Self::build(p) }
    }

    #[target_feature(enable = "avx2")]
    unsafe fn build(p: &Profile) -> OProfile {
        let m = p.m;
        let nq = 2usize.max((m - 1) / 8 + 1);
        let ninf = f32::NEG_INFINITY;
        let mut rfv = Vec::with_capacity(21 * nq);
        for x in 0..21 {
            for q in 0..nq {
                let k = q + 1;
                let mut t = [ninf; 8];
                for z in 0..8 {
                    if k + z * nq <= m { t[z] = p.msc_at(k + z * nq, x) as f32; }
                }
                rfv.push(avx_expf(t));
            }
        }
        let mut tfv = Vec::with_capacity(8 * nq);
        for q in 0..nq {
            let k = q + 1;
            for t in 0..7 {
                let (kb, arr): (usize, &Vec<f64>) = match t {
                    0 => (k - 1, &p.bm), 1 => (k - 1, &p.mm), 2 => (k - 1, &p.im), 3 => (k - 1, &p.dm),
                    4 => (k, &p.md), 5 => (k, &p.mi), _ => (k, &p.ii),
                };
                let mut v = [ninf; 8];
                for z in 0..8 {
                    if kb + z * nq < m { v[z] = arr[kb + z * nq] as f32; }
                }
                tfv.push(avx_expf(v));
            }
        }
        for q in 0..nq {
            let k = q + 1;
            let mut v = [ninf; 8];
            for z in 0..8 {
                if k + z * nq < m { v[z] = p.dd[k + z * nq] as f32; }
            }
            tfv.push(avx_expf(v));
        }
        OProfile { m, q: nq, rfv, tfv, e_loop: 0.0f32, e_move: 1.0f32 }
    }
}

/// Full DP storage, rows 0..=L.
#[derive(Default)]
pub struct Omx {
    dp: Vec<__m256>,
    xmx: Vec<f32>,
    totscale: f64,
    has_own_scales: bool,
}

impl Omx {
    fn grow(&mut self, rows: usize, q: usize) {
        self.grow2(rows, rows, q);
    }
    /// dp rows and special-state rows apart: a rolling Forward keeps two dp rows
    fn grow2(&mut self, rows: usize, xrows: usize, q: usize) {
        let n = rows * q * 3;
        if self.dp.len() < n { self.dp.resize(n, unsafe { _mm256_setzero_ps() }); }
        if self.xmx.len() < xrows * NX { self.xmx.resize(xrows * NX, 0.0); }
    }
}

#[derive(Default)]
pub struct AvxDp {
    fwd: Omx,
    bck: Omx,
    /// Backward's score of the last alignment, in bits as `AlnHit::bits` is.
    pub bck_bits: f32,
}

/// A matrix above this is dropped after its alignment, so no worker holds its largest to the end.
const KEEP_BYTES: usize = 32 << 20;

impl AvxDp {
    pub fn trim(&mut self) {
        for ox in [&mut self.fwd, &mut self.bck] {
            if ox.dp.capacity() * std::mem::size_of::<__m256>() > KEEP_BYTES { ox.dp = Vec::new(); }
        }
    }
}

struct Xf { n: [f32; 2], c: [f32; 2], j: [f32; 2], e: [f32; 2] }

fn xf_for(om: &OProfile, l: usize) -> Xf {
    let pmove = 2.0f32 / (l as f32 + 2.0f32);
    let ploop = 1.0f32 - pmove;
    Xf { n: [pmove, ploop], c: [pmove, ploop], j: [pmove, ploop], e: [om.e_move, om.e_loop] }
}

/// ROLL: two dp rows (row i in slot i & 1), for callers that want only the
/// score or the special states; the arithmetic is the same.
#[target_feature(enable = "avx2")]
unsafe fn forward<const ROLL: bool>(dsq: &[u8], om: &OProfile, xf: &Xf, ox: &mut Omx) -> f32 {
    let l = dsq.len() - 1; // dsq[1..=l]
    let q_ = om.q;
    let row = q_ * 3;
    let zerov = _mm256_setzero_ps();
    let dp = &mut ox.dp;
    for q in 0..q_ { dp[q * 3 + XM] = zerov; dp[q * 3 + XI] = zerov; dp[q * 3 + XD] = zerov; }
    let x = &mut ox.xmx;
    let mut xe = 0.0f32; let mut xn = 1.0f32; let mut xj = 0.0f32; let mut xb = xf.n[MOVE]; let mut xc = 0.0f32;
    x[SE] = xe; x[SN] = xn; x[SJ] = xj; x[SB] = xb; x[SC] = xc; x[SSCALE] = 1.0;
    ox.totscale = 0.0;
    for i in 1..=l {
        let (pc, pp) = if ROLL { ((i & 1) * row, ((i - 1) & 1) * row) } else { (i * row, (i - 1) * row) };
        let rp = dsq[i] as usize * q_;
        let mut tp = 0usize;
        let mut dcv = zerov;
        let mut xev = zerov;
        let xbv = _mm256_set1_ps(xb);
        let mut mpv = rightshiftz(dp[pp + (q_ - 1) * 3 + XM]);
        let mut dpv = rightshiftz(dp[pp + (q_ - 1) * 3 + XD]);
        let mut ipv = rightshiftz(dp[pp + (q_ - 1) * 3 + XI]);
        for q in 0..q_ {
            let t = &om.tfv;
            let mut sv = _mm256_mul_ps(xbv, t[tp]); tp += 1;
            sv = _mm256_add_ps(sv, _mm256_mul_ps(mpv, t[tp])); tp += 1;
            sv = _mm256_add_ps(sv, _mm256_mul_ps(ipv, t[tp])); tp += 1;
            sv = _mm256_add_ps(sv, _mm256_mul_ps(dpv, t[tp])); tp += 1;
            sv = _mm256_mul_ps(sv, om.rfv[rp + q]);
            xev = _mm256_add_ps(xev, sv);
            mpv = dp[pp + q * 3 + XM];
            dpv = dp[pp + q * 3 + XD];
            ipv = dp[pp + q * 3 + XI];
            dp[pc + q * 3 + XM] = sv;
            dp[pc + q * 3 + XD] = dcv;
            dcv = _mm256_mul_ps(sv, t[tp]); tp += 1;
            let sv2 = _mm256_mul_ps(mpv, t[tp]); tp += 1;
            dp[pc + q * 3 + XI] = _mm256_add_ps(sv2, _mm256_mul_ps(ipv, t[tp])); tp += 1;
        }
        dcv = rightshiftz(dcv);
        dp[pc + XD] = zerov;
        let dd0 = 7 * q_;
        for q in 0..q_ {
            dp[pc + q * 3 + XD] = _mm256_add_ps(dcv, dp[pc + q * 3 + XD]);
            dcv = _mm256_mul_ps(dp[pc + q * 3 + XD], om.tfv[dd0 + q]);
        }
        if om.m < 100 {
            for _ in 1..4 {
                dcv = rightshiftz(dcv);
                for q in 0..q_ {
                    dp[pc + q * 3 + XD] = _mm256_add_ps(dcv, dp[pc + q * 3 + XD]);
                    dcv = _mm256_mul_ps(dcv, om.tfv[dd0 + q]);
                }
            }
        } else {
            for _ in 1..4 {
                dcv = rightshiftz(dcv);
                let mut cv = zerov;
                for q in 0..q_ {
                    let sv = _mm256_add_ps(dcv, dp[pc + q * 3 + XD]);
                    cv = _mm256_or_ps(cv, _mm256_cmp_ps::<_CMP_GT_OS>(sv, dp[pc + q * 3 + XD]));
                    dp[pc + q * 3 + XD] = sv;
                    dcv = _mm256_mul_ps(dcv, om.tfv[dd0 + q]);
                }
                if _mm256_movemask_ps(cv) == 0 { break; }
            }
        }
        for q in 0..q_ { xev = _mm256_add_ps(dp[pc + q * 3 + XD], xev); }
        xe = hsum(xev);
        xn *= xf.n[LOOP];
        xc = (xc * xf.c[LOOP]) + (xe * xf.e[MOVE]);
        xj = (xj * xf.j[LOOP]) + (xe * xf.e[LOOP]);
        xb = (xj * xf.j[MOVE]) + (xn * xf.n[MOVE]);
        let xr = i * NX;
        if xe > 1.0e4f32 {
            xn /= xe; xc /= xe; xj /= xe; xb /= xe;
            let inv = _mm256_set1_ps(1.0f32 / xe);
            for q in 0..q_ {
                dp[pc + q * 3 + XM] = _mm256_mul_ps(dp[pc + q * 3 + XM], inv);
                dp[pc + q * 3 + XD] = _mm256_mul_ps(dp[pc + q * 3 + XD], inv);
                dp[pc + q * 3 + XI] = _mm256_mul_ps(dp[pc + q * 3 + XI], inv);
            }
            x[xr + SSCALE] = xe;
            ox.totscale += (xe as f64).ln();
            xe = 1.0;
        } else {
            x[xr + SSCALE] = 1.0;
        }
        x[xr + SE] = xe; x[xr + SN] = xn; x[xr + SJ] = xj; x[xr + SB] = xb; x[xr + SC] = xc;
    }
    (ox.totscale + ((xc * xf.c[MOVE]) as f64).ln()) as f32
}

/// Forward's row times Backward's (M and I): all that decoding needs of the two.
#[target_feature(enable = "avx2")]
unsafe fn fold_row(f: &mut [__m256], b: &[__m256], q_: usize) {
    for q in 0..q_ {
        f[q * 3 + XM] = _mm256_mul_ps(f[q * 3 + XM], b[q * 3 + XM]);
        f[q * 3 + XI] = _mm256_mul_ps(f[q * 3 + XI], b[q * 3 + XI]);
    }
}

/// Backward keeps two dp rows (row i in slot i & 1); each finished row is folded into Forward's.
#[target_feature(enable = "avx2")]
unsafe fn backward(dsq: &[u8], om: &OProfile, xf: &Xf, fwd: &mut Omx, bck: &mut Omx) {
    let l = dsq.len() - 1;
    let q_ = om.q;
    let row = q_ * 3;
    let zerov = _mm256_setzero_ps();
    let t = &om.tfv;
    bck.has_own_scales = false;
    let dp = &mut bck.dp;
    let x = &mut bck.xmx;
    let (fx, fdp) = (&fwd.xmx, &mut fwd.dp);
    let mut xj = 0.0f32; let mut xb = 0.0f32; let mut xn = 0.0f32;
    let mut xc = xf.c[MOVE];
    let mut xe = xc * xf.e[MOVE];
    let mut xev = _mm256_set1_ps(xe);
    let mut dcv = zerov;
    let pc = (l & 1) * row;
    for q in 0..q_ { dp[pc + q * 3 + XM] = xev; dp[pc + q * 3 + XD] = xev; }
    for q in 0..q_ { dp[pc + q * 3 + XI] = zerov; }
    // L row DD paths
    let mut tp = 8 * q_ - 1;
    let mut dpv = leftshiftz(dp[pc + (q_ - 1) * 3 + XD]);
    for q in (0..q_).rev() {
        dcv = _mm256_mul_ps(dpv, t[tp]); tp = tp.wrapping_sub(1);
        dp[pc + q * 3 + XD] = _mm256_add_ps(dp[pc + q * 3 + XD], dcv);
        dpv = dp[pc + q * 3 + XD];
    }
    for _ in 1..4 {
        tp = 8 * q_ - 1;
        dcv = leftshiftz(dcv);
        for q in (0..q_).rev() {
            dcv = _mm256_mul_ps(dcv, t[tp]); tp = tp.wrapping_sub(1);
            dp[pc + q * 3 + XD] = _mm256_add_ps(dp[pc + q * 3 + XD], dcv);
        }
    }
    // MD init
    let mut tpi: isize = 7 * q_ as isize - 3;
    dcv = leftshiftz(dp[pc + XD]);
    for q in (0..q_).rev() {
        dp[pc + q * 3 + XM] = _mm256_add_ps(dp[pc + q * 3 + XM], _mm256_mul_ps(dcv, t[tpi as usize])); tpi -= 7;
        dcv = dp[pc + q * 3 + XD];
    }
    let lr = l * NX;
    if fx[lr + SSCALE] > 1.0 {
        let sc = fx[lr + SSCALE];
        xe /= sc; xn /= sc; xc /= sc; xj /= sc; xb /= sc;
        xev = _mm256_set1_ps(1.0f32 / sc);
        for q in 0..q_ {
            dp[pc + q * 3 + XM] = _mm256_mul_ps(dp[pc + q * 3 + XM], xev);
            dp[pc + q * 3 + XD] = _mm256_mul_ps(dp[pc + q * 3 + XD], xev);
            dp[pc + q * 3 + XI] = _mm256_mul_ps(dp[pc + q * 3 + XI], xev);
        }
    }
    fold_row(&mut fdp[l * row..(l + 1) * row], &dp[pc..pc + row], q_);
    x[lr + SSCALE] = fx[lr + SSCALE];
    bck.totscale = (x[lr + SSCALE] as f64).ln();
    x[lr + SE] = xe; x[lr + SN] = xn; x[lr + SJ] = xj; x[lr + SB] = xb; x[lr + SC] = xc;

    for i in (1..l).rev() {
        let pc = (i & 1) * row;
        let pp = ((i + 1) & 1) * row;
        let rbase = dsq[i + 1] as usize * q_;
        let mut rp = rbase + q_ - 1;
        let mut tp: isize = 7 * q_ as isize - 1;
        let mut tmmv = leftshiftz(t[1]);
        let mut timv = leftshiftz(t[2]);
        let mut tdmv = leftshiftz(t[3]);
        let mut mpv = _mm256_mul_ps(dp[pp + XM], om.rfv[rbase]);
        mpv = leftshiftz(mpv);
        let mut xbv = zerov;
        for q in (0..q_).rev() {
            let ipv = dp[pp + q * 3 + XI];
            dp[pc + q * 3 + XI] = _mm256_add_ps(_mm256_mul_ps(ipv, t[tp as usize]), _mm256_mul_ps(mpv, timv)); tp -= 1;
            dp[pc + q * 3 + XD] = _mm256_mul_ps(mpv, tdmv);
            let mcv = _mm256_add_ps(_mm256_mul_ps(ipv, t[tp as usize]), _mm256_mul_ps(mpv, tmmv)); tp -= 2;
            mpv = _mm256_mul_ps(dp[pp + q * 3 + XM], om.rfv[rp]); rp = rp.wrapping_sub(1);
            dp[pc + q * 3 + XM] = mcv;
            tdmv = t[tp as usize]; tp -= 1;
            timv = t[tp as usize]; tp -= 1;
            tmmv = t[tp as usize]; tp -= 1;
            xbv = _mm256_add_ps(xbv, _mm256_mul_ps(mpv, t[tp as usize])); tp -= 1;
        }
        xb = hsum(xbv);
        xc *= xf.c[LOOP];
        xj = (xb * xf.j[MOVE]) + (xj * xf.j[LOOP]);
        xn = (xb * xf.n[MOVE]) + (xn * xf.n[LOOP]);
        xe = (xc * xf.e[MOVE]) + (xj * xf.e[LOOP]);
        xev = _mm256_set1_ps(xe);
        let mut tp2 = 8 * q_ - 1;
        let mut dpv = _mm256_add_ps(dp[pc + XD], xev);
        dpv = leftshiftz(dpv);
        for q in (0..q_).rev() {
            dcv = _mm256_mul_ps(dpv, t[tp2]); tp2 = tp2.wrapping_sub(1);
            dp[pc + q * 3 + XD] = _mm256_add_ps(dp[pc + q * 3 + XD], _mm256_add_ps(dcv, xev));
            dpv = dp[pc + q * 3 + XD];
            dp[pc + q * 3 + XM] = _mm256_add_ps(dp[pc + q * 3 + XM], xev);
        }
        for _ in 1..4 {
            dcv = leftshiftz(dcv);
            let mut tp3 = 8 * q_ - 1;
            for q in (0..q_).rev() {
                dcv = _mm256_mul_ps(dcv, t[tp3]); tp3 = tp3.wrapping_sub(1);
                dp[pc + q * 3 + XD] = _mm256_add_ps(dp[pc + q * 3 + XD], dcv);
            }
        }
        dcv = leftshiftz(dp[pc + XD]);
        let mut tp4: isize = 7 * q_ as isize - 3;
        for q in (0..q_).rev() {
            dp[pc + q * 3 + XM] = _mm256_add_ps(dp[pc + q * 3 + XM], _mm256_mul_ps(dcv, t[tp4 as usize])); tp4 -= 7;
            dcv = dp[pc + q * 3 + XD];
        }
        if xb > 1.0e16f32 { bck.has_own_scales = true; }
        let xr = i * NX;
        x[xr + SSCALE] = if bck.has_own_scales { if xb > 1.0e4f32 { xb } else { 1.0 } } else { fx[xr + SSCALE] };
        if x[xr + SSCALE] > 1.0 {
            let sc = x[xr + SSCALE];
            xe /= sc; xn /= sc; xj /= sc; xb /= sc; xc /= sc;
            let inv = _mm256_set1_ps(1.0f32 / sc);
            for q in 0..q_ {
                dp[pc + q * 3 + XM] = _mm256_mul_ps(dp[pc + q * 3 + XM], inv);
                dp[pc + q * 3 + XD] = _mm256_mul_ps(dp[pc + q * 3 + XD], inv);
                dp[pc + q * 3 + XI] = _mm256_mul_ps(dp[pc + q * 3 + XI], inv);
            }
            bck.totscale += (sc as f64).ln();
        }
        fold_row(&mut fdp[i * row..(i + 1) * row], &dp[pc..pc + row], q_);
        x[xr + SE] = xe; x[xr + SN] = xn; x[xr + SJ] = xj; x[xr + SB] = xb; x[xr + SC] = xc;
    }
    // i = 0
    let pp = row;
    let rbase = dsq[1] as usize * q_;
    let mut xbv = zerov;
    for q in 0..q_ {
        let mut mpv = _mm256_mul_ps(dp[pp + q * 3 + XM], om.rfv[rbase + q]);
        mpv = _mm256_mul_ps(mpv, t[7 * q]);
        xbv = _mm256_add_ps(xbv, mpv);
    }
    xb = hsum(xbv);
    xn = (xb * xf.n[MOVE]) + (xn * xf.n[LOOP]);
    x[SB] = xb; x[SC] = 0.0; x[SJ] = 0.0; x[SN] = xn; x[SE] = 0.0; x[SSCALE] = 1.0;
}

/// Posterior decoding and the optimal accuracy fill, in place in fwd (Forward
/// times Backward), one row at a time. false on overflow.
#[target_feature(enable = "avx2")]
unsafe fn decode_optacc(l: usize, om: &OProfile, xf: &Xf, fwd: &mut Omx, bck: &mut Omx) -> bool {
    let q_ = om.q;
    let row = q_ * 3;
    let mut scaleproduct: f32 = (1.0f64 / bck.xmx[SN] as f64) as f32;
    let has_own = bck.has_own_scales;
    for s in 0..5 { bck.xmx[s] = 0.0; }
    // Forward's N, J, C of the previous row: the fill overwrites them
    let mut fprev = [fwd.xmx[SN], fwd.xmx[SJ], fwd.xmx[SC]];
    let zerov = _mm256_setzero_ps();
    let infv = _mm256_set1_ps(f32::NEG_INFINITY);
    let rs = |v: __m256| _mm256_blend_ps::<0x01>(rightshiftz(v), infv);
    let t = &om.tfv;
    let gt = |v: __m256| _mm256_cmp_ps::<_CMP_GT_OQ>(v, zerov);
    let (pp, ox) = (bck, fwd);
    for q in 0..q_ { ox.dp[q * 3 + XM] = infv; ox.dp[q * 3 + XI] = infv; ox.dp[q * 3 + XD] = infv; }
    ox.xmx[SE] = f32::NEG_INFINITY; ox.xmx[SN] = 0.0; ox.xmx[SJ] = f32::NEG_INFINITY; ox.xmx[SB] = 0.0; ox.xmx[SC] = f32::NEG_INFINITY;
    for i in 1..=l {
        let xr = i * NX;
        let pr = (i - 1) * NX;
        // decoding
        let totrv = _mm256_set1_ps(scaleproduct * ox.xmx[xr + SSCALE]);
        let base = i * row;
        for q in 0..q_ {
            let mi = base + q * 3;
            ox.dp[mi + XM] = _mm256_mul_ps(ox.dp[mi + XM], totrv);
            ox.dp[mi + XI] = _mm256_mul_ps(ox.dp[mi + XI], totrv);
        }
        pp.xmx[xr + SE] = 0.0;
        pp.xmx[xr + SN] = fprev[0] * pp.xmx[xr + SN] * xf.n[LOOP] * scaleproduct;
        pp.xmx[xr + SJ] = fprev[1] * pp.xmx[xr + SJ] * xf.j[LOOP] * scaleproduct;
        pp.xmx[xr + SC] = fprev[2] * pp.xmx[xr + SC] * xf.c[LOOP] * scaleproduct;
        pp.xmx[xr + SB] = 0.0;
        fprev = [ox.xmx[xr + SN], ox.xmx[xr + SJ], ox.xmx[xr + SC]];
        if has_own { scaleproduct *= ox.xmx[xr + SSCALE] / pp.xmx[xr + SSCALE]; }
        // optimal accuracy
        let pc = i * row;
        let pv = (i - 1) * row;
        let mut tp = 0usize;
        let mut dcv = infv;
        let mut xev = infv;
        let xbv = _mm256_set1_ps(ox.xmx[pr + SB]);
        let mut mpv = rs(ox.dp[pv + (q_ - 1) * 3 + XM]);
        let mut dpv = rs(ox.dp[pv + (q_ - 1) * 3 + XD]);
        let mut ipv = rs(ox.dp[pv + (q_ - 1) * 3 + XI]);
        for q in 0..q_ {
            let mut sv = _mm256_and_ps(gt(t[tp]), xbv); tp += 1;
            sv = _mm256_max_ps(sv, _mm256_and_ps(gt(t[tp]), mpv)); tp += 1;
            sv = _mm256_max_ps(sv, _mm256_and_ps(gt(t[tp]), ipv)); tp += 1;
            sv = _mm256_max_ps(sv, _mm256_and_ps(gt(t[tp]), dpv)); tp += 1;
            sv = _mm256_add_ps(sv, ox.dp[pc + q * 3 + XM]);
            xev = _mm256_max_ps(xev, sv);
            mpv = ox.dp[pv + q * 3 + XM];
            dpv = ox.dp[pv + q * 3 + XD];
            ipv = ox.dp[pv + q * 3 + XI];
            ox.dp[pc + q * 3 + XM] = sv;
            ox.dp[pc + q * 3 + XD] = dcv;
            dcv = _mm256_and_ps(gt(t[tp]), sv); tp += 1;
            let mut s2 = _mm256_and_ps(gt(t[tp]), mpv); tp += 1;
            s2 = _mm256_max_ps(s2, _mm256_and_ps(gt(t[tp]), ipv)); tp += 1;
            ox.dp[pc + q * 3 + XI] = _mm256_add_ps(s2, ox.dp[pc + q * 3 + XI]);
        }
        dcv = rs(dcv);
        let dd0 = 7 * q_;
        for q in 0..q_ {
            ox.dp[pc + q * 3 + XD] = _mm256_max_ps(dcv, ox.dp[pc + q * 3 + XD]);
            dcv = _mm256_and_ps(gt(t[dd0 + q]), ox.dp[pc + q * 3 + XD]);
        }
        for _ in 1..8 {
            dcv = rs(dcv);
            for q in 0..q_ {
                ox.dp[pc + q * 3 + XD] = _mm256_max_ps(dcv, ox.dp[pc + q * 3 + XD]);
                dcv = _mm256_and_ps(gt(t[dd0 + q]), dcv);
            }
        }
        for q in 0..q_ { xev = _mm256_max_ps(xev, ox.dp[pc + q * 3 + XD]); }
        ox.xmx[xr + SE] = hmax(xev);
        let x = &mut ox.xmx;
        let t1 = if xf.j[LOOP] == 0.0 { 0.0 } else { x[pr + SJ] + pp.xmx[xr + SJ] };
        let t2 = if xf.e[LOOP] == 0.0 { 0.0 } else { x[xr + SE] };
        x[xr + SJ] = t1.max(t2);
        let t1 = if xf.c[LOOP] == 0.0 { 0.0 } else { x[pr + SC] + pp.xmx[xr + SC] };
        let t2 = if xf.e[MOVE] == 0.0 { 0.0 } else { x[xr + SE] };
        x[xr + SC] = t1.max(t2);
        x[xr + SN] = if xf.n[LOOP] == 0.0 { 0.0 } else { x[pr + SN] + pp.xmx[xr + SN] };
        let t1 = if xf.n[MOVE] == 0.0 { 0.0 } else { x[xr + SN] };
        let t2 = if xf.j[MOVE] == 0.0 { 0.0 } else { x[xr + SJ] };
        x[xr + SB] = t1.max(t2);
    }
    !scaleproduct.is_infinite()
}

#[derive(Clone, Copy, PartialEq)]
enum St { M, D, I, N, C, J, E, B, S }

fn esl_max(a: f32, b: f32) -> f32 { if a > b { a } else { b } }

/// Optimal accuracy traceback: matched (k in slice, i) pairs in sequence order.
#[target_feature(enable = "avx2")]
unsafe fn oatrace(l: usize, om: &OProfile, xf: &Xf, pp: &Omx, ox: &Omx) -> Vec<(usize, usize)> {
    let q_ = om.q;
    let row = q_ * 3;
    let t = &om.tfv;
    let cell = |i: usize, q: usize, s: usize| ox.dp[i * row + q * 3 + s];
    let (mut i, mut k) = (l, 0usize);
    let mut s0 = St::C;
    let mut out = Vec::new();
    let _ = esl_max;
    let mut guard = 0usize;
    while s0 != St::S {
        guard += 1;
        if guard > 4 * (l + 2) * (om.m + 2) { break; }
        let s1 = match s0 {
            St::M => {
                let q = (k - 1) % q_; let r = (k - 1) / q_;
                let tb = 7 * q;
                let xbv = _mm256_set1_ps(ox.xmx[(i - 1) * NX + SB]);
                let (mpv, dpv, ipv) = if q > 0 {
                    (cell(i - 1, q - 1, XM), cell(i - 1, q - 1, XD), cell(i - 1, q - 1, XI))
                } else {
                    (rightshiftz(cell(i - 1, q_ - 1, XM)), rightshiftz(cell(i - 1, q_ - 1, XD)), rightshiftz(cell(i - 1, q_ - 1, XI)))
                };
                let pick = |v: __m256, tv: __m256| { let tl = lanes(tv)[r]; if tl == 0.0 { f32::NEG_INFINITY } else { lanes(v)[r] } };
                let path = [pick(mpv, t[tb + 1]), pick(ipv, t[tb + 2]), pick(dpv, t[tb + 3]), pick(xbv, t[tb])];
                let mut a = 0;
                for j in 1..4 { if path[j] > path[a] { a = j; } }
                k -= 1; i -= 1;
                [St::M, St::I, St::D, St::B][a]
            }
            St::D => {
                let q = (k - 1) % q_; let r = (k - 1) / q_;
                let (mpv, dpv, tmd, tdd) = if q > 0 {
                    (cell(i, q - 1, XM), cell(i, q - 1, XD), t[7 * (q - 1) + 4], t[7 * q_ + (q - 1)])
                } else {
                    (rightshiftz(cell(i, q_ - 1, XM)), rightshiftz(cell(i, q_ - 1, XD)), rightshiftz(t[7 * (q_ - 1) + 4]), rightshiftz(t[8 * q_ - 1]))
                };
                let p0 = if lanes(tmd)[r] == 0.0 { f32::NEG_INFINITY } else { lanes(mpv)[r] };
                let p1 = if lanes(tdd)[r] == 0.0 { f32::NEG_INFINITY } else { lanes(dpv)[r] };
                k -= 1;
                if p0 >= p1 { St::M } else { St::D }
            }
            St::I => {
                let q = (k - 1) % q_; let r = (k - 1) / q_;
                let tb = 7 * q + 5;
                let p0 = if lanes(t[tb])[r] == 0.0 { f32::NEG_INFINITY } else { lanes(cell(i - 1, q, XM))[r] };
                let p1 = if lanes(t[tb + 1])[r] == 0.0 { f32::NEG_INFINITY } else { lanes(cell(i - 1, q, XI))[r] };
                i -= 1;
                if p0 >= p1 { St::M } else { St::I }
            }
            St::N => if i == 0 { St::S } else { St::N },
            St::C => {
                let p0 = if xf.c[LOOP] == 0.0 { f32::NEG_INFINITY } else { ox.xmx[(i - 1) * NX + SC] + pp.xmx[i * NX + SC] };
                let p1 = if xf.e[MOVE] == 0.0 { f32::NEG_INFINITY } else { ox.xmx[i * NX + SE] };
                if p0 > p1 { St::C } else { St::E }
            }
            St::J => {
                let p0 = if xf.j[LOOP] == 0.0 { f32::NEG_INFINITY } else { ox.xmx[(i - 1) * NX + SJ] + pp.xmx[i * NX + SJ] };
                let p1 = if xf.e[LOOP] == 0.0 { f32::NEG_INFINITY } else { ox.xmx[i * NX + SE] };
                if p0 > p1 { St::J } else { St::E }
            }
            St::E => {
                let (mut mx, mut smax, mut kmax) = (f32::NEG_INFINITY, St::M, 1usize);
                for q in 0..q_ {
                    let mv = lanes(cell(i, q, XM));
                    for r in 0..8 { if mv[r] >= mx { mx = mv[r]; smax = St::M; kmax = r * q_ + q + 1; } }
                    let dv = lanes(cell(i, q, XD));
                    for r in 0..8 { if dv[r] > mx { mx = dv[r]; smax = St::D; kmax = r * q_ + q + 1; } }
                }
                k = kmax;
                smax
            }
            St::B => {
                let p0 = if xf.n[MOVE] == 0.0 { f32::NEG_INFINITY } else { ox.xmx[i * NX + SN] };
                let p1 = if xf.j[MOVE] == 0.0 { f32::NEG_INFINITY } else { ox.xmx[i * NX + SJ] };
                if p0 > p1 { St::N } else { St::J }
            }
            St::S => St::S,
        };
        if s1 == St::M { out.push((k, i)); }
        if (s1 == St::N || s1 == St::J || s1 == St::C) && s1 == s0 { i -= 1; }
        s0 = s1;
    }
    out.reverse();
    out
}

/// One unilocal alignment, AVX2. `None` when decoding overflows (the caller
/// falls back to the generic path, as BATH does).
pub fn align(p: &Profile, om: &OProfile, seq: &[u8], dp: &mut AvxDp) -> Option<(AlnHit, Vec<(usize, usize)>)> {
    let l = seq.len();
    if l == 0 { return Some((AlnHit::default(), Vec::new())); }
    let mut dsq = Vec::with_capacity(l + 1);
    dsq.push(0u8);
    dsq.extend_from_slice(seq);
    let xf = xf_for(om, l);
    dp.fwd.grow(l + 1, om.q);
    dp.bck.grow2(2, l + 1, om.q);
    unsafe {
        let fwdsc = forward::<false>(&dsq, om, &xf, &mut dp.fwd);
        backward(&dsq, om, &xf, &mut dp.fwd, &mut dp.bck);
        let p1 = l as f32 / (l as f32 + 1.0);
        let nullsc = ((l as f32) as f64 * (p1 as f64).ln() + (1.0 - p1 as f64).ln()) as f32;
        let bits = ((fwdsc - nullsc) as f64 / std::f64::consts::LN_2) as f32;
        let bcksc = (dp.bck.totscale + (dp.bck.xmx[SN] as f64).ln()) as f32;
        dp.bck_bits = ((bcksc - nullsc) as f64 / std::f64::consts::LN_2) as f32;
        if !decode_optacc(l, om, &xf, &mut dp.fwd, &mut dp.bck) { return None; }
        let tr = oatrace(l, om, &xf, &dp.bck, &dp.fwd);
        let mut h = AlnHit { bits, ..Default::default() };
        let trace: Vec<(usize, usize)> = tr.iter().map(|&(kk, ii)| (kk + p.a - 1, ii)).collect();
        for (n, &(node, ii)) in trace.iter().enumerate() {
            if n == 0 { h.kf = node as i32; h.rf = ii as i32; }
            h.kl = node as i32; h.rl = ii as i32;
            h.nmatch += 1;
        }
        Some((h, trace))
    }
}

/// ln C(i) of `seq` per residue (the C state after residue i, unscaled;
/// entry 0 is 0), from Forward alone.
pub fn fwd_lnc(om: &OProfile, seq: &[u8], dp: &mut AvxDp) -> Vec<f64> {
    let l = seq.len();
    let mut lnc = vec![0.0f64; l + 1];
    if l == 0 { return lnc; }
    let mut dsq = Vec::with_capacity(l + 1);
    dsq.push(0u8);
    dsq.extend_from_slice(seq);
    let xf = xf_for(om, l);
    dp.fwd.grow2(2, l + 1, om.q);
    unsafe { forward::<true>(&dsq, om, &xf, &mut dp.fwd); }
    let mut cum = 0.0f64;
    for i in 1..=l {
        cum += (dp.fwd.xmx[i * NX + SSCALE] as f64).ln();
        lnc[i] = cum + (dp.fwd.xmx[i * NX + SC] as f64).ln();
    }
    lnc
}

/// Forward bits of `seq` (null-corrected), without Backward or decoding.
pub fn fwd_bits(om: &OProfile, seq: &[u8], dp: &mut AvxDp) -> f32 {
    let l = seq.len();
    if l == 0 { return 0.0; }
    let mut dsq = Vec::with_capacity(l + 1);
    dsq.push(0u8);
    dsq.extend_from_slice(seq);
    let xf = xf_for(om, l);
    dp.fwd.grow2(2, l + 1, om.q);
    unsafe {
        let fwdsc = forward::<true>(&dsq, om, &xf, &mut dp.fwd);
        let p1 = l as f32 / (l as f32 + 1.0);
        let nullsc = ((l as f32) as f64 * (p1 as f64).ln() + (1.0 - p1 as f64).ln()) as f32;
        ((fwdsc - nullsc) as f64 / std::f64::consts::LN_2) as f32
    }
}

// ---- Forward with an exponent in every cell: a cell holds m * 2^(64 e).
// The rescaled Forward keeps one scale per row and rounds a far-behind path to
// zero; here none is lost, so the score of a whole chain is exact.

const XEMIN: i32 = i32::MIN / 2;
const X2P64: f32 = 18446744073709551616.0;
const X2M64: f32 = 1.0 / X2P64;
const XLN: f64 = 64.0 * std::f64::consts::LN_2;

#[derive(Clone, Copy)]
struct Xv { m: __m256, e: __m256i }

/// A row of that Forward and its special states, to resume from.
#[derive(Clone)]
pub struct XRow { m: Vec<__m256>, e: Vec<__m256i>, xn: f64, xb: f64, xc: f64 }

/// 1 where the exponents agree, 2^-64 one block below, else 0.
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn xfactor(d: __m256i) -> __m256 {
    let is0 = _mm256_castsi256_ps(_mm256_cmpeq_epi32(d, _mm256_setzero_si256()));
    let is1 = _mm256_castsi256_ps(_mm256_cmpeq_epi32(d, _mm256_set1_epi32(1)));
    _mm256_or_ps(_mm256_and_ps(is0, _mm256_set1_ps(1.0)), _mm256_and_ps(is1, _mm256_set1_ps(X2M64)))
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn xadd(a: Xv, b: Xv) -> Xv {
    if _mm256_movemask_epi8(_mm256_cmpeq_epi32(a.e, b.e)) == -1 {
        return Xv { m: _mm256_add_ps(a.m, b.m), e: a.e };
    }
    // a zero has no exponent of its own
    let (z, emin) = (_mm256_setzero_ps(), _mm256_set1_epi32(XEMIN));
    let ea = _mm256_blendv_epi8(a.e, emin, _mm256_castps_si256(_mm256_cmp_ps::<_CMP_EQ_OQ>(a.m, z)));
    let eb = _mm256_blendv_epi8(b.e, emin, _mm256_castps_si256(_mm256_cmp_ps::<_CMP_EQ_OQ>(b.m, z)));
    let e = _mm256_max_epi32(ea, eb);
    Xv { m: _mm256_add_ps(_mm256_mul_ps(a.m, xfactor(_mm256_sub_epi32(e, ea))), _mm256_mul_ps(b.m, xfactor(_mm256_sub_epi32(e, eb)))), e }
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn xmul(a: Xv, t: __m256) -> Xv { Xv { m: _mm256_mul_ps(a.m, t), e: a.e } }

/// Mantissas back into [2^-32, 2^32).
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn xnorm(v: Xv) -> Xv {
    let big = _mm256_cmp_ps::<_CMP_GE_OQ>(v.m, _mm256_set1_ps(4294967296.0));
    let small = _mm256_and_ps(_mm256_cmp_ps::<_CMP_LT_OQ>(v.m, _mm256_set1_ps(1.0 / 4294967296.0)),
                              _mm256_cmp_ps::<_CMP_GT_OQ>(v.m, _mm256_setzero_ps()));
    if _mm256_movemask_ps(_mm256_or_ps(big, small)) == 0 { return v; }
    let m = _mm256_blendv_ps(_mm256_blendv_ps(v.m, _mm256_mul_ps(v.m, _mm256_set1_ps(X2M64)), big),
                             _mm256_mul_ps(v.m, _mm256_set1_ps(X2P64)), small);
    // a set lane is -1 as an integer: +1 where big, -1 where small
    Xv { m, e: _mm256_add_epi32(_mm256_sub_epi32(v.e, _mm256_castps_si256(big)), _mm256_castps_si256(small)) }
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn xshift(v: Xv) -> Xv {
    Xv { m: rightshiftz(v.m), e: _mm256_castps_si256(rightshiftz(_mm256_castsi256_ps(v.e))) }
}

/// ln of the sum of the lanes.
#[target_feature(enable = "avx2")]
unsafe fn xhsum_ln(v: Xv) -> f64 {
    let m: [f32; 8] = std::mem::transmute(v.m);
    let e: [i32; 8] = std::mem::transmute(v.e);
    let top = (0..8).filter(|&i| m[i] > 0.0).map(|i| e[i]).max();
    let Some(top) = top else { return f64::NEG_INFINITY };
    let mut s = 0.0f64;
    for i in 0..8 {
        if m[i] > 0.0 && top - e[i] <= 2 { s += m[i] as f64 * (2.0f64).powi(-64 * (top - e[i])); }
    }
    s.ln() + top as f64 * XLN
}

fn lse(a: f64, b: f64) -> f64 {
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    if lo == f64::NEG_INFINITY { hi } else { hi + (lo - hi).exp().ln_1p() }
}

/// Rows after `start` of the Forward of dsq[1..=l], resumed from `init` (None: row 0).
/// The state after each row of `marks` (ascending) goes to `snaps`. Score in nats.
#[target_feature(enable = "avx2")]
unsafe fn xforward(dsq: &[u8], om: &OProfile, start: usize, init: Option<&XRow>, marks: &[usize], snaps: &mut Vec<XRow>) -> f32 {
    let l = dsq.len() - 1;
    let q_ = om.q;
    let row = q_ * 3;
    let xf = xf_for(om, l);
    debug_assert!(xf.e[LOOP] == 0.0); // unilocal: J is never entered
    let (nloop, nmove) = ((xf.n[LOOP] as f64).ln(), (xf.n[MOVE] as f64).ln());
    let (cloop, cmove, emove) = ((xf.c[LOOP] as f64).ln(), (xf.c[MOVE] as f64).ln(), (xf.e[MOVE] as f64).ln());
    let zx = Xv { m: _mm256_setzero_ps(), e: _mm256_set1_epi32(XEMIN) };
    let mut dm = vec![zx.m; 2 * row];
    let mut de = vec![zx.e; 2 * row];
    let (mut xn, mut xb, mut xc) = (0.0f64, nmove, f64::NEG_INFINITY);
    if let Some(s) = init {
        let pc = (start & 1) * row;
        dm[pc..pc + row].copy_from_slice(&s.m);
        de[pc..pc + row].copy_from_slice(&s.e);
        xn = s.xn; xb = s.xb; xc = s.xc;
    }
    let mut mk = 0usize;
    let t = &om.tfv;
    let dd0 = 7 * q_;
    for i in start..=l {
        if i > start {
            let (pc, pp) = ((i & 1) * row, ((i - 1) & 1) * row);
            let rp = dsq[i] as usize * q_;
            let xbv = if xb == f64::NEG_INFINITY { zx } else {
                let e = (xb / XLN).round();
                Xv { m: _mm256_set1_ps((xb - e * XLN).exp() as f32), e: _mm256_set1_epi32(e as i32) }
            };
            let at = |s: usize| Xv { m: dm[s], e: de[s] };
            let mut mpv = xshift(at(pp + (q_ - 1) * 3 + XM));
            let mut dpv = xshift(at(pp + (q_ - 1) * 3 + XD));
            let mut ipv = xshift(at(pp + (q_ - 1) * 3 + XI));
            let (mut dcv, mut xev) = (zx, zx);
            let mut tp = 0usize;
            for q in 0..q_ {
                let mut sv = xadd(xadd(xadd(xmul(xbv, t[tp]), xmul(mpv, t[tp + 1])), xmul(ipv, t[tp + 2])), xmul(dpv, t[tp + 3]));
                sv = xnorm(xmul(sv, om.rfv[rp + q]));
                xev = xadd(xev, sv);
                let c = pc + q * 3;
                mpv = Xv { m: dm[pp + q * 3 + XM], e: de[pp + q * 3 + XM] };
                dpv = Xv { m: dm[pp + q * 3 + XD], e: de[pp + q * 3 + XD] };
                ipv = Xv { m: dm[pp + q * 3 + XI], e: de[pp + q * 3 + XI] };
                dm[c + XM] = sv.m; de[c + XM] = sv.e;
                dm[c + XD] = dcv.m; de[c + XD] = dcv.e;
                dcv = xnorm(xmul(sv, t[tp + 4]));
                let iv = xnorm(xadd(xmul(mpv, t[tp + 5]), xmul(ipv, t[tp + 6])));
                dm[c + XI] = iv.m; de[c + XI] = iv.e;
                tp += 7;
            }
            // D paths: along the nodes, then around the stripes until nothing changes
            dcv = xshift(dcv);
            dm[pc + XD] = zx.m; de[pc + XD] = zx.e;
            for q in 0..q_ {
                let c = pc + q * 3 + XD;
                let d = xnorm(xadd(dcv, Xv { m: dm[c], e: de[c] }));
                dm[c] = d.m; de[c] = d.e;
                dcv = xnorm(xmul(d, t[dd0 + q]));
            }
            for _ in 1..8 {
                dcv = xshift(dcv);
                let mut changed = 0i32;
                for q in 0..q_ {
                    let c = pc + q * 3 + XD;
                    let d = xnorm(xadd(dcv, Xv { m: dm[c], e: de[c] }));
                    changed |= _mm256_movemask_ps(_mm256_cmp_ps::<_CMP_NEQ_UQ>(d.m, dm[c])) | !_mm256_movemask_epi8(_mm256_cmpeq_epi32(d.e, de[c]));
                    dm[c] = d.m; de[c] = d.e;
                    dcv = xnorm(xmul(dcv, t[dd0 + q]));
                }
                if changed == 0 { break; }
            }
            for q in 0..q_ { xev = xadd(xev, Xv { m: dm[pc + q * 3 + XD], e: de[pc + q * 3 + XD] }); }
            let xe = xhsum_ln(xev);
            xn += nloop;
            xc = lse(xc + cloop, xe + emove);
            xb = xn + nmove;
        }
        while mk < marks.len() && marks[mk] == i {
            let pc = (i & 1) * row;
            snaps.push(XRow { m: dm[pc..pc + row].to_vec(), e: de[pc..pc + row].to_vec(), xn, xb, xc });
            mk += 1;
        }
    }
    (xc + cmove) as f32
}

/// Forward bits of `seq`, exact for a whole chain. `init` and `start` resume a call
/// on a sequence with the same first `start` residues and the same length.
pub fn fwd_bits_x(om: &OProfile, seq: &[u8], start: usize, init: Option<&XRow>, marks: &[usize], snaps: &mut Vec<XRow>) -> f32 {
    let l = seq.len();
    if l == 0 { return 0.0; }
    let mut dsq = Vec::with_capacity(l + 1);
    dsq.push(0u8);
    dsq.extend_from_slice(seq);
    unsafe {
        let fwdsc = xforward(&dsq, om, start, init, marks, snaps);
        let p1 = l as f32 / (l as f32 + 1.0);
        let nullsc = ((l as f32) as f64 * (p1 as f64).ln() + (1.0 - p1 as f64).ln()) as f32;
        ((fwdsc - nullsc) as f64 / std::f64::consts::LN_2) as f32
    }
}
