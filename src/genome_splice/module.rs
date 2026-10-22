//! Mutually exclusive alternatives of an internal exon: copies of it in the
//! introns on either side. Each intron is translated in three frames with stops
//! read as X, so a degraded copy aligns as one piece; its stops are left to the
//! pseudogene stage, not held against it. A hit is an alternative by
//! the module rule (similar length, the same model nodes) when it
//! beats its own reversal clearly and per residue: repeats give long hits
//! whose score is spread thin. The same scan finds the copies of an
//! N-terminal row's tail for the rebase step. seq_copy finds a sibling row's exon
//! again in a row's gap by the exon's own residues (BLOSUM62), for sibling_fill.

use super::chain::{Chain, ChainExon, SRC_ALT, SRC_SIBFILL};
use super::hmm::Hmm;
use super::orf::{is_repeat, translate_frame, Aligner};
use super::sites::{revcomp, Code};
use crate::blosum_tables::{BLOSUM, INVALID, XLAT};

/// node overlap of the smaller span (the length ratio and the overlap of the
/// larger span are the caller's: on the hit, or later on the refined exons)
const NODES_SMALL: f64 = 0.80;
/// bits a stretch's best alignment must score to be searched further
const MIN_BITS: f32 = 5.0;
/// bits a hit must score above its reversal, in all and per residue
const MIN_MARGIN: f32 = 9.0;
const MIN_DENSITY: f32 = 0.2;
/// a short hit may pass on less, if denser
const SHORT_MARGIN: f32 = 6.0;
const SHORT_DENSITY: f32 = 0.25;
/// nt between an alternative and the exon it replaces
const MIN_DIST: i64 = 30;
/// nodes past the exon's own range an alternative may reach
const MARGIN: i64 = 5;
/// hits tried per frame of one intron
const MAX_HITS: usize = 256;

fn oriented(sc: &[u8], lo: i64, hi: i64, strand: u8) -> Vec<u8> {
    let s = &sc[(lo - 1) as usize..hi as usize];
    if strand == b'-' { revcomp(s) } else { s.to_vec() }
}

/// Residues of exon x in its best frame (aligner set to its nodes).
fn exon_len(al: &mut Aligner, code: Code, sc: &[u8], x: &ChainExon, strand: u8) -> usize {
    let xnt = oriented(sc, x.start, x.end, strand);
    let (mut xb, mut xlen) = (f32::NEG_INFINITY, 0usize);
    for f in 0..3 {
        let t: Vec<u8> = translate_frame(code, &xnt, f).into_iter().map(|q| if q == b'*' { b'X' } else { q }).collect();
        if t.is_empty() { continue; }
        let bits = al.fwd_bits(&t);
        if bits > xb { xb = bits; xlen = t.len(); }
    }
    xlen
}

/// Copies of exon x (xlen residues) in lo..hi, each with its bits over its reversal (aligner set to x's nodes).
#[allow(clippy::too_many_arguments)]
fn scan(al: &mut Aligner, hmm: &Hmm, code: Code, strand: u8, x: &ChainExon, xlen: usize, sc: &[u8], lo: i64, hi: i64, size: f64, nodes_large: f64,
        out: &mut Vec<(ChainExon, f32)>) {
    if hi - lo + 1 < 3 * (size * xlen as f64) as i64 || hi as usize > sc.len() || lo < 1 { return; }
    let xspan = x.k2 - x.k1 + 1;
    let s = oriented(sc, lo, hi, strand);
    for f in 0..3usize {
        let raw = translate_frame(code, &s, f);
        let aa: Vec<u8> = raw.iter().map(|&q| if q == b'*' { b'X' } else { q }).collect();
        let mut todo = vec![(0usize, aa.len())];
        let mut tried = 0;
        while let Some((s0, e0)) = todo.pop() {
            if tried >= MAX_HITS || ((e0 - s0) as f64) < size * xlen as f64 { continue; }
            tried += 1;
            // Forward alone settles most stretches (the bits a full run reports)
            if al.fwd_bits(&aa[s0..e0]) < MIN_BITS { continue; }
            let h = al.one(&aa[s0..e0]);
            if h.nmatch == 0 || h.bits < MIN_BITS { continue; }
            let (hs, he) = (s0 + h.rf as usize - 1, s0 + h.rl as usize - 1);
            if hs > s0 { todo.push((s0, hs)); }
            if he + 1 < e0 { todo.push((he + 1, e0)); }
            let hlen = he + 1 - hs;
            let (kf, kl) = (h.kf as i64, h.kl as i64);
            let span = kl - kf + 1;
            let ov = kl.min(x.k2) - kf.max(x.k1) + 1;
            let ratio = hlen.min(xlen) as f64 / hlen.max(xlen) as f64;
            if ratio < size || (ov as f64) < NODES_SMALL * span.min(xspan) as f64 || (ov as f64) < nodes_large * span.max(xspan) as f64 {
                continue;
            }
            let hit = &aa[hs..=he];
            let rev: Vec<u8> = hit.iter().rev().copied().collect();
            let margin = al.fwd_bits(hit) - al.fwd_bits(&rev);
            let dens = margin / hlen as f32;
            if dens < MIN_DENSITY || (margin < MIN_MARGIN && (margin < SHORT_MARGIN || dens < SHORT_DENSITY)) { continue; }
            if is_repeat(hit, hmm, kf, kl) { continue; }
            let (o1, o2) = (f as i64 + 3 * hs as i64, f as i64 + 3 * he as i64 + 2);
            let (gs, ge) = if strand == b'-' { (hi - o2, hi - o1) } else { (lo + o1, lo + o2) };
            let dist = if ge < x.start { x.start - ge - 1 } else { gs - x.end - 1 };
            if dist < MIN_DIST { continue; }
            out.push((ChainExon { start: gs, end: ge, k1: kf, k2: kl, src: SRC_ALT, codon: 0, bits: 0.0, run: false, owes: false }, margin));
        }
    }
}

/// Alternatives of exon j (0 < j < last) of chain c.
/// size: least length ratio of hit to exon; nodes_large: least node overlap of the larger span.
/// stop: other modules' block in the intron before and after the exon; each scan ends there.
#[allow(clippy::too_many_arguments)]
pub fn alternatives(al: &mut Aligner, hid: usize, hmm: &Hmm, code: Code, c: &Chain, j: usize, sc: &[u8], size: f64, nodes_large: f64,
                    stop: [Option<(i64, i64)>; 2]) -> Vec<ChainExon> {
    let m = hmm.m as i64;
    let x = c.ex[j];
    let (a, b) = (1.max(x.k1 - MARGIN), m.min(x.k2 + MARGIN));
    if x.start < 1 || x.end as usize > sc.len() || b < a { return Vec::new(); }
    al.set(hid, hmm, a as usize, b as usize);
    let xlen = exon_len(al, code, sc, &x, c.strand);
    if xlen == 0 { return Vec::new(); }
    let mut out = Vec::new();
    for (nb, stop) in [(j - 1, stop[0]), (j + 1, stop[1])] {
        let p = c.ex[nb];
        let (mut lo, mut hi) = if p.end < x.start { (p.end + 1, x.start - 1) } else { (x.end + 1, p.start - 1) };
        if let Some((s0, s1)) = stop {
            if p.end < x.start { lo = lo.max(s1 + 1); } else { hi = hi.min(s0 - 1); }
        }
        scan(al, hmm, code, c.strand, &x, xlen, sc, lo, hi, size, nodes_large, &mut out);
    }
    out.into_iter().map(|(e, _)| e).collect()
}

/// Copies of exon x anywhere in lo..hi, each with its bits over its reversal.
#[allow(clippy::too_many_arguments)]
pub fn copies(al: &mut Aligner, hid: usize, hmm: &Hmm, code: Code, strand: u8, x: &ChainExon, sc: &[u8], lo: i64, hi: i64,
              size: f64, nodes_large: f64) -> Vec<(ChainExon, f32)> {
    let (a, b) = (1.max(x.k1 - MARGIN), (hmm.m as i64).min(x.k2 + MARGIN));
    let mut out = Vec::new();
    if x.start < 1 || x.end as usize > sc.len() || b < a { return out; }
    al.set(hid, hmm, a as usize, b as usize);
    let xlen = exon_len(al, code, sc, x, strand);
    if xlen > 0 { scan(al, hmm, code, strand, x, xlen, sc, lo, hi, size, nodes_large, &mut out); }
    out
}

/// BLOSUM62 a sequence copy must score over its reversal, in all and per residue
const SEQ_MIN_MARGIN: i32 = 60;
const SEQ_MIN_DENSITY: f32 = 0.8;
const SEQ_GAP_OPEN: i32 = 10;
const SEQ_GAP_EXT: i32 = 1;

fn b62(a: u8, b: u8) -> i32 {
    let (i, j) = (XLAT[a as usize], XLAT[b as usize]);
    if i == INVALID || j == INVALID { -4 } else { BLOSUM[27 * i as usize + j as usize] }
}

/// Best local alignment of q in t (affine gaps): score and t's first and last residue.
fn local(q: &[u8], t: &[u8]) -> (i32, usize, usize) {
    let n = q.len();
    let mut h = vec![0i32; n + 1];
    let mut hs = vec![0usize; n + 1];
    let mut e = vec![i32::MIN / 2; n + 1];
    let mut es = vec![0usize; n + 1];
    let mut best = (0, 0, 0);
    for (j, &tc) in t.iter().enumerate() {
        let (mut diag, mut diag_s) = (0, j);
        let (mut f, mut fs) = (i32::MIN / 2, j);
        for i in 1..=n {
            let (up, up_s) = (h[i], hs[i]);
            // gap in q (t advances): from the cell left of this one, kept in e
            if up - SEQ_GAP_OPEN > e[i] - SEQ_GAP_EXT { e[i] = up - SEQ_GAP_OPEN; es[i] = up_s; } else { e[i] -= SEQ_GAP_EXT; }
            // gap in t (q advances): from the cell above
            if h[i - 1] - SEQ_GAP_OPEN > f - SEQ_GAP_EXT { f = h[i - 1] - SEQ_GAP_OPEN; fs = hs[i - 1]; } else { f -= SEQ_GAP_EXT; }
            let (mut v, mut vs) = (diag + b62(q[i - 1], tc), if diag > 0 { diag_s } else { j });
            if e[i] > v { v = e[i]; vs = es[i]; }
            if f > v { v = f; vs = fs; }
            if v <= 0 { v = 0; vs = j + 1; }
            diag = up; diag_s = up_s;
            h[i] = v; hs[i] = vs;
            if v > best.0 { best = (v, vs, j); }
        }
    }
    best
}

/// The copy of exon x in lo..hi (spans in mask read as N) found by x's own residues: its best
/// local BLOSUM62 hit in three frames, kept when it beats its reversal clearly and per residue.
#[allow(clippy::too_many_arguments)]
pub fn seq_copy(code: Code, strand: u8, x: &ChainExon, sc: &[u8], lo: i64, hi: i64, mask: &[(i64, i64)], size: f64) -> Option<(ChainExon, f32)> {
    if x.start < 1 || x.end as usize > sc.len() || lo < 1 || hi as usize > sc.len() || hi - lo + 1 < 60 { return None; }
    let xnt = oriented(sc, x.start, x.end, strand);
    // x's frame: fewest stops, then longest
    let xp = (0..3).map(|f| translate_frame(code, &xnt, f)).min_by_key(|t| (t.iter().filter(|&&c| c == b'*').count(), usize::MAX - t.len()))?;
    let mut nt = sc[(lo - 1) as usize..hi as usize].to_vec();
    for &(s, e) in mask {
        for p in s.max(lo)..=e.min(hi) { nt[(p - lo) as usize] = b'N'; }
    }
    let s = if strand == b'-' { revcomp(&nt) } else { nt };
    let mut best: Option<(i32, usize, usize, usize, Vec<u8>)> = None;
    for f in 0..3usize {
        let t = translate_frame(code, &s, f);
        let (v, ts, te) = local(&xp, &t);
        if v > 0 && best.as_ref().map_or(true, |b| v > b.0) { best = Some((v, f, ts, te, t)); }
    }
    let (v, f, ts, te, t) = best?;
    let hit = &t[ts..=te];
    let ratio = hit.len().min(xp.len()) as f64 / hit.len().max(xp.len()) as f64;
    let rev: Vec<u8> = hit.iter().rev().copied().collect();
    let margin = v - local(&xp, &rev).0;
    if ratio < size || margin < SEQ_MIN_MARGIN || (margin as f32) / (hit.len() as f32) < SEQ_MIN_DENSITY { return None; }
    let (o1, o2) = (f as i64 + 3 * ts as i64, f as i64 + 3 * te as i64 + 2);
    let (gs, ge) = if strand == b'-' { (hi - o2, hi - o1) } else { (lo + o1, lo + o2) };
    Some((ChainExon { start: gs, end: ge, k1: x.k1, k2: x.k2, src: SRC_SIBFILL, codon: 0, bits: 0.0, run: false, owes: false }, margin as f32))
}

#[cfg(test)]
#[path = "../tests/genome_module.rs"]
mod tests;
