//! Mutually exclusive alternatives of an internal exon: copies of it in the
//! introns on either side. Each intron is translated in three frames with stops
//! read as X, so a degraded copy aligns as one piece; its stops are left to the
//! pseudogene stage, not held against it. A hit is an alternative by
//! exonfinder's module rule (similar length, the same model nodes) when it
//! beats its own reversal clearly and per residue: repeats give long hits
//! whose score is spread thin.

use super::chain::{Chain, ChainExon, SRC_ALT};
use super::hmm::Hmm;
use super::orf::{translate_frame, Aligner};
use super::sites::revcomp;

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

/// Alternatives of exon j (0 < j < last) of chain c.
/// size: least length ratio of hit to exon; nodes_large: least node overlap of the larger span.
/// stop: other modules' block in the intron before and after the exon; each scan ends there.
#[allow(clippy::too_many_arguments)]
pub fn alternatives(al: &mut Aligner, hid: usize, hmm: &Hmm, c: &Chain, j: usize, sc: &[u8], size: f64, nodes_large: f64,
                    stop: [Option<(i64, i64)>; 2]) -> Vec<ChainExon> {
    let m = hmm.m as i64;
    let x = c.ex[j];
    let (a, b) = (1.max(x.k1 - MARGIN), m.min(x.k2 + MARGIN));
    if x.start < 1 || x.end as usize > sc.len() || b < a { return Vec::new(); }
    al.set(hid, hmm, a as usize, b as usize);
    // the exon itself: its length in its best frame
    let xnt = oriented(sc, x.start, x.end, c.strand);
    let (mut xb, mut xlen) = (f32::NEG_INFINITY, 0usize);
    for f in 0..3 {
        let t: Vec<u8> = translate_frame(&xnt, f).into_iter().map(|q| if q == b'*' { b'X' } else { q }).collect();
        if t.is_empty() { continue; }
        let bits = al.fwd_bits(&t);
        if bits > xb { xb = bits; xlen = t.len(); }
    }
    if xlen == 0 { return Vec::new(); }
    let xspan = x.k2 - x.k1 + 1;
    let mut out = Vec::new();
    for (nb, stop) in [(j - 1, stop[0]), (j + 1, stop[1])] {
        let p = c.ex[nb];
        let (mut lo, mut hi) = if p.end < x.start { (p.end + 1, x.start - 1) } else { (x.end + 1, p.start - 1) };
        if let Some((s0, s1)) = stop {
            if p.end < x.start { lo = lo.max(s1 + 1); } else { hi = hi.min(s0 - 1); }
        }
        if hi - lo + 1 < 3 * (size * xlen as f64) as i64 || hi as usize > sc.len() || lo < 1 { continue; }
        let s = oriented(sc, lo, hi, c.strand);
        for f in 0..3usize {
            let raw = translate_frame(&s, f);
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
                let (o1, o2) = (f as i64 + 3 * hs as i64, f as i64 + 3 * he as i64 + 2);
                let (gs, ge) = if c.strand == b'-' { (hi - o2, hi - o1) } else { (lo + o1, lo + o2) };
                let dist = if ge < x.start { x.start - ge - 1 } else { gs - x.end - 1 };
                if dist < MIN_DIST { continue; }
                out.push(ChainExon { start: gs, end: ge, k1: kf, k2: kl, src: SRC_ALT });
            }
        }
    }
    out
}
