//! Repeat-copy scoring: for each unit, its best near-identical match among a base row's exons.

use ahash::AHashMap;
use pyo3::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

const K: usize = 7;
const MIN_IDENTITY: f64 = 0.90;
const MIN_COVER: f64 = 0.90;
const SAME_MAX: f64 = 1.2;
const MATCH: i32 = 1;
const MISMATCH: i32 = -1;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = -1;
const NEG: i32 = i32::MIN / 4;

fn code(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// 2-bit k-mer at each position; None where it holds a non-ACGT base.
fn kmers(c: &[u8]) -> Vec<Option<u32>> {
    if c.len() < K {
        return Vec::new();
    }
    (0..=c.len() - K)
        .map(|i| {
            let w = &c[i..i + K];
            if w.iter().any(|&x| x > 3) { None } else { Some(w.iter().fold(0u32, |v, &x| (v << 2) | x as u32)) }
        })
        .collect()
}

#[derive(Clone, Copy)]
struct Cell {
    s: i32,
    ident: u32,
    cols: u32,
    start: u32,
}

const DEAD: Cell = Cell { s: NEG, ident: 0, cols: 0, start: 0 };

/// Local alignment of u against e on diagonals d - w..=d + w: (identity, u start, u end); affine gaps.
fn banded(u: &[u8], e: &[u8], d: i64, w: i64) -> (f64, usize, usize) {
    let (n, m) = (u.len() as i64, e.len() as i64);
    let width = (2 * w + 1) as usize;
    // per band offset b, column j = i + d - w + b
    let mut h_prev = vec![DEAD; width + 2];
    let mut f_prev = vec![DEAD; width + 2];
    let mut h_cur = vec![DEAD; width + 2];
    let mut f_cur = vec![DEAD; width + 2];
    let mut best = Cell { s: 0, ident: 0, cols: 1, start: 0 };
    let mut best_end = 0usize;
    for i in 1..=n {
        let mut e_left = DEAD;   // gap in u, moving along e within this row
        let mut h_left = DEAD;
        for b in 0..width {
            let j = i + d - w + b as i64;
            if j < 1 || j > m {
                h_cur[b + 1] = DEAD;
                f_cur[b + 1] = DEAD;
                e_left = DEAD;
                h_left = DEAD;
                continue;
            }
            let (ui, ej) = (u[(i - 1) as usize], e[(j - 1) as usize]);
            let hit = ui == ej && ui < 4;
            // diagonal: (i-1, j-1) is offset b in the previous row
            let p = h_prev[b + 1];
            let diag = if p.s > 0 {
                Cell { s: p.s + if hit { MATCH } else { MISMATCH }, ident: p.ident + hit as u32, cols: p.cols + 1, start: p.start }
            } else {
                Cell { s: if hit { MATCH } else { MISMATCH }, ident: hit as u32, cols: 1, start: (i - 1) as u32 }
            };
            // gap in e: (i-1, j) is offset b+1 in the previous row
            let (hu, fu) = (h_prev[b + 2], f_prev[b + 2]);
            let f = if hu.s + GAP_OPEN >= fu.s + GAP_EXTEND {
                Cell { s: hu.s + GAP_OPEN, cols: hu.cols + 1, ..hu }
            } else {
                Cell { s: fu.s + GAP_EXTEND, cols: fu.cols + 1, ..fu }
            };
            // gap in u: (i, j-1) is offset b-1 in this row
            let el = if h_left.s + GAP_OPEN >= e_left.s + GAP_EXTEND {
                Cell { s: h_left.s + GAP_OPEN, cols: h_left.cols + 1, ..h_left }
            } else {
                Cell { s: e_left.s + GAP_EXTEND, cols: e_left.cols + 1, ..e_left }
            };
            let mut h = diag;
            if f.s > h.s { h = f; }
            if el.s > h.s { h = el; }
            if h.s <= 0 { h = Cell { s: 0, ident: 0, cols: 0, start: i as u32 }; }
            if h.s > best.s {
                best = h;
                best_end = i as usize;
            }
            h_cur[b + 1] = h;
            f_cur[b + 1] = if f.s < NEG / 2 { DEAD } else { f };
            e_left = if el.s < NEG / 2 { DEAD } else { el };
            h_left = h;
        }
        std::mem::swap(&mut h_prev, &mut h_cur);
        std::mem::swap(&mut f_prev, &mut f_cur);
    }
    if best.s <= 0 {
        return (0.0, 0, 0);
    }
    (best.ident as f64 / best.cols.max(1) as f64, best.start as usize, best_end)
}

/// Best match of each unit among the exons: (exon index or -1, identity, cover, length ratio).
fn score_group(exons: &[Vec<u8>], units: &[Vec<u8>]) -> Vec<(i64, f64, f64, f64)> {
    let mut index: AHashMap<u32, Vec<(u32, u32)>> = AHashMap::new();
    for (x, e) in exons.iter().enumerate() {
        for (y, k) in kmers(e).into_iter().enumerate() {
            if let Some(k) = k { index.entry(k).or_default().push((x as u32, y as u32)); }
        }
    }
    units
        .iter()
        .map(|u| {
            let mut best = (-1i64, 0.0f64, 0.0f64, 0.0f64);
            if u.len() < K {
                return best;
            }
            let mut votes: Vec<AHashMap<i64, u32>> = vec![AHashMap::new(); exons.len()];
            let mut total = vec![0u32; exons.len()];
            for (x, k) in kmers(u).into_iter().enumerate() {
                if let Some(k) = k {
                    if let Some(hits) = index.get(&k) {
                        for &(ex, y) in hits {
                            *votes[ex as usize].entry(y as i64 - x as i64).or_insert(0) += 1;
                            total[ex as usize] += 1;
                        }
                    }
                }
            }
            // q-gram lemma: k-mers any copy at MIN_IDENTITY over MIN_COVER must share
            let span = (MIN_COVER * u.len() as f64).ceil();
            let errors = ((1.0 - MIN_IDENTITY) * span * 1.2).ceil();
            let need = ((span + 1.0 - (errors + 1.0) * K as f64) as i64).max(1) as u32;
            let w = ((0.15 * u.len() as f64).ceil() as i64).max(16);
            for (ex, e) in exons.iter().enumerate() {
                if total[ex] < need || u.len() as f64 > SAME_MAX * e.len() as f64 {
                    continue;
                }
                let (&d, _) = votes[ex].iter().max_by_key(|(&d, &c)| (c, -d)).unwrap();
                // a copy's shared k-mers sit within the band around its diagonal
                if votes[ex].iter().filter(|(&x, _)| (x - d).abs() <= w).map(|(_, &c)| c).sum::<u32>() < need {
                    continue;
                }
                let (ident, s, t) = banded(u, e, d, w);
                let cover = (t - s) as f64 / u.len() as f64;
                if ident * cover.min(1.0) > best.1 * best.2.min(1.0) {
                    best = (ex as i64, ident, cover, u.len() as f64 / e.len() as f64);
                }
            }
            best
        })
        .collect()
}

/// Per group (base exons, units): each unit's best match as (exon or -1, identity, cover, length ratio).
#[pyfunction]
#[pyo3(signature = (groups, threads = 1))]
pub fn repeat_copy_matches(py: Python<'_>, groups: Vec<(Vec<String>, Vec<String>)>, threads: usize) -> Vec<Vec<(i64, f64, f64, f64)>> {
    let coded: Vec<(Vec<Vec<u8>>, Vec<Vec<u8>>)> = groups
        .iter()
        .map(|(ex, us)| (ex.iter().map(|s| s.bytes().map(code).collect()).collect(), us.iter().map(|s| s.bytes().map(code).collect()).collect()))
        .collect();
    py.detach(|| {
        let n = coded.len();
        let next = AtomicUsize::new(0);
        let mut out: Vec<(usize, Vec<(i64, f64, f64, f64)>)> = Vec::with_capacity(n);
        std::thread::scope(|sc| {
            let hs: Vec<_> = (0..threads.max(1).min(n.max(1)))
                .map(|_| sc.spawn(|| {
                    let mut got = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= n { break; }
                        got.push((i, score_group(&coded[i].0, &coded[i].1)));
                    }
                    got
                }))
                .collect();
            for h in hs { out.extend(h.join().expect("repeat_copy worker panicked")); }
        });
        out.sort_by_key(|x| x.0);
        out.into_iter().map(|x| x.1).collect()
    })
}

#[cfg(test)]
#[path = "tests/repeat_copy.rs"]
mod tests;
