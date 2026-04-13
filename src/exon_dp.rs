use flate2::read::GzDecoder;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use rocksdb::{Options, DB};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::Read;
use std::path::Path;

const INTRON_OPEN_PENALTY: f64 = -12.0;
const MIN_INTRON_NT: usize = 30;
const MAX_INTRON_NT: usize = 200_000;
const MIN_EXON_AA_DEFAULT: usize = 10;
const SCORE_FLOOR_FRAC_DEFAULT: f64 = 0.4;
const SEED_THRESHOLD: usize = 50_000;
const _MINIMUM_GAP_BP: usize = 15;
const _MAX_INTERNAL_GAP_BP: usize = 15_000;
const REF_COVERAGE_THRESH: f64 = 0.7;
const REF_GAP_THRESH: f64 = 0.7;
const MIN_CONSEC_CHAR: usize = 5;

const NEG_INF: f64 = f64::NEG_INFINITY;
const TR_START: u8 = 0;
const TR_EXTEND: u8 = 1;
const TR_INTRON: u8 = 2;
const TR_SPLICE_START: u8 = 3;

const BLOSUM_AAS: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";
const NUM_AA: usize = 20;

#[rustfmt::skip]
const BLOSUM62: [[i8; 20]; 20] = [
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],
];

const BLOSUM_DEFAULT: f64 = -4.0;

fn build_aa_idx_table() -> [u8; 256] {
    let mut t = [255u8; 256];
    for (i, &aa) in BLOSUM_AAS.iter().enumerate() {
        t[aa as usize] = i as u8;
    }
    t
}

fn codon_to_aa(c1: u8, c2: u8, c3: u8) -> u8 {
    match (c1, c2, c3) {
        (b'T', b'T', b'T') | (b'T', b'T', b'C') => b'F',
        (b'T', b'T', b'A') | (b'T', b'T', b'G') => b'L',
        (b'C', b'T', b'T') | (b'C', b'T', b'C') | (b'C', b'T', b'A') | (b'C', b'T', b'G') => b'L',
        (b'A', b'T', b'T') | (b'A', b'T', b'C') | (b'A', b'T', b'A') => b'I',
        (b'A', b'T', b'G') => b'M',
        (b'G', b'T', b'T') | (b'G', b'T', b'C') | (b'G', b'T', b'A') | (b'G', b'T', b'G') => b'V',
        (b'T', b'C', b'T') | (b'T', b'C', b'C') | (b'T', b'C', b'A') | (b'T', b'C', b'G') => b'S',
        (b'C', b'C', b'T') | (b'C', b'C', b'C') | (b'C', b'C', b'A') | (b'C', b'C', b'G') => b'P',
        (b'A', b'C', b'T') | (b'A', b'C', b'C') | (b'A', b'C', b'A') | (b'A', b'C', b'G') => b'T',
        (b'G', b'C', b'T') | (b'G', b'C', b'C') | (b'G', b'C', b'A') | (b'G', b'C', b'G') => b'A',
        (b'T', b'A', b'T') | (b'T', b'A', b'C') => b'Y',
        (b'T', b'A', b'A') | (b'T', b'A', b'G') => b'*',
        (b'C', b'A', b'T') | (b'C', b'A', b'C') => b'H',
        (b'C', b'A', b'A') | (b'C', b'A', b'G') => b'Q',
        (b'A', b'A', b'T') | (b'A', b'A', b'C') => b'N',
        (b'A', b'A', b'A') | (b'A', b'A', b'G') => b'K',
        (b'G', b'A', b'T') | (b'G', b'A', b'C') => b'D',
        (b'G', b'A', b'A') | (b'G', b'A', b'G') => b'E',
        (b'T', b'G', b'T') | (b'T', b'G', b'C') => b'C',
        (b'T', b'G', b'A') => b'*',
        (b'T', b'G', b'G') => b'W',
        (b'C', b'G', b'T') | (b'C', b'G', b'C') | (b'C', b'G', b'A') | (b'C', b'G', b'G') => b'R',
        (b'A', b'G', b'T') | (b'A', b'G', b'C') => b'S',
        (b'A', b'G', b'A') | (b'A', b'G', b'G') => b'R',
        (b'G', b'G', b'T') | (b'G', b'G', b'C') | (b'G', b'G', b'A') | (b'G', b'G', b'G') => b'G',
        _ => 0,
    }
}

fn build_codon_table(aa_idx: &[u8; 256]) -> HashMap<[u8; 3], u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut m = HashMap::new();
    for &b1 in &bases {
        for &b2 in &bases {
            for &b3 in &bases {
                let aa = codon_to_aa(b1, b2, b3);
                if aa != 0 && aa != b'*' {
                    let idx = aa_idx[aa as usize];
                    if (idx as usize) < NUM_AA {
                        m.insert([b1, b2, b3], idx);
                    }
                }
            }
        }
    }
    m
}

#[inline(always)]
fn translate_codon_byte(g: &[u8], pos: usize) -> u8 {
    if pos + 2 >= g.len() {
        return 0;
    }
    codon_to_aa(g[pos], g[pos + 1], g[pos + 2])
}

fn translate_seq(nt: &[u8]) -> String {
    let mut out = String::with_capacity(nt.len() / 3);
    let mut i = 0;
    while i + 2 < nt.len() {
        let aa = codon_to_aa(nt[i], nt[i + 1], nt[i + 2]);
        out.push(if aa == 0 { 'X' } else { aa as char });
        i += 3;
    }
    out
}

fn bio_revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            other => other,
        })
        .collect()
}

#[inline(always)]
fn score_donor(g: &[u8], p: usize) -> f64 {
    if p + 1 >= g.len() {
        return NEG_INF;
    }
    match (g[p], g[p + 1]) {
        (b'G', b'T') => 2.0,
        (b'G', b'C') => 1.0,
        (b'A', b'T') => 0.5,
        _ => NEG_INF,
    }
}

#[inline(always)]
fn score_acceptor(g: &[u8], p: usize) -> f64 {
    if p + 1 >= g.len() {
        return NEG_INF;
    }
    match (g[p], g[p + 1]) {
        (b'A', b'G') => 2.0,
        (b'A', b'C') => 0.5,
        _ => NEG_INF,
    }
}

#[inline(always)]
fn donor_val(g: &[u8], p: usize) -> f64 {
    if p + 1 >= g.len() {
        return 0.0;
    }
    match (g[p], g[p + 1]) {
        (b'G', b'T') => 2.0,
        (b'G', b'C') => 1.0,
        (b'A', b'T') => 0.5,
        _ => 0.0,
    }
}

#[inline(always)]
fn acceptor_val(g: &[u8], p: usize) -> f64 {
    if p + 1 >= g.len() {
        return 0.0;
    }
    match (g[p], g[p + 1]) {
        (b'A', b'G') => 2.0,
        (b'A', b'C') => 0.5,
        _ => 0.0,
    }
}

fn build_blosum_pssm(
    consensus: &HashMap<usize, Vec<u8>>,
    ref_cols: &[usize],
    aa_idx: &[u8; 256],
) -> (Vec<[f64; NUM_AA]>, Vec<f64>) {
    let nc = ref_cols.len();
    let mut pssm = vec![[0.0f64; NUM_AA]; nc];
    let mut maxes = vec![0.0f64; nc];
    for (ci, &col) in ref_cols.iter().enumerate() {
        let ra = match consensus.get(&col) {
            Some(v) if !v.is_empty() => v,
            _ => {
                pssm[ci] = [NEG_INF; NUM_AA];
                continue;
            }
        };
        let n = ra.len() as f64;
        let mut acc = [0.0f64; NUM_AA];
        for &r in ra {
            let ri = aa_idx[r as usize];
            if (ri as usize) < NUM_AA {
                for qi in 0..NUM_AA {
                    acc[qi] += BLOSUM62[qi][ri as usize] as f64;
                }
            } else {
                for qi in 0..NUM_AA {
                    acc[qi] += BLOSUM_DEFAULT;
                }
            }
        }
        let mut mx = NEG_INF;
        for qi in 0..NUM_AA {
            let v = acc[qi] / n;
            pssm[ci][qi] = v;
            if v > mx {
                mx = v;
            }
        }
        maxes[ci] = mx;
    }
    (pssm, maxes)
}

#[derive(Clone, Debug)]
struct ExonHit {
    pssm_start: usize,
    pssm_end: usize,
    genome_start: usize,
    genome_end: usize,
    frame: usize,
    score: f64,
    aa_seq: String,
}

type DpRow = HashMap<usize, (f64, usize, u8)>;

fn dp_align(
    genome: &[u8],
    pssm: &[[f64; NUM_AA]],
    _col_maxes: &[f64],
    ct: &HashMap<[u8; 3], u8>,
) -> (f64, isize, isize, Vec<DpRow>) {
    let g = genome.len();
    let n = pssm.len();
    if g < 6 || n < 1 {
        return (NEG_INF, -1, -1, Vec::new());
    }

    let mut cidx: Vec<i8> = vec![-1; g];
    for j in 0..g.saturating_sub(2) {
        if let Some(&idx) = ct.get(&[genome[j], genome[j + 1], genome[j + 2]]) {
            cidx[j] = idx as i8;
        }
    }

    let mut darr = vec![0.0f64; g];
    let mut aarr = vec![0.0f64; g];
    for j in 0..g.saturating_sub(1) {
        darr[j] = donor_val(genome, j);
        aarr[j] = acceptor_val(genome, j);
    }

    let mut rows: Vec<DpRow> = (0..n).map(|_| HashMap::new()).collect();
    let mut bs = NEG_INF;
    let mut bi: isize = -1;
    let mut bj: isize = -1;

    {
        let pr = &pssm[0];
        let row = &mut rows[0];
        for j in 0..g.saturating_sub(2) {
            let ai = cidx[j];
            if ai < 0 {
                continue;
            }
            let ms = pr[ai as usize];
            if ms == NEG_INF {
                continue;
            }
            let mut sb = NEG_INF;
            let mut tr = TR_START;
            if ms > 0.0 {
                sb = ms;
                tr = TR_START;
            }
            if j >= 2 {
                let ac = aarr[j - 2];
                if ac > 0.0 {
                    let ss = ac + ms;
                    if ss > sb {
                        sb = ss;
                        tr = TR_SPLICE_START;
                    }
                }
            }
            if sb > NEG_INF {
                row.insert(j, (sb, usize::MAX, tr));
                if sb > bs {
                    bs = sb;
                    bi = 0;
                    bj = j as isize;
                }
            }
        }
    }

    for i in 1..n {
        let pr = &pssm[i];
        let prev_empty = rows[i - 1].is_empty();

        if prev_empty {
            let row = &mut rows[i];
            for j in 0..g.saturating_sub(2) {
                let ai = cidx[j];
                if ai < 0 {
                    continue;
                }
                let ms = pr[ai as usize];
                if ms == NEG_INF {
                    continue;
                }
                let mut sb = NEG_INF;
                let mut tr = TR_START;
                if ms > 0.0 {
                    sb = ms;
                    tr = TR_START;
                }
                if j >= 2 {
                    let ac = aarr[j - 2];
                    if ac > 0.0 {
                        let ss = ac + ms;
                        if ss > sb {
                            sb = ss;
                            tr = TR_SPLICE_START;
                        }
                    }
                }
                if sb > NEG_INF {
                    row.insert(j, (sb, usize::MAX, tr));
                    if sb > bs {
                        bs = sb;
                        bi = i as isize;
                        bj = j as isize;
                    }
                }
            }
            continue;
        }

        let mut eps: Vec<(usize, f64)> = Vec::new();
        for (&jp, &(psc, _, _)) in &rows[i - 1] {
            let dp = jp + 3;
            if dp < g {
                let ds = darr[dp];
                if ds > 0.0 {
                    eps.push((jp, psc + INTRON_OPEN_PENALTY + ds));
                }
            }
        }
        eps.sort_unstable_by_key(|&(jp, _)| jp);

        let mut epi = 0usize;
        let epl = eps.len();
        let mut rbs = NEG_INF;
        let mut rbj = 0usize;

        let prev: HashMap<usize, f64> = rows[i - 1].iter().map(|(&k, &(s, _, _))| (k, s)).collect();
        let row = &mut rows[i];

        for j in 0..g.saturating_sub(2) {
            let ai = cidx[j];
            if ai < 0 {
                continue;
            }
            let ms = pr[ai as usize];
            if ms == NEG_INF {
                continue;
            }

            let mut sb = NEG_INF;
            let mut src = usize::MAX;
            let mut trb = TR_START;

            if ms > 0.0 {
                sb = ms;
                src = usize::MAX;
                trb = TR_START;
            }

            if j >= 2 {
                let ac = aarr[j - 2];
                if ac > 0.0 {
                    let ss = ac + ms;
                    if ss > sb {
                        sb = ss;
                        src = usize::MAX;
                        trb = TR_SPLICE_START;
                    }
                }
            }

            if j >= 3 {
                if let Some(&psc) = prev.get(&(j - 3)) {
                    let sc = psc + ms;
                    if sc > sb {
                        sb = sc;
                        src = j - 3;
                        trb = TR_EXTEND;
                    }
                }
            }

            if j >= 2 {
                let ac = aarr[j - 2];
                if ac > 0.0 && j >= MIN_INTRON_NT + 3 {
                    let mep = j - MIN_INTRON_NT - 3;
                    while epi < epl && eps[epi].0 <= mep {
                        let (ej, es) = eps[epi];
                        if es > rbs {
                            rbs = es;
                            rbj = ej;
                        }
                        epi += 1;
                    }
                    if rbs > NEG_INF {
                        let il = j - rbj - 3;
                        if il <= MAX_INTRON_NT {
                            let sc = rbs + ac + ms;
                            if sc > sb {
                                sb = sc;
                                src = rbj;
                                trb = TR_INTRON;
                            }
                        }
                    }
                }
            }

            if sb > NEG_INF {
                row.insert(j, (sb, src, trb));
                if sb > bs {
                    bs = sb;
                    bi = i as isize;
                    bj = j as isize;
                }
            }
        }
    }

    (bs, bi, bj, rows)
}

fn traceback(bi: isize, bj: isize, rows: &[DpRow]) -> Vec<(usize, usize, u8)> {
    if bi < 0 || bj < 0 {
        return Vec::new();
    }
    let mut path = Vec::new();
    let (mut i, mut j) = (bi as usize, bj as usize);
    loop {
        let &(_, src, tr) = match rows[i].get(&j) {
            Some(entry) => entry,
            None => break,
        };
        path.push((i, j, tr));
        if tr == TR_START || tr == TR_SPLICE_START {
            break;
        }
        if (tr == TR_EXTEND || tr == TR_INTRON) && i > 0 {
            i -= 1;
            j = src;
        } else {
            break;
        }
    }
    path.reverse();
    path
}

fn extract_exons(
    path: &[(usize, usize, u8)],
    genome: &[u8],
    pssm: &[[f64; NUM_AA]],
    aa_idx: &[u8; 256],
    min_aa: usize,
) -> Vec<ExonHit> {
    if path.is_empty() {
        return Vec::new();
    }

    let mut exons = Vec::new();
    let mut cols: Vec<usize> = Vec::new();
    let mut gps: Vec<usize> = Vec::new();
    let mut sc = 0.0;
    let mut aas: Vec<u8> = Vec::new();

    let flush = |c: &[usize], g: &[usize], s: f64, a: &[u8], out: &mut Vec<ExonHit>| {
        if c.len() >= min_aa {
            out.push(ExonHit {
                pssm_start: c[0],
                pssm_end: c[c.len() - 1] + 1,
                genome_start: g[0],
                genome_end: g[g.len() - 1] + 3,
                frame: g[0] % 3,
                score: s,
                aa_seq: String::from_utf8_lossy(a).into_owned(),
            });
        }
    };

    for &(ci, gp, tr) in path {
        let aa = translate_codon_byte(genome, gp);
        let ai = aa_idx[aa as usize];
        let ms = if (ai as usize) < NUM_AA {
            pssm[ci][ai as usize]
        } else {
            0.0
        };

        if tr == TR_INTRON || tr == TR_START || tr == TR_SPLICE_START {
            if !cols.is_empty() {
                flush(&cols, &gps, sc, &aas, &mut exons);
            }
            cols = vec![ci];
            gps = vec![gp];
            sc = ms;
            aas = vec![aa];
        } else if tr == TR_EXTEND {
            cols.push(ci);
            gps.push(gp);
            sc += ms;
            aas.push(aa);
        } else {
            if !cols.is_empty() {
                flush(&cols, &gps, sc, &aas, &mut exons);
            }
            cols = vec![ci];
            gps = vec![gp];
            sc = ms;
            aas = vec![aa];
        }
    }
    if !cols.is_empty() {
        flush(&cols, &gps, sc, &aas, &mut exons);
    }
    exons
}

fn extend_to_splice(
    exons: &[ExonHit],
    genome: &[u8],
    pssm: &[[f64; NUM_AA]],
    _col_maxes: &[f64],
    aa_idx: &[u8; 256],
    max_ext: usize,
) -> Vec<ExonHit> {
    if exons.is_empty() {
        return Vec::new();
    }
    let g = genome.len();
    let n = pssm.len();
    let mut out = Vec::with_capacity(exons.len());

    for exon in exons {
        let (mut nps, mut npe, mut ngs, mut nge) = (
            exon.pssm_start,
            exon.pssm_end,
            exon.genome_start,
            exon.genome_end,
        );
        let mut xsc = 0.0;
        let mut front: Vec<u8> = Vec::new();
        let mut back: Vec<u8> = Vec::new();

        let has_acc = ngs >= 2 && genome[ngs - 2] == b'A' && genome[ngs - 1] == b'G';
        if !has_acc && nps > 0 {
            let (mut bg, mut bp, mut bsc, mut baas, mut bacc) = (ngs, nps, 0.0, Vec::new(), -1.0);
            let (mut cg, mut cp, mut rsc) = (ngs, nps, 0.0);
            let mut raas: Vec<u8> = Vec::new();
            for _ in 1..=max_ext {
                if cg < 3 || cp == 0 {
                    break;
                }
                let (ng, np) = (cg - 3, cp - 1);
                let aa = translate_codon_byte(genome, ng);
                if aa == 0 || aa == b'X' || aa == b'*' {
                    break;
                }
                let ai = aa_idx[aa as usize];
                let ms = if (ai as usize) < NUM_AA {
                    pssm[np][ai as usize]
                } else {
                    0.0
                };
                rsc += ms;
                raas.insert(0, aa);
                cg = ng;
                cp = np;
                if cg >= 2 {
                    let asc = score_acceptor(genome, cg - 2);
                    if asc > NEG_INF && asc > bacc {
                        bg = cg;
                        bp = cp;
                        bsc = rsc;
                        baas = raas.clone();
                        bacc = asc;
                        if genome[cg - 2] == b'A' && genome[cg - 1] == b'G' {
                            break;
                        }
                    }
                }
            }
            if bg < ngs {
                ngs = bg;
                nps = bp;
                xsc += bsc;
                front = baas;
            }
        }

        let has_don = nge + 1 < g && genome[nge] == b'G' && genome[nge + 1] == b'T';
        if !has_don && npe < n {
            let (mut bg, mut bp, mut bsc, mut baas, mut bdon) = (nge, npe, 0.0, Vec::new(), -1.0);
            let (mut cg, mut cp, mut rsc) = (nge, npe, 0.0);
            let mut raas: Vec<u8> = Vec::new();
            for _ in 1..=max_ext {
                let (ng, np) = (cg + 3, cp + 1);
                if ng > g || np > n || cp >= n {
                    break;
                }
                let aa = translate_codon_byte(genome, cg);
                if aa == 0 || aa == b'X' || aa == b'*' {
                    break;
                }
                let ai = aa_idx[aa as usize];
                let ms = if (ai as usize) < NUM_AA {
                    pssm[cp][ai as usize]
                } else {
                    0.0
                };
                rsc += ms;
                raas.push(aa);
                cg = ng;
                cp = np;
                if cg + 1 < g {
                    let dsc = score_donor(genome, cg);
                    if dsc > NEG_INF && dsc > bdon {
                        bg = cg;
                        bp = cp;
                        bsc = rsc;
                        baas = raas.clone();
                        bdon = dsc;
                        if genome[cg] == b'G' && genome[cg + 1] == b'T' {
                            break;
                        }
                    }
                }
            }
            if bg > nge {
                nge = bg;
                npe = bp;
                xsc += bsc;
                back = baas;
            }
        }

        if ngs != exon.genome_start || nge != exon.genome_end {
            let mut aa = String::from_utf8_lossy(&front).into_owned();
            aa.push_str(&exon.aa_seq);
            aa.push_str(&String::from_utf8_lossy(&back));
            out.push(ExonHit {
                pssm_start: nps,
                pssm_end: npe,
                genome_start: ngs,
                genome_end: nge,
                frame: ngs % 3,
                score: exon.score + xsc,
                aa_seq: aa,
            });
        } else {
            out.push(exon.clone());
        }
    }
    out
}

fn seed_regions(
    genome: &[u8],
    pssm: &[[f64; NUM_AA]],
    col_maxes: &[f64],
    aa_idx: &[u8; 256],
) -> Vec<(usize, usize)> {
    let g = genome.len();
    let n = pssm.len();
    let wa = 15.min(n);
    let mut ps = vec![0, n.saturating_sub(wa)];
    if n >= wa {
        ps.push(n / 2 - wa / 2);
    }
    ps.sort_unstable();
    ps.dedup();

    let mut caa = vec![0u8; g];
    for j in 0..g.saturating_sub(2) {
        caa[j] = codon_to_aa(
            genome[j].to_ascii_uppercase(),
            genome[j + 1].to_ascii_uppercase(),
            genome[j + 2].to_ascii_uppercase(),
        );
    }

    let mut seeds: Vec<(f64, usize)> = Vec::new();
    for &pstart in &ps {
        let pend = (pstart + wa).min(n);
        let plen = pend - pstart;
        let pp = &pssm[pstart..pend];
        let pm: f64 = col_maxes[pstart..pend].iter().sum();
        let thr = pm * SCORE_FLOOR_FRAC_DEFAULT;
        for frame in 0..3 {
            let mut j = frame;
            while j + plen * 3 <= g {
                let mut sc = 0.0;
                for ci in 0..plen {
                    let aa = caa[j + ci * 3];
                    if aa != 0 && aa != b'*' {
                        let ai = aa_idx[aa as usize];
                        if (ai as usize) < NUM_AA {
                            sc += pp[ci][ai as usize];
                        }
                    }
                }
                if sc >= thr {
                    seeds.push((sc, j + (plen * 3) / 2));
                }
                j += 3;
            }
        }
    }
    if seeds.is_empty() {
        return Vec::new();
    }
    seeds.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    let mut regs: Vec<(usize, usize)> = Vec::new();
    for &(_, center) in &seeds {
        let rs = center.saturating_sub(5000);
        let re = (center + 5000).min(g);
        if regs.iter().any(|&(s, e)| rs < e && re > s) {
            continue;
        }
        regs.push((rs, re));
        if regs.len() >= 20 {
            break;
        }
    }
    regs.sort();
    regs
}

type ExonTuple = (usize, usize, usize, usize, usize, f64, String);

fn refine_boundaries(exons: &[ExonTuple], genome: &[u8]) -> Vec<ExonTuple> {
    if exons.len() < 2 {
        return exons.to_vec();
    }
    let mut r = exons.to_vec();
    for idx in 0..r.len() - 1 {
        let a_end = r[idx].3;
        let b_start = r[idx + 1].2;
        if b_start <= a_end {
            continue;
        }
        let cd = score_donor(genome, a_end);
        let ca = if b_start >= 2 {
            score_acceptor(genome, b_start - 2)
        } else {
            NEG_INF
        };
        let cs = (if cd > NEG_INF { cd } else { 0.0 }) + (if ca > NEG_INF { ca } else { 0.0 });
        let mut bsh: i32 = 0;
        let mut bsp = cs;
        for sh in [-2i32, -1, 1, 2] {
            let nae = a_end as i64 + sh as i64;
            let nbs = b_start as i64 + sh as i64;
            if nae < r[idx].2 as i64 + 3 || nbs + 3 > r[idx + 1].3 as i64 {
                continue;
            }
            let (nae, nbs) = (nae as usize, nbs as usize);
            if nbs < 2 {
                continue;
            }
            let il = (nbs - 2) as isize - nae as isize;
            if il < MIN_INTRON_NT as isize {
                continue;
            }
            let d = score_donor(genome, nae);
            let a = score_acceptor(genome, nbs - 2);
            let sp = (if d > NEG_INF { d } else { 0.0 }) + (if a > NEG_INF { a } else { 0.0 });
            if sp > bsp {
                bsp = sp;
                bsh = sh;
            }
        }
        if bsh != 0 {
            let (psa, pea, gsa, gea, _, sca, _) = r[idx].clone();
            let (psb, peb, gsb, geb, _, scb, _) = r[idx + 1].clone();
            let ngea = (gea as i64 + bsh as i64) as usize;
            let ngsb = (gsb as i64 + bsh as i64) as usize;
            let sc: Vec<u8> = if bsh.abs() == 1 {
                if bsh > 0 {
                    let mut c = genome[gea..gea + 1].to_vec();
                    c.extend_from_slice(&genome[ngsb..ngsb + 2]);
                    c
                } else {
                    let mut c = genome[ngea - 1..ngea].to_vec();
                    c.extend_from_slice(&genome[gsb..gsb + 2]);
                    c
                }
            } else if bsh.abs() == 2 {
                if bsh > 0 {
                    let mut c = genome[gea..gea + 2].to_vec();
                    c.extend_from_slice(&genome[ngsb..ngsb + 1]);
                    c
                } else {
                    let mut c = genome[ngea - 2..ngea].to_vec();
                    c.extend_from_slice(&genome[gsb..gsb + 1]);
                    c
                }
            } else {
                Vec::new()
            };
            if sc.len() == 3 {
                let saa = codon_to_aa(sc[0], sc[1], sc[2]);
                if saa != b'*' && saa != 0 {
                    // Re-translate each exon over its new (possibly non
                    // codon-aligned) span, then snap the NT range back to
                    // an exact codon-aligned subrange of length aa.len()*3.
                    //
                    // Why: bsh ∈ {-2,-1,1,2} shifts the splice boundary by
                    // 1-2 nt, which is biologically valid (splice sites are
                    // not constrained to codon boundaries) but leaves a
                    // partial codon at one end of each exon. translate_seq
                    // ignores those trailing nt, so the AA length stays the
                    // same while the NT span gains 1-2 untranslated nt.
                    // The "spanning codon" saa above accounts for them
                    // logically but cannot be appended to either exon's NT
                    // (the 3 nt straddle the intron and are not contiguous).
                    // Trim them so downstream consumers (pn2codon, etc.)
                    // see NT length == aa.len() * 3.
                    let aa_left = translate_seq(&genome[gsa..ngea]);
                    let nt_end_left = gsa + aa_left.len() * 3;
                    r[idx] = (
                        psa,
                        pea,
                        gsa,
                        nt_end_left,
                        gsa % 3,
                        sca,
                        aa_left,
                    );
                    let aa_right = translate_seq(&genome[ngsb..geb]);
                    let nt_end_right = ngsb + aa_right.len() * 3;
                    r[idx + 1] = (
                        psb,
                        peb,
                        ngsb,
                        nt_end_right,
                        ngsb % 3,
                        scb,
                        aa_right,
                    );
                }
            }
        }
    }
    r
}

fn solve_gap_dp(
    splice: &[u8],
    consensus: &HashMap<usize, Vec<u8>>,
    ref_cols: &[usize],
    _ref_gaps: &HashSet<usize>,
    aa_idx: &[u8; 256],
    ct: &HashMap<[u8; 3], u8>,
    log: &mut Vec<String>,
    score_thr: f64,
    min_aa: usize,
) -> Vec<ExonTuple> {
    let nc = ref_cols.len();
    if nc < min_aa {
        log.push(format!("DP: too few PSSM columns ({}<{})", nc, min_aa));
        return Vec::new();
    }

    let (pssm, cmx) = build_blosum_pssm(consensus, ref_cols, aa_idx);
    let total_max: f64 = cmx.iter().sum();
    if total_max <= 0.0 {
        log.push("DP: PSSM has no positive scoring columns".into());
        return Vec::new();
    }

    let genome: Vec<u8> = splice.iter().map(|b| b.to_ascii_uppercase()).collect();
    log.push(format!(
        "DP: {} PSSM cols, genome={} nt, max_possible={:.1}",
        nc,
        genome.len(),
        total_max
    ));

    let regions = if genome.len() > SEED_THRESHOLD {
        log.push(format!(
            "DP: genome > {} nt, using seeded approach",
            SEED_THRESHOLD
        ));
        let r = seed_regions(&genome, &pssm, &cmx, aa_idx);
        if r.is_empty() {
            log.push("DP: no seed regions found".into());
            return Vec::new();
        }
        log.push(format!("DP: {} seed region(s)", r.len()));
        r
    } else {
        vec![(0, genome.len())]
    };

    let mut all_exons: Vec<ExonHit> = Vec::new();

    // Track DP data from the best-scoring region for runner-up analysis
    let mut best_region_score = NEG_INF;
    let mut best_region_data: Option<(Vec<u8>, Vec<DpRow>, usize)> = None;

    for &(rs, re) in &regions {
        let sub = &genome[rs..re];
        let (score, bi, bj, rows) = dp_align(sub, &pssm, &cmx, ct);
        if score <= 0.0 || bi < 0 {
            continue;
        }

        if score > best_region_score {
            best_region_score = score;
            best_region_data = Some((sub.to_vec(), rows.clone(), rs));
        }

        let path = traceback(bi, bj, &rows);
        let exons = extract_exons(&path, sub, &pssm, aa_idx, min_aa);
        let exons = extend_to_splice(&exons, sub, &pssm, &cmx, aa_idx, 15);
        for e in exons {
            all_exons.push(ExonHit {
                genome_start: e.genome_start + rs,
                genome_end: e.genome_end + rs,
                ..e
            });
        }
    }

    if all_exons.is_empty() {
        log.push("DP: no exons recovered".into());
        return Vec::new();
    }

    // Track status for each exon candidate
    let mut exon_statuses: Vec<(ExonHit, String)> = Vec::new();

    // Quality filter
    let mut quality_passed: Vec<ExonHit> = Vec::new();
    for e in all_exons {
        let em: f64 = cmx[e.pssm_start..e.pssm_end].iter().sum();
        if em > 0.0 && e.score / em >= score_thr {
            quality_passed.push(e);
        } else {
            let frac = if em > 0.0 { e.score / em } else { 0.0 };
            exon_statuses.push((e, format!("FAIL:quality({:.1}%<{:.0}%)", frac * 100.0, score_thr * 100.0)));
        }
    }

    // Dedup overlapping PSSM regions
    let mut all_exons = quality_passed;
    if all_exons.len() > 1 {
        all_exons.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
        let mut used: HashSet<usize> = HashSet::new();
        let mut dd: Vec<ExonHit> = Vec::new();
        for e in &all_exons {
            let ec: HashSet<usize> = (e.pssm_start..e.pssm_end).collect();
            if ec.intersection(&used).count() > 8 {
                exon_statuses.push((e.clone(), "FAIL:dedup_overlap".into()));
                continue;
            }
            used.extend(e.pssm_start..e.pssm_end);
            dd.push(e.clone());
        }
        dd.sort_by_key(|e| e.pssm_start);
        all_exons = dd;
    }

    let mut cov: HashSet<usize> = HashSet::new();
    let (mut ts, mut tem) = (0.0, 0.0);
    for e in &all_exons {
        for c in e.pssm_start..e.pssm_end {
            cov.insert(c);
        }
        ts += e.score;
        tem += cmx[e.pssm_start..e.pssm_end].iter().sum::<f64>();
    }
    let sf = if tem > 0.0 { ts / tem } else { 0.0 };

    if all_exons.is_empty() {
        log.push("DP: all main exons filtered by quality".into());
    } else {
        log.push(format!(
            "DP result: {} exon(s), {}/{} cols ({:.1}%), score {:.1}/{:.1} ({:.1}%)",
            all_exons.len(),
            cov.len(),
            nc,
            cov.len() as f64 / nc as f64 * 100.0,
            ts,
            tem,
            sf * 100.0
        ));
    }

    // Mark surviving exons with their status
    for e in &all_exons {
        if sf >= score_thr {
            exon_statuses.push((e.clone(), "PASS".into()));
        } else {
            exon_statuses.push((e.clone(), format!("FAIL:chain({:.1}%)", sf * 100.0)));
        }
    }

    // Log all exon candidates with their status, aligned to gap region
    let gap_len = splice.len();
    exon_statuses.sort_by_key(|(e, _)| e.pssm_start);
    for (i, (e, status)) in exon_statuses.iter().enumerate() {
        let em: f64 = cmx[e.pssm_start..e.pssm_end].iter().sum();
        let frac = if em > 0.0 { e.score / em } else { 0.0 };
        log.push(format!(
            ">exon{} {} pssm={}-{} ({} AA) frame={} score={:.1}/{:.1} ({:.1}%) genome_nt={}-{} aa={}",
            i + 1, status, e.pssm_start, e.pssm_end,
            e.pssm_end - e.pssm_start, e.frame, e.score, em, frac * 100.0,
            e.genome_start, e.genome_end, e.aa_seq
        ));
        let nt_end = e.genome_end.min(gap_len);
        let nt_slice = String::from_utf8_lossy(&splice[e.genome_start..nt_end]);
        log.push(format!("nt={}", nt_slice));
    }

    // Find and log runner-up candidates, then pick the winner by best score fraction
    // Candidates: index 0 = main hit, 1..=N = runner-ups
    struct Candidate {
        score_frac: f64,
        exons: Vec<ExonHit>, // global coords
    }

    // If main has no exons (all filtered), use NEG_INF so runner-ups can win
    let main_frac = if all_exons.is_empty() { NEG_INF } else { sf };
    let mut candidates: Vec<Candidate> = vec![Candidate {
        score_frac: main_frac,
        exons: all_exons.clone(),
    }];

    if let Some((ref sub_genome, ref rows, r_start)) = best_region_data {
        let best_relative: Vec<ExonHit> = all_exons
            .iter()
            .map(|e| ExonHit {
                genome_start: e.genome_start - r_start,
                genome_end: e.genome_end - r_start,
                ..e.clone()
            })
            .collect();
        let runners = find_runner_ups(
            rows, sub_genome, &pssm, &cmx, &best_relative, aa_idx, min_aa, score_thr, 2,
        );
        for (ri, (r_score, r_frac, r_exons)) in runners.iter().enumerate() {
            log.push(format!(
                "Runner-up #{}: total_score={:.1} ({:.1}% of max), {} exon(s)",
                ri + 1, r_score, r_frac * 100.0, r_exons.len()
            ));
            // Convert runner-up exons to global coords and log them
            let mut global_exons: Vec<ExonHit> = Vec::new();
            for (ei, exon) in r_exons.iter().enumerate() {
                let exon_max: f64 = cmx[exon.pssm_start..exon.pssm_end].iter().sum();
                let efrac = if exon_max > 0.0 { exon.score / exon_max } else { 0.0 };
                let g_start_global = exon.genome_start + r_start;
                let g_end_global = exon.genome_end + r_start;
                log.push(format!(
                    ">exon{} pssm={}-{} ({} AA) frame={} score={:.1}/{:.1} ({:.1}%) genome_nt={}-{} aa={}",
                    ei + 1, exon.pssm_start, exon.pssm_end,
                    exon.pssm_end - exon.pssm_start, exon.frame, exon.score, exon_max, efrac * 100.0,
                    g_start_global, g_end_global, exon.aa_seq
                ));
                let nt_end = g_end_global.min(gap_len);
                let nt_slice = String::from_utf8_lossy(&splice[g_start_global..nt_end]);
                let leading = "-".repeat(g_start_global);
                let trailing = "-".repeat(gap_len.saturating_sub(nt_end));
                log.push(format!("{}{}{}", leading, nt_slice, trailing));

                global_exons.push(ExonHit {
                    genome_start: g_start_global,
                    genome_end: g_end_global,
                    ..exon.clone()
                });
            }
            candidates.push(Candidate {
                score_frac: *r_frac,
                exons: global_exons,
            });
        }
        if runners.is_empty() {
            log.push("Runner-up candidates: none found".into());
        }
    }

    // Pick winner: highest score fraction among all candidates that pass threshold
    let mut best_idx = 0usize;
    let mut best_frac = candidates[0].score_frac;
    for (ci, cand) in candidates.iter().enumerate().skip(1) {
        if cand.score_frac > best_frac {
            best_frac = cand.score_frac;
            best_idx = ci;
        }
    }

    let winner = &candidates[best_idx];

    if winner.score_frac < score_thr {
        log.push(format!(
            "DP: chain rejected -- best candidate {:.1}% < {:.0}%",
            winner.score_frac * 100.0,
            score_thr * 100.0
        ));
        return Vec::new();
    }

    // Log the winner
    let winner_label = if best_idx == 0 {
        "main".to_string()
    } else {
        format!("runner-up #{}", best_idx)
    };
    log.push(format!("Winner (from {}):", winner_label));
    for (ei, exon) in winner.exons.iter().enumerate() {
        let exon_max: f64 = cmx[exon.pssm_start..exon.pssm_end].iter().sum();
        let efrac = if exon_max > 0.0 { exon.score / exon_max } else { 0.0 };
        log.push(format!(
            ">exon{} pssm={}-{} ({} AA) frame={} score={:.1}/{:.1} ({:.1}%) genome_nt={}-{}",
            ei + 1, exon.pssm_start, exon.pssm_end,
            exon.pssm_end - exon.pssm_start, exon.frame, exon.score, exon_max, efrac * 100.0,
            exon.genome_start, exon.genome_end
        ));
        let nt_end = exon.genome_end.min(gap_len);
        let nt_slice = String::from_utf8_lossy(&splice[exon.genome_start..nt_end]);
        log.push(nt_slice.into_owned());
        log.push(exon.aa_seq.clone());
    }

    winner
        .exons
        .iter()
        .map(|e| {
            (
                e.pssm_start,
                e.pssm_end,
                e.genome_start,
                e.genome_end,
                e.frame,
                e.score,
                e.aa_seq.clone(),
            )
        })
        .collect()
}

fn find_runner_ups(
    rows: &[DpRow],
    genome: &[u8],
    pssm: &[[f64; NUM_AA]],
    col_maxes: &[f64],
    best_exons: &[ExonHit],
    aa_idx: &[u8; 256],
    min_aa: usize,
    score_thr: f64,
    n_runners: usize,
) -> Vec<(f64, f64, Vec<ExonHit>)> {
    if rows.is_empty() {
        return Vec::new();
    }

    // Collect genomic intervals already taken by the best chain
    let taken: Vec<(usize, usize)> = best_exons
        .iter()
        .map(|e| (e.genome_start, e.genome_end))
        .collect();

    // Gather all scored endpoints sorted by descending score
    let mut endpoints: Vec<(f64, usize, usize)> = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        for (&j, &(sc, _, _)) in row {
            endpoints.push((sc, i, j));
        }
    }
    endpoints.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

    let mut runners: Vec<(f64, f64, Vec<ExonHit>)> = Vec::new();

    for &(_, ei, ej) in &endpoints {
        if runners.len() >= n_runners {
            break;
        }

        let path = traceback(ei as isize, ej as isize, rows);
        let exons = extract_exons(&path, genome, pssm, aa_idx, min_aa);
        if exons.is_empty() {
            continue;
        }

        // Check that this candidate is NOT a substring of any accepted chain
        let mut all_intervals: Vec<(usize, usize)> = taken.clone();
        for (_, _, rexons) in &runners {
            for e in rexons {
                all_intervals.push((e.genome_start, e.genome_end));
            }
        }

        let is_substring = exons.iter().all(|exon| {
            all_intervals
                .iter()
                .any(|&(ts, te)| exon.genome_start >= ts && exon.genome_end <= te)
        });
        if is_substring {
            continue;
        }

        let total_score: f64 = exons.iter().map(|e| e.score).sum();
        let total_max: f64 = exons
            .iter()
            .map(|e| col_maxes[e.pssm_start..e.pssm_end].iter().sum::<f64>())
            .sum();
        let score_frac = if total_max > 0.0 {
            total_score / total_max
        } else {
            0.0
        };

        runners.push((total_score, score_frac, exons));
    }

    runners
}

fn count_ref_gaps(
    gs: usize,
    ge: usize,
    rc: &HashMap<usize, Vec<u8>>,
) -> (HashSet<usize>, Vec<usize>, usize) {
    let (mut rg, mut ia) = (HashSet::new(), Vec::new());
    let (mut longest, mut cur) = (0usize, 0usize);
    for (ri, col) in (gs..ge).enumerate() {
        let chars = match rc.get(&col) {
            Some(c) => c,
            None => {
                rg.insert(col);
                ia.push(ri);
                longest = longest.max(cur);
                cur = 0;
                continue;
            }
        };
        let gf = if chars.is_empty() {
            1.0
        } else {
            chars.iter().filter(|&&c| c == b'-').count() as f64 / chars.len() as f64
        };
        if gf >= REF_GAP_THRESH {
            rg.insert(col);
            ia.push(ri);
            longest = longest.max(cur);
            cur = 0;
        } else {
            cur += 1;
        }
    }
    longest = longest.max(cur);
    (rg, ia, longest)
}

fn filter_ref_consensus(
    rcount: usize,
    rc: &HashMap<usize, Vec<u8>>,
    rg: &HashSet<usize>,
    gs: usize,
    ge: usize,
    coverage_thresh: f64,
) -> (HashMap<usize, Vec<u8>>, Vec<usize>) {
    // Note: the Python call site passes data_cols_required=(None, None),
    // meaning no edge trimming is performed.  We match that behaviour here
    // by skipping edge trimming entirely.
    let usable = (ge - gs) - rg.len();
    if usable == 0 {
        return (HashMap::new(), Vec::new());
    }
    let mut covs = Vec::with_capacity(rcount);
    for y in 0..rcount {
        let (mut tot, mut ng) = (0usize, 0usize);
        for x in gs..ge {
            if rg.contains(&x) {
                continue;
            }
            tot += 1;
            if let Some(ch) = rc.get(&x) {
                if y < ch.len() && ch[y] != b'-' && ch[y] != b' ' {
                    ng += 1;
                }
            }
        }
        covs.push(if tot > 0 { ng as f64 / tot as f64 } else { 0.0 });
    }
    let best = covs.iter().cloned().fold(0.0f64, f64::max);
    let target = coverage_thresh * best;
    let ri: Vec<usize> = covs
        .iter()
        .enumerate()
        .filter(|(_, c)| **c >= target)
        .map(|(i, _)| i)
        .collect();
    let mut fc: HashMap<usize, Vec<u8>> = HashMap::new();
    for col in gs..ge {
        if rg.contains(&col) {
            continue;
        }
        if let Some(chars) = rc.get(&col) {
            fc.insert(
                col,
                ri.iter()
                    .filter(|&&y| y < chars.len() && chars[y] != b' ' && chars[y] != b'-')
                    .map(|&y| chars[y])
                    .collect(),
            );
        }
    }
    (fc, ri)
}

fn prepare_gap_consensus(
    gs: usize,
    ge: usize,
    rc: &HashMap<usize, Vec<u8>>,
    rcount: usize,
    cc: Option<&HashMap<usize, Vec<u8>>>,
    log: &mut Vec<String>,
) -> Option<(
    HashMap<usize, Vec<u8>>,
    Vec<usize>,
    HashSet<usize>,
    Vec<usize>,
)> {
    let (rg, ia, longest) = count_ref_gaps(gs, ge, rc);
    if longest < MIN_CONSEC_CHAR {
        log.push("Not enough consecutive non-gap columns".into());
        return None;
    }
    let (tc, _ri) = filter_ref_consensus(rcount, rc, &rg, gs, ge, REF_COVERAGE_THRESH);
    if tc.is_empty() {
        log.push("No usable consensus columns after filtering".into());
        return None;
    }
    let rcols: Vec<usize> = (gs..ge).filter(|c| !rg.contains(c)).collect();
    let fc = if let Some(cc) = cc {
        let mut f: HashMap<usize, Vec<u8>> = HashMap::new();
        for &col in &rcols {
            let chars: Vec<u8> = cc.get(&col).map_or(Vec::new(), |v| {
                v.iter()
                    .filter(|&&c| c != b'-' && c != b' ')
                    .copied()
                    .collect()
            });
            f.insert(
                col,
                if chars.is_empty() {
                    tc.get(&col).cloned().unwrap_or_default()
                } else {
                    chars
                },
            );
        }
        f
    } else {
        tc
    };
    Some((fc, rcols, rg, ia))
}

fn find_index_pair(seq: &[u8], gap: u8) -> (usize, usize) {
    let s = seq.iter().position(|&c| c != gap).unwrap_or(0);
    let e = seq.iter().rposition(|&c| c != gap).unwrap_or(0);
    (s, e)
}

fn parse_node_field(h: &str) -> &str {
    h.split('|').nth(3).unwrap_or("")
}

// ---------------------------------------------------------------------------
// Motif scan: kmer-based flank extension
// ---------------------------------------------------------------------------

struct MotifNode {
    header: String,
    frame: i32,
    sequence: String,
    start: usize,
    end: usize,
}

struct MotifHit {
    header: String,
    aa_seq: String,
    nt_seq: String,
    region: String,
    score: usize,
    node_name: String,
    is_leading: bool,
}

/// Build sorted exon interval lists per scaffold from GFF entries.
fn build_scaffold_intervals(
    gff: &HashMap<String, GffEntry>,
) -> HashMap<String, Vec<(usize, usize, String)>> {
    let mut intervals: HashMap<String, Vec<(usize, usize, String)>> = HashMap::new();
    for (name, entry) in gff {
        intervals
            .entry(entry.scaffold.clone())
            .or_default()
            .push((entry.start, entry.end, name.clone()));
    }
    for v in intervals.values_mut() {
        v.sort();
    }
    intervals
}

/// Find the end of the nearest exon LEFT of `position` on this scaffold.
/// Intervals are sorted by start, so binary search for the insertion point.
fn find_left_bound(
    scaffold: &str,
    position: usize,
    intervals: &HashMap<String, Vec<(usize, usize, String)>>,
) -> usize {
    if let Some(ivs) = intervals.get(scaffold) {
        // Find the last interval whose start < position
        let idx = ivs.partition_point(|&(s, _, _)| s < position);
        // Walk backwards from idx to find the best end < position
        let mut best = 0usize;
        for i in (0..idx).rev() {
            if ivs[i].1 < position {
                best = best.max(ivs[i].1);
                break; // sorted by start, so earlier intervals have smaller ends in practice
            }
        }
        best
    } else {
        0
    }
}

/// Find the start of the nearest exon RIGHT of `position` on this scaffold.
/// Intervals are sorted by start, so binary search for first start > position.
fn find_right_bound(
    scaffold: &str,
    position: usize,
    intervals: &HashMap<String, Vec<(usize, usize, String)>>,
    scaffold_len: usize,
) -> usize {
    if let Some(ivs) = intervals.get(scaffold) {
        let idx = ivs.partition_point(|&(s, _, _)| s <= position);
        if idx < ivs.len() {
            return ivs[idx].0;
        }
    }
    scaffold_len
}

/// Count reference gap columns in [gap_start..gap_end).
fn motif_count_ref_gaps(
    gap_start: usize,
    gap_end: usize,
    ref_consensus: &HashMap<usize, Vec<u8>>,
    ref_gap_thresh: f64,
) -> (HashSet<usize>, Vec<usize>, usize) {
    let mut ref_gaps = HashSet::new();
    let mut insert_at = Vec::new();
    let mut longest_consecutive: usize = 0;
    let mut consecutive_non_gap: usize = 0;
    for (i, x) in (gap_start..gap_end).enumerate() {
        let col = ref_consensus.get(&x);
        let is_gap = if let Some(col) = col {
            let gap_count = col.iter().filter(|&&c| c == b'-' || c == b' ').count();
            (gap_count as f64 / col.len() as f64) >= ref_gap_thresh
        } else {
            true
        };
        if is_gap {
            ref_gaps.insert(x);
            insert_at.push(i);
            longest_consecutive = longest_consecutive.max(consecutive_non_gap);
            consecutive_non_gap = 0;
        } else {
            consecutive_non_gap += 1;
        }
    }
    longest_consecutive = longest_consecutive.max(consecutive_non_gap);
    (ref_gaps, insert_at, longest_consecutive)
}

/// Filter reference consensus columns, trimming low-coverage edges.
fn motif_filter_ref_consensus(
    ref_count: usize,
    ref_consensus: &HashMap<usize, Vec<u8>>,
    ref_gaps: &HashSet<usize>,
    gap_start: usize,
    gap_end: usize,
    coverage_thresh: f64,
    left_data_cols: Option<f64>,
    right_data_cols: Option<f64>,
) -> (HashMap<usize, Vec<u8>>, HashSet<usize>) {
    let mut trimmed_cols = HashSet::new();

    if let Some(thresh) = left_data_cols {
        for x in gap_start..gap_end {
            if let Some(col) = ref_consensus.get(&x) {
                let gap_count = col.iter().filter(|&&c| c == b'-' || c == b' ').count();
                let data_present = 1.0 - (gap_count as f64 / ref_count as f64);
                if data_present < thresh {
                    trimmed_cols.insert(x);
                } else {
                    break;
                }
            } else {
                trimmed_cols.insert(x);
            }
        }
    }

    if let Some(thresh) = right_data_cols {
        for x in (gap_start..gap_end).rev() {
            if let Some(col) = ref_consensus.get(&x) {
                let gap_count = col.iter().filter(|&&c| c == b'-' || c == b' ').count();
                let data_present = 1.0 - (gap_count as f64 / ref_count as f64);
                if data_present < thresh {
                    trimmed_cols.insert(x);
                } else {
                    break;
                }
            } else {
                trimmed_cols.insert(x);
            }
        }
    }

    let total_skipped = trimmed_cols.len() + ref_gaps.len();
    if total_skipped >= gap_end - gap_start {
        return (HashMap::new(), trimmed_cols);
    }

    // Build per-ref sequences to find coverage, then filter refs
    let mut ref_coverages: Vec<f64> = Vec::with_capacity(ref_count);
    for y in 0..ref_count {
        let mut data = 0usize;
        let mut total = 0usize;
        for x in gap_start..gap_end {
            if ref_gaps.contains(&x) || trimmed_cols.contains(&x) {
                continue;
            }
            total += 1;
            if let Some(col) = ref_consensus.get(&x) {
                if y < col.len() && col[y] != b'-' && col[y] != b' ' {
                    data += 1;
                }
            }
        }
        ref_coverages.push(if total > 0 { data as f64 / total as f64 } else { 0.0 });
    }

    let max_cov = ref_coverages.iter().cloned().fold(0.0f64, f64::max);
    let target = coverage_thresh * max_cov;
    let ref_indices: Vec<usize> = ref_coverages
        .iter()
        .enumerate()
        .filter(|(_, &c)| c >= target)
        .map(|(i, _)| i)
        .collect();

    let mut this_consensus: HashMap<usize, Vec<u8>> = HashMap::new();
    for x in gap_start..gap_end {
        if ref_gaps.contains(&x) || trimmed_cols.contains(&x) {
            continue;
        }
        if let Some(col) = ref_consensus.get(&x) {
            let filtered: Vec<u8> = ref_indices
                .iter()
                .filter_map(|&y| {
                    if y < col.len() && col[y] != b' ' && col[y] != b'-' {
                        Some(col[y])
                    } else {
                        None
                    }
                })
                .collect();
            this_consensus.insert(x, filtered);
        }
    }

    (this_consensus, trimmed_cols)
}

/// Translate nucleotide bytes to amino acid string (standard code).
fn motif_translate(nt: &[u8]) -> Vec<u8> {
    let mut aa = Vec::with_capacity(nt.len() / 3);
    let mut i = 0;
    while i + 2 < nt.len() {
        let c = codon_to_aa(
            nt[i].to_ascii_uppercase(),
            nt[i + 1].to_ascii_uppercase(),
            nt[i + 2].to_ascii_uppercase(),
        );
        aa.push(if c == 0 { b'X' } else { c });
        i += 3;
    }
    aa
}

/// Trim sequence to multiple of 3.
fn trim_to_codon(seq: &[u8]) -> &[u8] {
    let r = seq.len() % 3;
    if r > 0 { &seq[..seq.len() - r] } else { seq }
}

/// Insert gap dashes into an amino acid or nucleotide string at given positions.
fn insert_gaps(input: &str, positions: &[usize], offset: usize, is_nt: bool) -> String {
    if positions.is_empty() {
        return input.to_string();
    }
    if is_nt {
        let mut triplets: Vec<String> = input
            .as_bytes()
            .chunks(3)
            .map(|c| String::from_utf8_lossy(c).to_string())
            .collect();
        for &pos in positions {
            let idx = offset + pos;
            if idx <= triplets.len() {
                triplets.insert(idx, "---".to_string());
            }
        }
        triplets.join("")
    } else {
        let mut chars: Vec<char> = input.chars().collect();
        for &pos in positions {
            let idx = offset + pos;
            if idx <= chars.len() {
                chars.insert(idx, '-');
            }
        }
        chars.into_iter().collect()
    }
}

/// Dense per-column count table for O(1) amino acid count lookups.
/// Uses a flat array [256] per column to avoid HashMap overhead in the hot loop.
struct DenseColTable {
    /// For each ref_col index (0..ref_cols.len()), a 256-byte count array.
    counts: Vec<[u16; 256]>,
    /// Per ref_col index: max count across all AAs, capped at max_score.
    col_max: Vec<usize>,
    /// The max_score cap used when building.
    max_score: usize,
}

impl DenseColTable {
    fn build(
        consensus: &HashMap<usize, Vec<u8>>,
        ref_cols: &[usize],
        max_score: usize,
    ) -> Self {
        let n = ref_cols.len();
        let mut counts = vec![[0u16; 256]; n];
        let mut col_max = vec![0usize; n];
        for (i, &col_idx) in ref_cols.iter().enumerate() {
            if let Some(col) = consensus.get(&col_idx) {
                for &a in col {
                    counts[i][a as usize] += 1;
                }
                let mx = counts[i].iter().copied().max().unwrap_or(0) as usize;
                col_max[i] = mx.min(max_score);
            }
        }
        DenseColTable { counts, col_max, max_score }
    }

    /// O(1) lookup: count of `aa` in ref_col at dense index `i`.
    #[inline(always)]
    fn count(&self, i: usize, aa: u8) -> usize {
        (self.counts[i][aa as usize] as usize).min(self.max_score)
    }

    /// O(1) lookup: per-column max at dense index `i`.
    #[inline(always)]
    fn max(&self, i: usize) -> usize {
        self.col_max[i]
    }
}

/// Scan result: no kmer allocation, just indices.
struct ScanResult {
    score: usize,
    qstart: usize,
    qend: usize,
    frame: usize,
}

/// Pre-translated protein sequences for all 3 frames.
struct FrameTranslations {
    proteins: [Vec<u8>; 3],
}

impl FrameTranslations {
    fn build(seq: &[u8]) -> Self {
        let mut proteins = [Vec::new(), Vec::new(), Vec::new()];
        for frame in 0..3 {
            if frame < seq.len() {
                proteins[frame] = motif_translate(trim_to_codon(&seq[frame..]));
            }
        }
        FrameTranslations { proteins }
    }

    /// Reconstruct kmer bytes from a scan result.
    fn kmer(&self, r: &ScanResult) -> &[u8] {
        let protein_i = (r.qstart - r.frame) / 3;
        let kmer_len = (r.qend - r.qstart) / 3;
        &self.proteins[r.frame][protein_i..protein_i + kmer_len]
    }
}

/// Core kmer scanning across 3 reading frames.
/// Returns index-only results (no kmer allocation) plus precomputed tables.
fn scan_kmer(
    amount: usize,
    log: &mut Vec<String>,
    splice_region: &[u8],
    skip_cols: &HashSet<usize>,
    flex: usize,
    consensus: &HashMap<usize, Vec<u8>>,
    gap_start: usize,
    gap_end: usize,
    max_score: usize,
    stop_penalty: usize,
) -> (usize, Vec<ScanResult>, Vec<usize>, DenseColTable, FrameTranslations) {
    let ref_cols: Vec<usize> = (gap_start..gap_end)
        .filter(|x| !skip_cols.contains(x))
        .collect();
    let kmer_size = (amount / 3).saturating_sub(skip_cols.len()).saturating_sub(flex);

    let ft = FrameTranslations::build(splice_region);

    if kmer_size == 0 || ref_cols.is_empty() {
        let empty_ct = DenseColTable { counts: Vec::new(), col_max: Vec::new(), max_score };
        return (0, Vec::new(), ref_cols, empty_ct, ft);
    }

    let ct = DenseColTable::build(consensus, &ref_cols, max_score);

    let highest_possible_score: usize = (0..ref_cols.len()).map(|i| ct.max(i)).sum();

    let mut results: Vec<ScanResult> = Vec::new();

    for frame in 0..3usize {
        let protein_seq = &ft.proteins[frame];
        if protein_seq.len() <= kmer_size {
            continue;
        }

        for i in 0..protein_seq.len() - kmer_size {
            let kmer = &protein_seq[i..i + kmer_size];
            for offset in 0..=flex {
                if offset + kmer_size > ref_cols.len() {
                    continue;
                }
                let mut kmer_score: usize = 0;
                for ki in 0..kmer_size {
                    kmer_score += ct.count(ki + offset, kmer[ki]);
                }
                let stops = kmer.iter().filter(|&&c| c == b'*').count();
                kmer_score = kmer_score.saturating_sub(stop_penalty * stops);

                let best_qstart = (i * 3) + frame;
                let best_qend = best_qstart + (kmer_size * 3);
                results.push(ScanResult {
                    score: kmer_score,
                    qstart: best_qstart,
                    qend: best_qend,
                    frame,
                });
            }
        }
    }

    log.push(format!(
        "Translated splice ({} bp, {} ref cols, kmer_size={})",
        splice_region.len(),
        ref_cols.len(),
        kmer_size
    ));

    (highest_possible_score, results, ref_cols, ct, ft)
}

/// Expand a kmer hit by up to `max_expand` AA in each direction.
/// Uses pre-translated frames and precomputed dense column table.
fn expand_hit(
    kmer: &[u8],
    score: usize,
    qstart: usize,
    frame: usize,
    protein_seq: &[u8],
    ref_cols: &[usize],
    ct: &DenseColTable,
    max_score: usize,
    stop_penalty: usize,
    max_expand: usize,
) -> (Vec<u8>, usize, usize, usize, usize, f64) {
    let kmer_size = kmer.len();
    let protein_i = (qstart - frame) / 3;

    // Find original offset
    let raw_score = score + stop_penalty * kmer.iter().filter(|&&c| c == b'*').count();
    let mut best_offset = 0;
    for offset in 0..2.min(ref_cols.len().saturating_sub(kmer_size) + 1) {
        if offset + kmer_size > ref_cols.len() {
            continue;
        }
        let sc: usize = (0..kmer_size)
            .map(|k| ct.count(offset + k, kmer[k]))
            .sum();
        if sc == raw_score {
            best_offset = offset;
            break;
        }
    }

    let orig_mx: usize = (0..kmer_size).map(|k| ct.max(best_offset + k)).sum();
    let orig_ratio = if orig_mx > 0 { score as f64 / orig_mx as f64 } else { 0.0 };

    let mut best_ratio = orig_ratio;
    let mut best_result = (kmer_size, score, qstart, qstart + kmer_size * 3, frame, protein_i, best_offset);

    for left in 0..=max_expand {
        for right in 0..=max_expand {
            if left == 0 && right == 0 {
                continue;
            }
            if protein_i < left || best_offset < left {
                continue;
            }
            let ni = protein_i - left;
            let ns = kmer_size + left + right;
            let no = best_offset - left;

            if ni + ns > protein_seq.len() || no + ns > ref_cols.len() {
                continue;
            }

            let expanded = &protein_seq[ni..ni + ns];
            let mut sc: usize = 0;
            let mut mx: usize = 0;
            for k in 0..ns {
                sc += ct.count(no + k, expanded[k]);
                mx += ct.max(no + k);
            }
            let stops = expanded.iter().filter(|&&c| c == b'*').count();
            sc = sc.saturating_sub(stop_penalty * stops);

            let ratio = if mx > 0 { sc as f64 / mx as f64 } else { 0.0 };
            if ratio > best_ratio {
                best_ratio = ratio;
                let new_qstart = ni * 3 + frame;
                best_result = (ns, sc, new_qstart, new_qstart + ns * 3, frame, ni, no);
            }
        }
    }

    let (_, b_score, b_qstart, b_qend, b_frame, b_ni, _) = best_result;
    let b_len = (b_qend - b_qstart) / 3;
    let result_kmer = protein_seq[b_ni..b_ni + b_len].to_vec();

    (result_kmer, b_score, b_qstart, b_qend, b_frame, best_ratio)
}

/// Scan the trailing (right) flank of the last node in a cluster.
fn motif_scan_flank(
    gap_start: usize,
    gap_end: usize,
    node: &MotifNode,
    is_leading: bool,
    id_check: &mut HashSet<String>,
    id_count: &mut HashMap<usize, usize>,
    log: &mut Vec<String>,
    ref_consensus: &HashMap<usize, Vec<u8>>,
    ref_count: usize,
    gff: &HashMap<String, GffEntry>,
    genome: &HashMap<String, Vec<u8>>,
    scaffold_intervals: &HashMap<String, Vec<(usize, usize, String)>>,
    last_node_seq_len: usize,
    last_node_nt_seq_len: usize,
    results: &mut Vec<MotifHit>,
) {
    let amount = (gap_end - gap_start) * 3;
    let direction = if is_leading { "Left leading" } else { "Right trailing" };
    log.push(format!(
        "{} gap of {} with size {}",
        direction, node.header, amount
    ));

    const MINIMUM_GAP_BP: usize = 15;
    const MAX_GAP_BP: usize = 750;
    const REF_GAP_THR: f64 = 0.7;
    const LR_REF_COV: f64 = 0.8;
    const MIN_CONSEC: usize = 5;
    const MINIMUM_AA: usize = 5;
    const REQ_END_DATA_COLS: f64 = 0.75;
    const FLANK_BP: usize = 20000;
    const MAX_SCORE: usize = 100;
    const STOP_PENALTY: usize = 5;
    const FLEX: usize = 1;

    if amount < MINIMUM_GAP_BP || amount >= MAX_GAP_BP {
        log.push("Gap too small or too large\n".to_string());
        return;
    }

    let node_field = parse_node_field(&node.header).replace("NODE_", "");
    let node_name = if node_field.contains("&&") {
        if is_leading {
            node_field.split("&&").next().unwrap_or(&node_field)
        } else {
            node_field.split("&&").last().unwrap_or(&node_field)
        }
        .split('_')
        .next()
        .unwrap_or(&node_field)
    } else {
        node_field.split('_').next().unwrap_or(&node_field)
    };

    let gff_entry = gff.get(&node_field).or_else(|| gff.get(node_name));
    let gff_entry = match gff_entry {
        Some(e) => e,
        None => {
            log.push(format!("No GFF entry for {}\n", node_field));
            return;
        }
    };

    let scaffold_seq = match genome.get(&gff_entry.scaffold) {
        Some(s) => s,
        None => {
            log.push(format!("Scaffold {} not loaded\n", gff_entry.scaffold));
            return;
        }
    };
    let scaffold_len = scaffold_seq.len();

    let (seq, rs, re) = if is_leading {
        // Upstream of this exon
        if gff_entry.strand == "+" {
            let left = find_left_bound(&gff_entry.scaffold, gff_entry.start, scaffold_intervals);
            let s = left.max(gff_entry.start.saturating_sub(1).saturating_sub(FLANK_BP));
            let e = gff_entry.start.saturating_sub(1);
            (scaffold_seq[s..e].to_vec(), s, e)
        } else {
            let right = find_right_bound(
                &gff_entry.scaffold,
                gff_entry.end,
                scaffold_intervals,
                scaffold_len,
            );
            let s = gff_entry.end;
            let e = right.min(gff_entry.end + FLANK_BP);
            (bio_revcomp(&scaffold_seq[s..e]), s, e)
        }
    } else {
        // Downstream of this exon
        if gff_entry.strand == "+" {
            let right = find_right_bound(
                &gff_entry.scaffold,
                gff_entry.end,
                scaffold_intervals,
                scaffold_len,
            );
            let s = gff_entry.end;
            let e = right.min(gff_entry.end + FLANK_BP);
            (scaffold_seq[s..e].to_vec(), s, e)
        } else {
            let left = find_left_bound(&gff_entry.scaffold, gff_entry.start, scaffold_intervals);
            let s = left.max(gff_entry.start.saturating_sub(1).saturating_sub(FLANK_BP));
            let e = gff_entry.start.saturating_sub(1);
            (bio_revcomp(&scaffold_seq[s..e]), s, e)
        }
    };

    if seq.is_empty() {
        log.push("Empty flank region\n".to_string());
        return;
    }

    log.push(format!(
        "  Node: {}:{}-{} ({}) Search: {}:{}-{} len={}",
        gff_entry.scaffold,
        gff_entry.start,
        gff_entry.end,
        gff_entry.strand,
        gff_entry.scaffold,
        rs + 1,
        re,
        seq.len()
    ));

    let (ref_gaps, insert_at, longest_consecutive) =
        motif_count_ref_gaps(gap_start, gap_end, ref_consensus, REF_GAP_THR);
    if longest_consecutive < MIN_CONSEC {
        log.push(format!("Longest consecutive ref cols {} < {}\n", longest_consecutive, MIN_CONSEC));
        return;
    }

    let (left_dc, right_dc) = if is_leading {
        (Some(REQ_END_DATA_COLS), None)
    } else {
        (None, Some(REQ_END_DATA_COLS))
    };

    let (this_consensus, edge_trim_cols) = motif_filter_ref_consensus(
        ref_count,
        ref_consensus,
        &ref_gaps,
        gap_start,
        gap_end,
        LR_REF_COV,
        left_dc,
        right_dc,
    );
    if this_consensus.is_empty() {
        log.push("No consensus columns after filtering\n".to_string());
        return;
    }

    let skip_cols: HashSet<usize> = edge_trim_cols.union(&ref_gaps).cloned().collect();

    let (highest_possible_score, mut kmer_results, ref_cols, col_ct, frame_trans) = scan_kmer(
        amount,
        log,
        &seq,
        &skip_cols,
        FLEX,
        &this_consensus,
        gap_start,
        gap_end,
        MAX_SCORE,
        STOP_PENALTY,
    );

    if kmer_results.is_empty() {
        log.push("No suitable kmer found\n".to_string());
        return;
    }

    // Partial sort: only need the top candidates, not a full sort of ~40k entries
    let top_n = kmer_results.len().min(200);
    kmer_results.select_nth_unstable_by(top_n.saturating_sub(1), |a, b| b.score.cmp(&a.score));
    kmer_results.truncate(top_n);
    kmer_results.sort_unstable_by(|a, b| b.score.cmp(&a.score));

    let mut accepted = 0usize;
    let mut failed_logged = 0usize;
    let mut seen_positions: HashSet<(usize, usize)> = HashSet::new();

    for r in &kmer_results {
        if accepted >= 1 {
            break;
        }
        let kmer = frame_trans.kmer(r);
        if kmer.len() < MINIMUM_AA {
            continue;
        }

        let pos_key = (r.frame, r.qstart / 9);
        if seen_positions.contains(&pos_key) {
            continue;
        }

        let protein_seq = &frame_trans.proteins[r.frame];
        let (exp_kmer, exp_score, exp_qstart, exp_qend, exp_frame, _exp_ratio) = expand_hit(
            kmer, r.score, r.qstart, r.frame, protein_seq, &ref_cols, &col_ct, MAX_SCORE,
            STOP_PENALTY, 10,
        );

        let exp_size = exp_kmer.len();
        let n_cols = ref_cols.len().min(exp_size);
        let exp_max: usize = (0..n_cols).map(|i| col_ct.max(i)).sum();

        let threshold: f64 = if exp_size >= 10 { 0.5 } else { 0.85 };
        let threshold_score = (exp_max as f64 * threshold).round() as usize;
        if exp_max > 0 && exp_score < threshold_score {
            if accepted == 0 && failed_logged < 10 {
                log.push(format!(
                    "Failed score threshold: {} ({:.0}%) < {:.0}%",
                    exp_score,
                    exp_score as f64 / exp_max as f64 * 100.0,
                    threshold * 100.0
                ));
                failed_logged += 1;
            }
            continue;
        }

        seen_positions.insert(pos_key);

        // Finalise: build header, aa sequence, nt sequence
        let exp_kmer_str = String::from_utf8_lossy(&exp_kmer).to_string();

        let node_id: usize = node_name.parse().unwrap_or(0);

        let count = id_count.entry(node_id).or_insert(0);
        let mut new_node = if *count == 0 {
            node_id.to_string()
        } else {
            format!("{}_{}", node_id, count)
        };
        while id_check.contains(&new_node) {
            *count += 1;
            new_node = if *count == 0 {
                node_id.to_string()
            } else {
                format!("{}_{}", node_id, count)
            };
        }
        *count += 1;
        id_check.insert(new_node.clone());

        let mut best_frame_val = (exp_frame as i32) + 1;
        if node.frame < 0 {
            best_frame_val = -best_frame_val;
        }

        let mut header_fields: Vec<&str> = node.header.split('|').collect();
        let node_field_owned = format!("NODE_{}", new_node);
        let frame_str = best_frame_val.to_string();
        let score_str = "1".to_string();
        if header_fields.len() > 5 {
            header_fields[3] = &node_field_owned;
            header_fields[4] = &frame_str;
            header_fields[5] = &score_str;
        }
        let new_header = header_fields.join("|");

        let gapped_kmer = insert_gaps(&exp_kmer_str, &insert_at, 0, false);
        let mut new_aa = "-".repeat(gap_start) + &gapped_kmer;
        if new_aa.len() < last_node_seq_len {
            new_aa.push_str(&"-".repeat(last_node_seq_len - new_aa.len()));
        }

        let nt_seq_raw = String::from_utf8_lossy(&seq[exp_qstart..exp_qend]).to_string();
        let gapped_nt = insert_gaps(&nt_seq_raw, &insert_at, 0, true);
        let mut new_nt = "-".repeat(gap_start * 3) + &gapped_nt;
        if new_nt.len() < last_node_nt_seq_len {
            new_nt.push_str(&"-".repeat(last_node_nt_seq_len - new_nt.len()));
        }

        // Compute genomic coordinates for GFF
        let strand = &gff_entry.strand;
        let (hit_start, hit_end) = if strand == "+" {
            (rs + exp_qstart + 1, rs + exp_qend)
        } else {
            let slen = seq.len();
            (rs + slen - exp_qend + 1, rs + slen - exp_qstart)
        };
        let region = format!("{}:{}-{}({})", gff_entry.scaffold, hit_start, hit_end, strand);

        log.push(format!(
            "Best match: {} score={} max={} threshold={:.0}%",
            exp_kmer_str, exp_score, highest_possible_score, threshold * 100.0
        ));
        log.push(format!("Inserted sequence: {}", new_header));
        log.push(format!("nt={}\n", nt_seq_raw));

        results.push(MotifHit {
            header: new_header,
            aa_seq: new_aa,
            nt_seq: new_nt,
            region,
            score: exp_score,
            node_name: node_name.to_string(),
            is_leading,
        });
        accepted += 1;
    }

    log.push(String::new());
}

/// Find the start and end of non-gap content in a sequence.
fn motif_find_index_pair(seq: &str) -> (usize, usize) {
    let bytes = seq.as_bytes();
    let mut start = 0;
    let mut end = bytes.len();
    for (i, &b) in bytes.iter().enumerate() {
        if b != b'-' {
            start = i;
            break;
        }
    }
    for i in (0..bytes.len()).rev() {
        if bytes[i] != b'-' {
            end = i + 1;
            break;
        }
    }
    (start, end)
}

/// Run motif scan on a single gene.
fn motif_scan_gene(
    gene: &str,
    aa_entries: &[(String, String)],
    clusters_map: &HashMap<String, HashMap<String, (Vec<String>, String)>>,
    gff: &HashMap<String, GffEntry>,
    genome: &HashMap<String, Vec<u8>>,
    scaffold_intervals: &HashMap<String, Vec<(usize, usize, String)>>,
) -> (Vec<String>, Vec<MotifHit>) {
    let mut log: Vec<String> = Vec::new();
    let mut all_hits: Vec<MotifHit> = Vec::new();

    log.push(format!("=== Motif scan: {} ===", gene));

    // Build ref consensus and node list from AA entries
    let mut ref_consensus: HashMap<usize, Vec<u8>> = HashMap::new();
    let mut ref_starts: Vec<usize> = Vec::new();
    let mut ref_ends: Vec<usize> = Vec::new();
    let mut ref_count: usize = 0;
    let mut aa_nodes: Vec<MotifNode> = Vec::new();

    for (header, seq) in aa_entries {
        if header.ends_with('.') {
            ref_count += 1;
            let (start, end) = motif_find_index_pair(seq);
            ref_starts.push(start);
            ref_ends.push(end);
            for (i, b) in seq.bytes().enumerate() {
                let col = ref_consensus.entry(i).or_default();
                if i >= start && i <= end {
                    col.push(b);
                } else {
                    col.push(b'-');
                }
            }
            continue;
        }
        let fields: Vec<&str> = header.split('|').collect();
        let frame: i32 = fields.get(4).and_then(|s| s.parse().ok()).unwrap_or(1);
        let (start, end) = motif_find_index_pair(seq);
        aa_nodes.push(MotifNode {
            header: header.clone(),
            frame,
            sequence: seq.clone(),
            start,
            end,
        });
    }

    if aa_nodes.is_empty() || ref_count == 0 || ref_starts.is_empty() {
        return (log, all_hits);
    }

    // Ref median start/end (use select_nth_unstable instead of clone+sort)
    let median_idx = ref_starts.len() / 2;
    let ref_median_start = *ref_starts.select_nth_unstable(median_idx).1;
    let ref_median_end = *ref_ends.select_nth_unstable(median_idx).1;

    // Build id tracking
    let mut id_count: HashMap<usize, usize> = HashMap::new();
    let mut id_check: HashSet<String> = HashSet::new();
    for node in &aa_nodes {
        let nf = parse_node_field(&node.header).replace("NODE_", "");
        for part in nf.split("&&") {
            id_check.insert(part.to_string());
        }
        for part in nf.split("&&") {
            if let Ok(id) = part.split('_').next().unwrap_or("").parse::<usize>() {
                *id_count.entry(id).or_insert(0) += 1;
            }
        }
    }

    // Get cluster sets for this gene
    let clusters = match clusters_map.get(gene) {
        Some(c) => c,
        None => return (log, all_hits),
    };

    // Convert cluster map to sets of NODE_ prefixed tokens, keeping the cluster key
    let cluster_sets: Vec<(&str, HashSet<String>)> = clusters
        .iter()
        .map(|(key, (nodes, _))| {
            (
                key.as_str(),
                nodes.iter().map(|n| format!("NODE_{}", n)).collect(),
            )
        })
        .collect();

    // Get reference sequence length for padding
    let last_seq_len = aa_entries
        .iter()
        .filter(|(h, _)| !h.ends_with('.'))
        .map(|(_, s)| s.len())
        .max()
        .unwrap_or(0);
    let last_nt_seq_len = last_seq_len * 3;

    for (_cluster_key, cluster_set) in &cluster_sets {
        let mut aa_subset: Vec<&MotifNode> = aa_nodes
            .iter()
            .filter(|n| {
                let nf = n.header.split('|').nth(3).unwrap_or("");
                cluster_set.contains(nf)
            })
            .collect();

        if aa_subset.is_empty() {
            continue;
        }

        if aa_subset[0].frame < 0 {
            aa_subset.sort_by(|a, b| b.start.cmp(&a.start));
        } else {
            aa_subset.sort_by(|a, b| a.start.cmp(&b.start));
        }

        let (first_node, last_node) = if aa_subset[0].frame < 0 {
            (
                aa_subset[aa_subset.len() - 1],
                aa_subset[0],
            )
        } else {
            (aa_subset[0], aa_subset[aa_subset.len() - 1])
        };

        // Leading gap: from ref_median_start to first_node.start
        if first_node.start > ref_median_start {
            motif_scan_flank(
                ref_median_start,
                first_node.start,
                first_node,
                true,
                &mut id_check,
                &mut id_count,
                &mut log,
                &ref_consensus,
                ref_count,
                gff,
                genome,
                scaffold_intervals,
                last_seq_len,
                last_nt_seq_len,
                &mut all_hits,
            );
        }

        // Trailing gap: from last_node.end to ref_median_end
        if ref_median_end > last_node.end {
            motif_scan_flank(
                last_node.end,
                ref_median_end,
                last_node,
                false,
                &mut id_check,
                &mut id_count,
                &mut log,
                &ref_consensus,
                ref_count,
                gff,
                genome,
                scaffold_intervals,
                last_seq_len,
                last_nt_seq_len,
                &mut all_hits,
            );
        }
    }

    (log, all_hits)
}

// ---------------------------------------------------------------------------
// End motif scan
// ---------------------------------------------------------------------------

struct GffEntry {
    scaffold: String,
    start: usize,
    end: usize,
    strand: String,
}

fn collect_needed_scaffolds(
    gene_clusters: &HashMap<String, HashMap<String, (Vec<String>, String)>>,
    gff: &HashMap<String, GffEntry>,
) -> HashSet<String> {
    let mut needed = HashSet::new();
    for (_gene, clusters) in gene_clusters {
        for (_ck, (node_tokens, _iso_type)) in clusters {
            for token in node_tokens {
                if let Some(entry) = gff.get(token.as_str()) {
                    needed.insert(entry.scaffold.clone());
                }
            }
        }
    }
    needed
}

fn load_genome_from_rocksdb(
    rocksdb_path: &str,
    needed: &HashSet<String>,
) -> HashMap<String, Vec<u8>> {
    let mut opts = Options::default();
    opts.set_error_if_exists(false);
    let db = match DB::open_for_read_only(&opts, rocksdb_path, false) {
        Ok(db) => db,
        Err(e) => {
            eprintln!("Failed to open RocksDB at {}: {}", rocksdb_path, e);
            return HashMap::new();
        }
    };

    // Try indexed lookups first (new format: scaffold_index key)
    if let Ok(Some(idx_raw)) = db.get(b"scaffold_index") {
        let idx_str = String::from_utf8_lossy(&idx_raw);
        // Parse index: each line is "name\tbatch_index\tbyte_offset\tseq_length"
        // Group needed scaffolds by batch so we only load each batch once.
        let mut batch_requests: HashMap<String, Vec<(String, usize, usize)>> = HashMap::new();
        for line in idx_str.lines() {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 4 {
                continue;
            }
            let name = parts[0];
            if !needed.contains(name) {
                continue;
            }
            let batch_id = parts[1].to_string();
            let offset: usize = parts[2].parse().unwrap_or(0);
            let length: usize = parts[3].parse().unwrap_or(0);
            batch_requests
                .entry(batch_id)
                .or_default()
                .push((name.to_string(), offset, length));
        }

        let mut scaffolds = HashMap::with_capacity(needed.len());
        for (batch_id, entries) in &batch_requests {
            let key = format!("parentbatch:{}", batch_id);
            let data = match db.get(key.as_bytes()) {
                Ok(Some(v)) => v,
                _ => continue,
            };
            for (name, offset, length) in entries {
                let end = (*offset + *length).min(data.len());
                let start = (*offset).min(end);
                scaffolds.insert(name.clone(), data[start..end].to_vec());
            }
        }
        return scaffolds;
    }

    // Fallback: legacy batched format (scan all batches)
    let recipe = match db.get(b"getall:parentbatches") {
        Ok(Some(v)) => String::from_utf8_lossy(&v).to_string(),
        _ => {
            eprintln!("No parentbatches recipe in RocksDB");
            return HashMap::new();
        }
    };

    let batch_keys: Vec<&str> = recipe.split(',').collect();
    let mut scaffolds = HashMap::new();

    for bk in &batch_keys {
        let key = format!("parentbatch:{}", bk);
        let data = match db.get(key.as_bytes()) {
            Ok(Some(v)) => String::from_utf8_lossy(&v).to_string(),
            _ => continue,
        };
        let mut pos = 0;
        let bytes = data.as_bytes();
        let dlen = bytes.len();
        while pos < dlen {
            let nl = match data[pos..].find('\n') {
                Some(p) => pos + p,
                None => break,
            };
            let name = &data[pos + 1..nl];
            pos = nl + 1;
            let nl2 = match data[pos..].find('\n') {
                Some(p) => pos + p,
                None => {
                    if needed.contains(name) {
                        scaffolds.insert(name.to_string(), data[pos..].as_bytes().to_vec());
                    }
                    break;
                }
            };
            if needed.contains(name) {
                scaffolds.insert(name.to_string(), data[pos..nl2].as_bytes().to_vec());
            }
            pos = nl2 + 1;
        }
    }
    scaffolds
}

fn read_gff(path: &str) -> HashMap<String, GffEntry> {
    let mut nodes = HashMap::new();
    let data = match fs::read_to_string(path) {
        Ok(data) => data,
        Err(_) => return nodes,
    };
    for line in data.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 {
            continue;
        }
        let scaffold = parts[0].to_string();
        let start: usize = parts[3].parse().unwrap_or(0);
        let end: usize = parts[4].parse().unwrap_or(0);
        let strand = parts[6].to_string();
        let attrs = parts[8];
        let mut name = None;
        for attr in attrs.split(';') {
            if let Some(v) = attr.strip_prefix("Name=") {
                name = Some(v.to_string());
                break;
            }
        }
        if let Some(n) = name {
            nodes.insert(
                n,
                GffEntry {
                    scaffold,
                    start,
                    end,
                    strand,
                },
            );
        }
    }
    nodes
}

fn read_clusters(path: &str) -> HashMap<String, HashMap<String, (Vec<String>, String)>> {
    let mut gene_clusters: HashMap<String, HashMap<String, (Vec<String>, String)>> = HashMap::new();
    let data = match fs::read_to_string(path) {
        Ok(data) => data,
        Err(_) => return gene_clusters,
    };
    let mut first = true;
    for line in data.lines() {
        if first {
            first = false;
            continue;
        }
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 7 {
            continue;
        }
        let gene = fields[0].trim_end_matches(".gz").to_string();
        let cluster_key = fields[1].to_string();
        let nodes: Vec<String> = fields[6]
            .split(';')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect();
        let isoform_type = if fields.len() > 8 {
            fields[8].trim().to_string()
        } else {
            String::new()
        };
        gene_clusters
            .entry(gene)
            .or_default()
            .insert(cluster_key, (nodes, isoform_type));
    }
    gene_clusters
}

fn read_aa_fasta(path: &str) -> Vec<(String, String)> {
    let mut entries = Vec::new();
    let data = if path.ends_with(".gz") {
        let file = match fs::File::open(path) {
            Ok(f) => f,
            Err(_) => return entries,
        };
        let mut decoder = GzDecoder::new(file);
        let mut s = String::new();
        if decoder.read_to_string(&mut s).is_err() {
            return entries;
        }
        s
    } else {
        match fs::read_to_string(path) {
            Ok(data) => data,
            Err(_) => return entries,
        }
    };
    let mut header = String::new();
    let mut seq = String::new();
    for line in data.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if !header.is_empty() {
                entries.push((header.clone(), seq.clone()));
            }
            header = rest.to_string();
            seq.clear();
        } else {
            seq.push_str(line.trim());
        }
    }
    if !header.is_empty() {
        entries.push((header, seq));
    }
    entries
}

struct RecoveredExon {
    gene: String,
    header: String,
    aa_seq: String,
    nt_seq: String,
    region: String,
    score: f64,
    cluster_key: String,
    node_a_name: String,
    node_b_name: String,
}

struct GeneResult {
    log: Vec<String>,
    gaps: usize,
    hits: usize,
    gff_lines: Vec<String>,
    recovered: Vec<RecoveredExon>,
}

fn process_gene(
    aa_file: &str,
    entries: &[(String, String)],
    clusters_map: &HashMap<String, HashMap<String, (Vec<String>, String)>>,
    gff: &HashMap<String, GffEntry>,
    genome: &HashMap<String, Vec<u8>>,
    aa_idx: &[u8; 256],
    ct: &HashMap<[u8; 3], u8>,
    score_thr: f64,
    min_aa: usize,
) -> GeneResult {
    let gene = aa_file.trim_end_matches(".gz");
    let mut flog: Vec<String> = Vec::new();
    let mut gff_lines: Vec<String> = Vec::new();
    let mut recovered: Vec<RecoveredExon> = Vec::new();
    let (mut tgaps, mut thits) = (0usize, 0usize);

    let clusters = match clusters_map.get(gene) {
        Some(clusters) => clusters,
        None => {
            return GeneResult {
                log: flog,
                gaps: 0,
                hits: 0,
                gff_lines,
                recovered,
            };
        }
    };
    let gene_base = gene.split('.').next().unwrap_or(gene);

    #[allow(dead_code)]
    struct Cand {
        header: String,
        seq: Vec<u8>,
        frame: i32,
        start: usize,
        end: usize,
        nf: String,
    }

    let mut rc: HashMap<usize, Vec<u8>> = HashMap::new();
    let mut rcount = 0usize;
    let mut cands: Vec<Cand> = Vec::new();
    let mut cseqs: Vec<(usize, usize, Vec<u8>)> = Vec::new();

    for (header, seq) in entries {
        let sb = seq.as_bytes();
        if header.ends_with('.') {
            rcount += 1;
            let (s, e) = find_index_pair(sb, b'-');
            for (i, &bp) in sb.iter().enumerate() {
                let ent = rc.entry(i).or_default();
                ent.push(if i >= s && i <= e { bp } else { b'-' });
            }
        } else {
            let parts: Vec<&str> = header.split('|').collect();
            if parts.len() > 4 {
                if let Ok(frame) = parts[4].parse::<i32>() {
                    let (s, e) = find_index_pair(sb, b'-');
                    cands.push(Cand {
                        header: header.clone(),
                        seq: sb.to_vec(),
                        frame,
                        start: s,
                        end: e,
                        nf: parse_node_field(header).to_string(),
                    });
                    cseqs.push((s, e, sb.to_vec()));
                }
            }
        }
    }

    if rcount == 0 || cands.is_empty() {
        return GeneResult {
            log: flog,
            gaps: 0,
            hits: 0,
            gff_lines,
            recovered,
        };
    }

    // Build combined consensus (refs + candidates)
    let mut combined: HashMap<usize, Vec<u8>> = HashMap::new();
    for (&col, chars) in &rc {
        combined.insert(col, chars.clone());
    }
    for (cs, ce, cseq) in &cseqs {
        for i in *cs..*ce {
            if i < cseq.len() {
                let aa = cseq[i];
                if aa != b'-' {
                    combined.entry(i).or_default().push(aa);
                }
            }
        }
    }

    flog.push(format!(
        "=== {}: {} refs, {} candidates ===",
        gene, rcount, cands.len()
    ));

    let mut hit_id = 0usize;

    // Track recovered exon positions for MXE narrowing
    let mut base_recovered_gff: HashMap<String, Vec<(String, usize, usize)>> = HashMap::new();

    // Separate clusters into base and MXE, skip N-terminal
    let mut base_clusters: Vec<(&String, &Vec<String>, &String)> = Vec::new();
    let mut mxe_clusters: Vec<(&String, &Vec<String>, &String)> = Vec::new();
    for (ck, (node_tokens, iso_type)) in clusters.iter() {
        if iso_type == "N-terminal" {
            continue;
        }
        if iso_type == "MXE" {
            mxe_clusters.push((ck, node_tokens, iso_type));
        } else {
            base_clusters.push((ck, node_tokens, iso_type));
        }
    }
    base_clusters.sort_by_key(|(ck, _, _)| *ck);
    mxe_clusters.sort_by_key(|(ck, _, _)| *ck);

    let all_clusters: Vec<_> = base_clusters
        .iter()
        .chain(mxe_clusters.iter())
        .cloned()
        .collect();

    // Build a map of known GFF intervals for each MXE cluster group.
    // Key = base cluster key (e.g. "86"), Value = vec of (scaffold, start, end)
    // for every node across all sibling isoforms (including the base).
    // Used to narrow MXE gap search regions so the DP does not scan
    // genomic intervals already occupied by known cassette variants
    // from other isoforms.
    let mxe_sibling_gff: HashMap<String, Vec<(String, usize, usize)>> = {
        // Collect all MXE cluster keys grouped by base key
        let mut groups: HashMap<String, Vec<&Vec<String>>> = HashMap::new();
        for (ck, node_tokens, iso_type) in &all_clusters {
            if *iso_type != "MXE" {
                continue;
            }
            let base_key = ck
                .rsplit_once('_')
                .map(|(k, _)| k.to_string())
                .unwrap_or_else(|| (*ck).clone());
            groups.entry(base_key).or_default().push(node_tokens);
        }
        let mut result: HashMap<String, Vec<(String, usize, usize)>> = HashMap::new();
        for (base_key, token_lists) in &groups {
            let mut intervals = Vec::new();
            for tokens in token_lists {
                for token in *tokens {
                    if let Some(entry) = gff.get(token.as_str()) {
                        intervals.push((
                            entry.scaffold.clone(),
                            entry.start,
                            entry.end,
                        ));
                    }
                }
            }
            result.insert(base_key.clone(), intervals);
        }
        result
    };

    // For each MXE base cluster key, compute which of its nodes are
    // modular (i.e. get swapped out in at least one isoform).
    // A base node is modular if any isoform child does NOT contain it.
    let mxe_base_modular: HashMap<String, HashSet<String>> = {
        let mut result: HashMap<String, HashSet<String>> = HashMap::new();
        // Group isoform node-sets by base key
        let mut iso_sets: HashMap<String, Vec<HashSet<String>>> = HashMap::new();
        let mut base_sets: HashMap<String, HashSet<String>> = HashMap::new();
        for (ck, (node_tokens, iso_type)) in clusters.iter() {
            if iso_type != "MXE" {
                continue;
            }
            let fields: HashSet<String> = node_tokens
                .iter()
                .map(|t| format!("NODE_{}", t))
                .collect();
            if let Some(bk) = ck.rsplit_once('_').map(|(k, _)| k.to_string()) {
                iso_sets.entry(bk).or_default().push(fields);
            } else {
                base_sets.insert(ck.clone(), fields);
            }
        }
        for (bk, base_fields) in &base_sets {
            if let Some(children) = iso_sets.get(bk) {
                // A base node is modular if ANY child lacks it
                let modular: HashSet<String> = base_fields
                    .iter()
                    .filter(|node| children.iter().any(|child| !child.contains(*node)))
                    .cloned()
                    .collect();
                result.insert(bk.clone(), modular);
            }
        }
        result
    };

    for (ck, node_tokens, iso_type) in &all_clusters {
        let cluster_node_fields: HashSet<String> = node_tokens
            .iter()
            .map(|t| format!("NODE_{}", t))
            .collect();

        // For MXE clusters, identify modular nodes.
        // Isoform clusters: nodes present in the isoform but not the base.
        // Base clusters: nodes in the base that are replaced in at
        //   least one isoform child.
        let modular_node_fields: HashSet<String> = if *iso_type == "MXE" {
            if let Some(base_key) = ck.rsplit_once('_').map(|(k, _)| k.to_string()) {
                if let Some((base_tokens, _)) = clusters.get(&base_key) {
                    let base_fields: HashSet<String> = base_tokens
                        .iter()
                        .map(|t| format!("NODE_{}", t))
                        .collect();
                    cluster_node_fields
                        .difference(&base_fields)
                        .cloned()
                        .collect()
                } else {
                    HashSet::new()
                }
            } else {
                // Base cluster: use precomputed modular set
                mxe_base_modular
                    .get(*ck)
                    .cloned()
                    .unwrap_or_default()
            }
        } else {
            HashSet::new()
        };

        // Filter candidates by cluster membership
        let mut aa_subset: Vec<usize> = cands
            .iter()
            .enumerate()
            .filter(|(_, c)| cluster_node_fields.contains(&c.nf))
            .map(|(i, _)| i)
            .collect();

        if aa_subset.len() < 2 {
            continue;
        }

        let is_negative = cands[aa_subset[0]].frame < 0;
        if is_negative {
            aa_subset.sort_by(|&a, &b| cands[b].start.cmp(&cands[a].start));
        } else {
            aa_subset.sort_by(|&a, &b| cands[a].start.cmp(&cands[b].start));
        }

        for wi in 1..aa_subset.len() {
            let na_idx = aa_subset[wi - 1];
            let nb_idx = aa_subset[wi];
            let na = &cands[na_idx];
            let nb = &cands[nb_idx];

            let (left_end, right_start) = if is_negative {
                (nb.end, na.start)
            } else {
                (na.end, nb.start)
            };

            let gap_aa = if right_start > left_end {
                right_start - left_end
            } else {
                0
            };
            let gap_nt = gap_aa * 3;

            if gap_aa == 0
                || gap_nt < _MINIMUM_GAP_BP
                || gap_nt > _MAX_INTERNAL_GAP_BP
            {
                continue;
            }

            let gap_start = left_end;
            let gap_end = right_start;

            // Ref-presence filter: gap columns must have enough
            // consecutive reference-occupied positions to build a
            // usable PSSM.  Gaps in insertion columns are padding.
            {
                let mut longest_consec = 0usize;
                let mut cur_consec = 0usize;
                for col in gap_start..gap_end {
                    let has_ref = match rc.get(&col) {
                        Some(chars) if !chars.is_empty() => {
                            let gf = chars.iter().filter(|&&c| c == b'-').count() as f64
                                / chars.len() as f64;
                            gf < 0.67
                        }
                        _ => false,
                    };
                    if has_ref {
                        cur_consec += 1;
                    } else {
                        if cur_consec > longest_consec {
                            longest_consec = cur_consec;
                        }
                        cur_consec = 0;
                    }
                }
                if cur_consec > longest_consec {
                    longest_consec = cur_consec;
                }
                if longest_consec < MIN_CONSEC_CHAR {
                    continue;
                }
            }

            tgaps += 1;

            let node_a_name = na.nf.replace("NODE_", "");
            let node_b_name = nb.nf.replace("NODE_", "");

            // MXE isoforms: skip gaps between two shared (non-modular)
            // nodes.  These gaps duplicate work done by the base
            // cluster.  Only gaps touching a modular node may contain
            // undiscovered cassette variants in the flanking intron.
            // The base cluster itself must process shared-shared gaps
            // because no other cluster will.
            if *iso_type == "MXE" && ck.contains('_') && !modular_node_fields.is_empty() {
                let a_mod = modular_node_fields.contains(&na.nf);
                let b_mod = modular_node_fields.contains(&nb.nf);
                if !a_mod && !b_mod {
                    continue;
                }
            }

            let ga = match gff.get(&node_a_name) {
                Some(gff_entry) => gff_entry,
                None => continue,
            };
            let gb = match gff.get(&node_b_name) {
                Some(gff_entry) => gff_entry,
                None => continue,
            };

            if ga.scaffold != gb.scaffold {
                continue;
            }
            let strand = &ga.strand;

            let (region_start, region_end) = if ga.start < gb.start {
                (ga.end + 1, if gb.start > 0 { gb.start - 1 } else { 0 })
            } else {
                (gb.end + 1, if ga.start > 0 { ga.start - 1 } else { 0 })
            };

            // MXE gap narrowing.  Two complementary mechanisms:
            //
            // (A) All MXE clusters: narrow around pre-existing GFF
            //     nodes from sibling isoforms so the DP does not
            //     scan over known cassette exons and produce
            //     overlapping hits.
            //
            // (B) Isoform clusters only: further narrow around
            //     DP-recovered exons from the base cluster (which
            //     runs first) so each isoform searches only the
            //     sub-region belonging to its own cassette module.
            let (mut rs, mut re) = (region_start, region_end);
            if *iso_type == "MXE" {
                let base_key = ck
                    .rsplit_once('_')
                    .map(|(k, _)| k.to_string())
                    .unwrap_or_else(|| (*ck).clone());

                // (A) Narrow around sibling GFF nodes
                if let Some(sibling_intervals) = mxe_sibling_gff.get(&base_key) {
                    let foreign: Vec<(usize, usize)> = sibling_intervals
                        .iter()
                        .filter(|(sscaf, ss, se)| {
                            *sscaf == ga.scaffold
                                && *ss >= rs
                                && *se <= re
                                && !cluster_node_fields.iter().any(|nf| {
                                    let bare = nf.strip_prefix("NODE_").unwrap_or(nf);
                                    if let Some(entry) = gff.get(bare) {
                                        entry.scaffold == *sscaf
                                            && entry.start == *ss
                                            && entry.end == *se
                                    } else {
                                        false
                                    }
                                })
                        })
                        .map(|(_, ss, se)| (*ss, *se))
                        .collect();
                    if !foreign.is_empty() {
                        let a_mod = modular_node_fields.contains(&na.nf);
                        let b_mod = modular_node_fields.contains(&nb.nf);
                        let leftmost_foreign = foreign.iter().map(|(s, _)| *s).min().unwrap();
                        let rightmost_foreign_end =
                            foreign.iter().map(|(_, e)| *e).max().unwrap();
                        let mod_is_left = if a_mod && !b_mod {
                            ga.start < gb.start
                        } else if b_mod && !a_mod {
                            gb.start < ga.start
                        } else {
                            let left_size = leftmost_foreign.saturating_sub(rs);
                            let right_size = re.saturating_sub(rightmost_foreign_end);
                            left_size >= right_size
                        };
                        if mod_is_left {
                            re = re.min(leftmost_foreign - 1);
                        } else {
                            rs = rs.max(rightmost_foreign_end + 1);
                        }
                        if rs >= re {
                            continue;
                        }
                    }
                }

                // (B) Isoform clusters: further narrow around base DP hits
                if ck.contains('_') {
                    if let Some(recovered_list) = base_recovered_gff.get(&base_key) {
                        let gap_recovered: Vec<_> = recovered_list
                            .iter()
                            .filter(|(rscaf, rrs, rre)| {
                                *rscaf == ga.scaffold && *rrs >= rs && *rre <= re
                            })
                            .collect();
                        if !gap_recovered.is_empty() {
                            let a_mod = modular_node_fields.contains(&na.nf);
                            let b_mod = modular_node_fields.contains(&nb.nf);
                            if a_mod && !b_mod {
                                let leftmost = gap_recovered
                                    .iter()
                                    .min_by_key(|(_, s, _)| s)
                                    .unwrap();
                                re = leftmost.1 - 1;
                            } else if b_mod && !a_mod {
                                let rightmost = gap_recovered
                                    .iter()
                                    .max_by_key(|(_, _, e)| e)
                                    .unwrap();
                                rs = rightmost.2 + 1;
                            }
                            if rs >= re {
                                continue;
                            }
                        }
                    }
                }
            }

            let scaffold_seq = match genome.get(&ga.scaffold) {
                Some(sequence) => sequence,
                None => continue,
            };

            if re < rs || rs == 0 || re > scaffold_seq.len() {
                continue;
            }

            flog.push(format!(
                "  Gap {}-{}: {}:{}-{} ({})",
                node_a_name, node_b_name, ga.scaffold, rs, re, strand
            ));

            let sr: Vec<u8> = if strand == "+" {
                scaffold_seq[rs - 1..re].to_vec()
            } else {
                bio_revcomp(&scaffold_seq[rs - 1..re])
            };

            // Log the entire NT gap region in FASTA format
            flog.push(format!(
                "  >gap_region {}:{}-{}({}) len={}",
                ga.scaffold, rs, re, strand, sr.len()
            ));
            flog.push(format!("  {}", String::from_utf8_lossy(&sr)));

            if sr.len() < 6 {
                continue;
            }

            if gap_end <= gap_start || gap_end - gap_start < 3 {
                continue;
            }

            let mut lout: Vec<String> = Vec::new();
            let prep = prepare_gap_consensus(
                gap_start,
                gap_end,
                &rc,
                rcount,
                Some(&combined),
                &mut lout,
            );
            let (tcon, rcols, rgaps, _ia) = match prep {
                Some(prepared) => prepared,
                None => {
                    for line in &lout {
                        flog.push(format!("  {}", line));
                    }
                    continue;
                }
            };
            let hits = solve_gap_dp(
                &sr, &tcon, &rcols, &rgaps, aa_idx, ct, &mut lout, score_thr, min_aa,
            );
            for line in &lout {
                flog.push(format!("  {}", line));
            }

            if !hits.is_empty() {
                let mut hits = hits;
                if hits.len() > 1 {
                    let (_pssm, _maxes) = build_blosum_pssm(&tcon, &rcols, aa_idx);
                    hits = refine_boundaries(&hits, &sr);
                }
                thits += hits.len();
                flog.push(format!("  >>> FOUND {} exon(s)", hits.len()));
                let pa: Vec<&str> = na.header.split('|').collect();
                for &(ps, pe, g_s, g_e, fr, sc, ref aa) in &hits {
                    hit_id += 1;
                    let (hss, hse) = if strand == "+" {
                        (rs + g_s, rs + g_e - 1)
                    } else {
                        let rl = re - rs + 1;
                        (rs + (rl - g_e), rs + (rl - g_s) - 1)
                    };
                    let nt = if strand == "+" {
                        scaffold_seq[hss - 1..hse].to_vec()
                    } else {
                        bio_revcomp(&scaffold_seq[hss - 1..hse])
                    };
                    let fv = if strand == "+" {
                        fr as i32 + 1
                    } else {
                        -(fr as i32 + 1)
                    };
                    let dn = format!("DP_{}_{}", gene_base, hit_id);
                    let eh = format!(
                        "{}|{}|{}|NODE_{}|{}|1|{}-{}|{}",
                        pa.first().unwrap_or(&""),
                        pa.get(1).unwrap_or(&""),
                        pa.get(2).unwrap_or(&""),
                        dn,
                        fv,
                        hss,
                        hse,
                        strand
                    );
                    let attrs = format!(
                        "ID={};Name={};Parent={};Note=dp_recovered,frame={},pssm={}-{},score={:.1},aa_len={},between={}+{},cluster={};AA={}",
                        dn, dn, gene_base, fv, ps, pe, sc, aa.len(), node_a_name, node_b_name, ck, aa
                    );
                    gff_lines.push(format!(
                        "{}\tExonDP\texon\t{}\t{}\t{:.1}\t{}\t.\t{}",
                        ga.scaffold, hss, hse, sc, strand, attrs
                    ));
                    recovered.push(RecoveredExon {
                        gene: gene.into(),
                        header: eh,
                        aa_seq: aa.clone(),
                        nt_seq: String::from_utf8_lossy(&nt).into(),
                        region: format!("{}:{}-{}({})", ga.scaffold, hss, hse, strand),
                        score: sc,
                        cluster_key: ck.to_string(),
                        node_a_name: node_a_name.clone(),
                        node_b_name: node_b_name.clone(),
                    });

                    // Track for MXE narrowing
                    base_recovered_gff
                        .entry(ck.to_string())
                        .or_default()
                        .push((ga.scaffold.clone(), hss, hse));
                }
            }
        }
    }

    GeneResult {
        log: flog,
        gaps: tgaps,
        hits: thits,
        gff_lines,
        recovered,
    }
}

fn find_gff_path(folder: &Path, sub_dir: &str) -> Option<String> {
    let candidates = [
        folder.join("outlier").join(sub_dir),
        folder.join("outlier").join("resolve"),
        folder.join(sub_dir),
    ];
    for dir in &candidates {
        if dir.exists() {
            if let Ok(entries) = fs::read_dir(dir) {
                let mut gffs: Vec<String> = entries
                    .flatten()
                    .filter_map(|e| {
                        let name = e.file_name().to_string_lossy().to_string();
                        if name.ends_with("_coords.gff") {
                            Some(e.path().to_string_lossy().to_string())
                        } else {
                            None
                        }
                    })
                    .collect();
                gffs.sort();
                if let Some(path) = gffs.first() {
                    return Some(path.clone());
                }
            }
        }
    }
    None
}

fn find_input_folder(folder: &Path, sub_dir: &str) -> Option<std::path::PathBuf> {
    let p1 = folder.join("outlier").join(sub_dir);
    if p1.exists() {
        return Some(p1);
    }
    let p2 = folder.join(sub_dir);
    if p2.exists() {
        return Some(p2);
    }
    None
}

#[pyfunction]
#[pyo3(signature = (folder, sub_dir))]
pub fn exon_dp(py: Python<'_>, folder: String, sub_dir: String) -> PyResult<PyObject> {
    let folder_path = Path::new(&folder);

    let input_folder = find_input_folder(folder_path, &sub_dir).ok_or_else(|| {
        pyo3::exceptions::PyFileNotFoundError::new_err(format!(
            "Input folder not found for sub_dir={}",
            sub_dir
        ))
    })?;

    let aa_dir = input_folder.join("aa");

    let gff_path = find_gff_path(folder_path, &sub_dir);

    let mut gff_nodes = match &gff_path {
        Some(path) => read_gff(path),
        None => HashMap::new(),
    };

    let cluster_csv = input_folder.join("resolve_clusters.csv");
    let gene_clusters = if cluster_csv.exists() {
        read_clusters(&cluster_csv.to_string_lossy())
    } else {
        HashMap::new()
    };

    // Determine which scaffolds are actually needed before hitting RocksDB.
    // If no cluster nodes map to any GFF entry we have nothing to process.
    let needed_scaffolds = collect_needed_scaffolds(&gene_clusters, &gff_nodes);

    if needed_scaffolds.is_empty() {
        let output = PyDict::new(py);
        output.set_item("results", PyList::empty(py))?;
        output.set_item("gff_nodes", PyDict::new(py))?;
        output.set_item("genes", PyList::empty(py))?;
        output.set_item("genome_count", 0usize)?;
        return Ok(output.into());
    }

    let rocksdb_path = folder_path.join("rocksdb").join("sequences").join("nt");
    let genome = load_genome_from_rocksdb(&rocksdb_path.to_string_lossy(), &needed_scaffolds);

    if genome.is_empty() {
        return Err(pyo3::exceptions::PyRuntimeError::new_err(format!(
            "No parent sequences found in RocksDB at {:?}",
            rocksdb_path
        )));
    }

    let aa_idx = build_aa_idx_table();
    let ct = build_codon_table(&aa_idx);

    let mut gene_files: Vec<String> = Vec::new();
    if let Ok(entries) = fs::read_dir(&aa_dir) {
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().to_string();
            if name.ends_with(".fa") || name.contains(".fa.") {
                gene_files.push(name);
            }
        }
    }
    gene_files.sort();

    let score_thr = SCORE_FLOOR_FRAC_DEFAULT;
    let min_aa = MIN_EXON_AA_DEFAULT;

    let mut results_list: Vec<GeneResult> = Vec::new();
    let mut output_genes: Vec<String> = Vec::new();
    let mut gene_entries_cache: Vec<Vec<(String, String)>> = Vec::new();

    for gene_file in &gene_files {
        let gene_key = gene_file.trim_end_matches(".gz");
        if !gene_clusters.contains_key(gene_key) {
            continue;
        }

        let fpath = aa_dir.join(gene_file);
        let entries = read_aa_fasta(&fpath.to_string_lossy());

        if entries.is_empty() {
            continue;
        }

        output_genes.push(gene_file.clone());

        results_list.push(process_gene(
            gene_key,
            &entries,
            &gene_clusters,
            &gff_nodes,
            &genome,
            &aa_idx,
            &ct,
            score_thr,
            min_aa,
        ));
        gene_entries_cache.push(entries);
    }

    // -----------------------------------------------------------------------
    // Motif scan: update GFF with recovered exons, then scan flanks
    // -----------------------------------------------------------------------
    // Register recovered DP nodes into gff_nodes so motif scan can see them
    for result in &results_list {
        for rec in &result.recovered {
            let dp_name = parse_node_field(&rec.header).replace("NODE_", "");
            // Parse region string: "scaffold:start-end(strand)"
            if let Some(colon) = rec.region.find(':') {
                if let Some(paren) = rec.region.find('(') {
                    let scaffold = &rec.region[..colon];
                    let coords = &rec.region[colon + 1..paren];
                    let strand = &rec.region[paren + 1..paren + 2];
                    let parts: Vec<&str> = coords.split('-').collect();
                    if parts.len() == 2 {
                        if let (Ok(s), Ok(e)) = (parts[0].parse(), parts[1].parse()) {
                            gff_nodes.insert(
                                dp_name,
                                GffEntry {
                                    scaffold: scaffold.to_string(),
                                    start: s,
                                    end: e,
                                    strand: strand.to_string(),
                                },
                            );
                        }
                    }
                }
            }
        }
    }

    let scaffold_intervals = build_scaffold_intervals(&gff_nodes);

    let mut motif_log_all: Vec<String> = Vec::new();
    let mut motif_results: Vec<(String, Vec<MotifHit>)> = Vec::new();
    let mut motif_gff_lines: Vec<String> = Vec::new();

    for (gi, gene_file) in output_genes.iter().enumerate() {
        let gene_key = gene_file.trim_end_matches(".gz");
        let gene_base = gene_key.split('.').next().unwrap_or(gene_key);
        let entries = &gene_entries_cache[gi];

        let (log, hits) = motif_scan_gene(
            gene_key,
            entries,
            &gene_clusters,
            &gff_nodes,
            &genome,
            &scaffold_intervals,
        );

        if !log.is_empty() {
            motif_log_all.extend(log);
        }
        if !hits.is_empty() {
            // Generate GFF lines and register motif nodes in gff_nodes
            for hit in &hits {
                // Parse region string: "scaffold:start-end(strand)"
                if let Some(colon) = hit.region.find(':') {
                    if let Some(paren) = hit.region.find('(') {
                        let scaffold = &hit.region[..colon];
                        let coords = &hit.region[colon + 1..paren];
                        let strand = &hit.region[paren + 1..paren + 2];
                        let parts: Vec<&str> = coords.split('-').collect();
                        if parts.len() == 2 {
                            if let (Ok(s), Ok(e)) = (parts[0].parse::<usize>(), parts[1].parse::<usize>()) {
                                let motif_node_field = hit.header.split('|').nth(3).unwrap_or("").replace("NODE_", "");
                                let aa_clean: String = hit.aa_seq.chars().filter(|&c| c != '-').collect();
                                let attrs = format!(
                                    "ID={};Name={};Parent={};Note=motif_scan,score={},aa_len={}",
                                    motif_node_field, motif_node_field, gene_base, hit.score, aa_clean.len()
                                );
                                motif_gff_lines.push(format!(
                                    "{}\tMotifScan\texon\t{}\t{}\t{}\t{}\t.\t{}",
                                    scaffold, s, e, hit.score, strand, attrs
                                ));

                                // Register in gff_nodes so Python can sort them
                                gff_nodes.insert(
                                    motif_node_field,
                                    GffEntry {
                                        scaffold: scaffold.to_string(),
                                        start: s,
                                        end: e,
                                        strand: strand.to_string(),
                                    },
                                );
                            }
                        }
                    }
                }
            }
            motif_results.push((gene_key.to_string(), hits));
        }
    }

    // -----------------------------------------------------------------------
    // Build Python output
    // -----------------------------------------------------------------------
    let output = PyDict::new(py);

    let py_results = PyList::empty(py);
    for result in &results_list {
        let rd = PyDict::new(py);
        rd.set_item("log", &result.log)?;
        rd.set_item("gaps", result.gaps)?;
        rd.set_item("hits", result.hits)?;
        rd.set_item("gff_lines", &result.gff_lines)?;

        let py_recovered = PyList::empty(py);
        for rec in &result.recovered {
            let rec_dict = PyDict::new(py);
            rec_dict.set_item("gene", &rec.gene)?;
            rec_dict.set_item("header", &rec.header)?;
            rec_dict.set_item("aa_seq", &rec.aa_seq)?;
            rec_dict.set_item("nt_seq", &rec.nt_seq)?;
            rec_dict.set_item("region", &rec.region)?;
            rec_dict.set_item("score", rec.score)?;
            rec_dict.set_item("cluster_key", &rec.cluster_key)?;
            rec_dict.set_item("node_a_name", &rec.node_a_name)?;
            rec_dict.set_item("node_b_name", &rec.node_b_name)?;
            py_recovered.append(rec_dict)?;
        }
        rd.set_item("recovered", py_recovered)?;
        py_results.append(rd)?;
    }
    output.set_item("results", py_results)?;

    // Motif scan results
    let py_motif = PyList::empty(py);
    for (gene, hits) in &motif_results {
        let gd = PyDict::new(py);
        gd.set_item("gene", gene)?;
        let py_hits = PyList::empty(py);
        for hit in hits {
            let hd = PyDict::new(py);
            hd.set_item("header", &hit.header)?;
            hd.set_item("aa_seq", &hit.aa_seq)?;
            hd.set_item("nt_seq", &hit.nt_seq)?;
            hd.set_item("region", &hit.region)?;
            hd.set_item("score", hit.score)?;
            hd.set_item("node_name", &hit.node_name)?;
            hd.set_item("is_leading", hit.is_leading)?;
            py_hits.append(hd)?;
        }
        gd.set_item("hits", py_hits)?;
        py_motif.append(gd)?;
    }
    output.set_item("motif_results", py_motif)?;
    output.set_item("motif_log", motif_log_all)?;
    output.set_item("motif_gff_lines", motif_gff_lines)?;

    let py_gff = PyDict::new(py);
    for (name, entry) in &gff_nodes {
        let ed = PyDict::new(py);
        ed.set_item("scaffold", &entry.scaffold)?;
        ed.set_item("start", entry.start)?;
        ed.set_item("end", entry.end)?;
        ed.set_item("strand", &entry.strand)?;
        py_gff.set_item(name, ed)?;
    }
    output.set_item("gff_nodes", py_gff)?;
    output.set_item("genes", &output_genes)?;
    output.set_item("genome_count", genome.len())?;

    Ok(output.into())
}
