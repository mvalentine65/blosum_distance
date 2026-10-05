//! Per-exon scores of a chain. Each exon is translated in its best frame and
//! the chain's peptide is aligned once (Forward).
//! bits: BATH's p7_splice_ScoreExons score, the Forward gain across the exon's
//! residues (ln C at its last residue minus ln C before its first), with the
//! length model reset to the exon and its null score removed.
//! rev: the chain's Forward bits minus the same with the exon's residues
//! reversed in place; below 0 the exon fits no better than its reversal.
//! nodes: the model node each residue matches (0 for none), for MAP placement
//! in the MSA: from the chain peptide's optimal-accuracy path, or from the
//! exon's own path against its node range ±5 when that matches more residues
//! (a weak end exon can fall outside the chain's single domain).

use super::chain::Chain;
use super::hmm::Hmm;
use super::orf::{translate_frame, Aligner};
use super::sites::revcomp;

#[derive(Clone, Debug)]
pub struct ExonScore {
    pub frame: i64,
    pub aa: usize,
    pub bits: f64,
    pub rev: f64,
    pub nodes: Vec<u32>,
    /// genomic position (first base) of each in-frame stop codon in the exon
    pub stops: Vec<i64>,
}

impl Default for ExonScore {
    fn default() -> Self { ExonScore { frame: -1, aa: 0, bits: f64::NAN, rev: f64::NAN, nodes: Vec::new(), stops: Vec::new() } }
}

/// spans: (start, end) of each exon, genomic, in chain order. rev is computed
/// for exons where `want_rev` holds; residue nodes when `nodes` is set; bits
/// for all when `bits` is set. Frames are always found.
#[allow(clippy::too_many_arguments)]
pub fn score_chain(al: &mut Aligner, hmm_id: usize, hmm: &Hmm, c: &Chain, spans: &[(i64, i64)], sc: &[u8],
                   nodes: bool, bits: bool, want_rev: impl Fn(usize) -> bool) -> Vec<ExonScore> {
    let m = hmm.m as i64;
    let mut out = vec![ExonScore::default(); c.ex.len()];
    let mut pep: Vec<u8> = Vec::new();
    let mut res: Vec<(usize, usize)> = Vec::with_capacity(c.ex.len());
    let (mut a, mut b) = (i64::MAX, 0i64);
    for (j, e) in c.ex.iter().enumerate() {
        let (s, t) = spans[j];
        let (k1, k2) = (e.k1.clamp(1, m), e.k2.clamp(1, m));
        let lo = pep.len();
        if s >= 1 && t <= sc.len() as i64 && t - s + 1 >= 3 && k1 <= k2 {
            let fwd = &sc[(s - 1) as usize..t as usize];
            let nt = if c.strand == b'-' { revcomp(fwd) } else { fwd.to_vec() };
            // best frame by Forward score alone (the full decoding pipeline isn't needed here)
            al.set(hmm_id, hmm, k1 as usize, k2 as usize);
            let (mut f, mut best) = (0i64, f32::NEG_INFINITY);
            let mut frames: Vec<Vec<u8>> = Vec::with_capacity(3);
            for fr in 0..3 {
                let t: Vec<u8> = translate_frame(&nt, fr).into_iter().map(|x| if x == b'*' { b'X' } else { x }).collect();
                if !t.is_empty() {
                    let b = al.fwd_bits(&t);
                    if b > best { best = b; f = fr as i64; }
                }
                frames.push(t);
            }
            // in-frame stops, as genomic positions of the codon's first base;
            // the gene's own stop (the last codon of the chain's last exon) is not one
            let last = j + 1 == c.ex.len();
            for (i, &x) in translate_frame(&nt, f as usize).iter().enumerate() {
                if x == b'*' {
                    let o = f + 3 * i as i64; // offset in coding orientation
                    if last && o + 3 + 2 >= nt.len() as i64 { continue; }
                    out[j].stops.push(if c.strand == b'-' { t - o } else { s + o });
                }
            }
            pep.extend(frames.swap_remove(f as usize));
            out[j].frame = f;
            a = a.min(k1);
            b = b.max(k2);
        }
        res.push((lo, pep.len()));
        out[j].aa = pep.len() - lo;
    }
    let n = pep.len();
    if n == 0 { return out; }
    al.set(hmm_id, hmm, a as usize, b as usize);
    let tested: Vec<usize> = (0..res.len()).filter(|&j| res[j].1 > res[j].0 && want_rev(j)).collect();
    if !tested.is_empty() {
        // the chain with an exon reversed resumes from the row saved before that exon
        let marks: Vec<usize> = tested.iter().map(|&j| res[j].0).collect();
        let (full, rows) = al.fwd_bits_chain(&pep, &marks);
        for (t, &j) in tested.iter().enumerate() {
            let (lo, hi) = res[j];
            let mut r = pep.clone();
            r[lo..hi].reverse();
            out[j].rev = full as f64 - al.fwd_bits_chain_from(&r, rows.get(t), lo) as f64;
        }
    }
    if !nodes && !bits { return out; }
    for (j, &(lo, hi)) in res.iter().enumerate() { out[j].nodes = vec![0; hi - lo]; }
    for (node, i) in al.chain_trace(&pep) {
        if i == 0 || i > n { continue; }
        if let Some(j) = res.iter().position(|&(lo, hi)| i > lo && i <= hi) { out[j].nodes[i - 1 - res[j].0] = node as u32; }
    }
    for (j, e) in c.ex.iter().enumerate() {
        let (lo, hi) = res[j];
        let got = out[j].nodes.iter().filter(|&&k| k > 0).count();
        if hi == lo || got == hi - lo { continue; }
        let (k1, k2) = (e.k1.clamp(1, m), e.k2.clamp(1, m));
        al.set(hmm_id, hmm, (k1 - 5).max(1) as usize, (k2 + 5).min(m) as usize);
        let mut own = vec![0u32; hi - lo];
        for (node, i) in al.one_trace(&pep[lo..hi]) {
            if i >= 1 && i <= hi - lo { own[i - 1] = node as u32; }
        }
        if own.iter().filter(|&&k| k > 0).count() > got { out[j].nodes = own; }
    }
    al.set(hmm_id, hmm, a as usize, b as usize);
    if !bits { return out; }
    let lnc = al.fwd_lnc(&pep);
    let ln2 = std::f64::consts::LN_2;
    for (j, &(lo, hi)) in res.iter().enumerate() {
        let len = hi - lo;
        if len == 0 { continue; }
        let start = if lo == 0 { 0.0 } else { lnc[lo] };
        let mut s = lnc[hi] - start;
        s -= (2.0 / (n as f64 + 2.0)).ln();
        s += 2.0 * (2.0 / (len as f64 + 2.0)).ln();
        let p1 = len as f64 / (len as f64 + 1.0);
        let nullsc = len as f64 * p1.ln() + (1.0 - p1).ln();
        out[j].bits = (s - nullsc) / ln2;
    }
    out
}
