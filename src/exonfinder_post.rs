//! Rust port of `sapphyre/outlier/exonfinder/post_process.py` + the parts of
//! `recovery.py` it depends on.  The single PyO3 entry point
//! [`exonfinder_process_gene`] runs the entire per-gene worker without
//! returning to Python in between -- hmm_align, column cull, low-complexity
//! filter, PSSM scoring, dedup, cross-cluster genomic dedup, IMX late kick,
//! native supersession, tag assignment, and the final per-cluster overlap +
//! module split all happen inside this call.
//!
//! Logging and rejection strings are emitted byte-for-byte identical to the
//! Python implementation so existing log parsers and audit pipelines keep
//! working.

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::time::Instant;

use crate::aligner::hmm_align;
use crate::column_cull::cull_columns;

// ===========================================================================
// Constants -- mirror sapphyre/outlier/exonfinder/core.py
// ===========================================================================

const DEDUP_MAX_OVERLAP: usize = 8;
const MIN_MODULE_SIZE_RATIO: f64 = 0.75;
const MIN_MODULE_COL_OVERLAP: f64 = 0.80;
const MIN_MODULE_COL_OVERLAP_LARGER: f64 = 0.75;
const MIN_INTRON_BP: usize = 30;

const NATIVE_OVERLAP_FRAC: f64 = 0.80;

const MIN_GAP_RESIDUES: usize = 15;
const MIN_COVERED_PCT: f64 = 20.0;

const MIN_AA_AFTER_SPLIT: usize = 15;

const MODULE_OVERLAP_FRAC: f64 = 0.5;
const MODULE_SCORE_RATIO: f64 = 0.85;

const MODULE_MIN_SIZE_RATIO_SPLIT: f64 = 0.75;
const MODULE_MIN_INTRON_BP_SPLIT: usize = 15;
const MODULE_MIN_OVERLAP_SPLIT: f64 = 0.80;

const CROSS_CLUSTER_OVERLAP_FRAC: f64 = 0.5;

// PairKind tags returned by classify_pair.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum PairKind {
    Independent,
    MonotonicityBreak,
    GenomicOverlap,
    Module,
    SlotDuplicate,
}

// ===========================================================================
// Input structs extracted from Python dicts
// ===========================================================================

#[derive(FromPyObject, Clone)]
#[pyo3(from_item_all)]
pub struct FlankInput {
    pub header: String,
    pub tag: String,
    pub gene_key: String,
    pub aa_seq: String,
    pub nt_seq: String,
    pub scaffold: String,
    pub strand: String,
    pub hit_start: usize,
    pub hit_end: usize,
    pub frame: i32,
    pub gap_start: usize,
    pub gap_end: usize,
    pub node_name: String,
    pub is_leading: bool,
    pub cluster_key: String,
}

#[derive(FromPyObject, Clone)]
#[pyo3(from_item_all)]
pub struct GapInput {
    pub header: String,
    pub tag: String,
    pub gene_key: String,
    pub aa_seq: String,
    pub nt_seq: String,
    pub scaffold: String,
    pub strand: String,
    pub hit_start: usize,
    pub hit_end: usize,
    pub frame: i32,
    pub gap_start: usize,
    pub gap_end: usize,
    pub node_a_name: String,
    pub node_b_name: String,
    pub imx_slot_node: String,
    pub cluster_key: String,
}

#[derive(FromPyObject, Clone)]
#[pyo3(from_item_all)]
pub struct GffEntryInput {
    pub scaffold: String,
    pub start: usize,
    pub end: usize,
    pub strand: String,
}

// ===========================================================================
// Internal Recovery types
// ===========================================================================

#[derive(Clone, Debug)]
struct Locus {
    scaffold: String,
    strand: String,
    bp_start: usize,
    bp_end: usize,
    msa_cols: BTreeSet<usize>,
    frame: i32,
}

impl Locus {
    fn trim_aa(&self, left_aa: usize, right_aa: usize) -> Locus {
        if left_aa == 0 && right_aa == 0 {
            return self.clone();
        }
        let mut out = self.clone();
        if self.strand == "-" {
            out.bp_start = self.bp_start.saturating_add(right_aa * 3);
            out.bp_end = self.bp_end.saturating_sub(left_aa * 3);
        } else {
            out.bp_start = self.bp_start.saturating_add(left_aa * 3);
            out.bp_end = self.bp_end.saturating_sub(right_aa * 3);
        }
        out
    }
}

#[derive(Clone, Debug)]
struct ScoreInfo {
    full_blosum: f64,
    covered_blosum: f64,
    gap_residues: usize,
    win_start: usize,
    win_end: usize,
}

#[derive(Clone, Debug)]
struct FlankAnchor {
    adj_node: String,
    is_leading: bool,
}

#[derive(Clone, Debug)]
struct GapAnchor {
    node_a: String,
    node_b: String,
}

#[derive(Clone, Debug)]
struct ImxAnchor {
    slot_node: String,
    base_cluster_key: String,
    left_flank: String,
    right_flank: String,
}

#[derive(Clone, Debug)]
struct ModuleRel {
    partner_tag: String,
    module_type: String,
    slot_cols: BTreeSet<usize>,
}

#[derive(Clone, Debug, Copy, PartialEq, Eq)]
enum RecoveryKind {
    Flank,
    Gap,
}

impl RecoveryKind {
    fn as_str(self) -> &'static str {
        match self {
            RecoveryKind::Flank => "FLANK",
            RecoveryKind::Gap => "GAP",
        }
    }
}

#[derive(Clone, Debug)]
struct Recovery {
    tag: String,
    header: String,
    aa_seq: String,
    nt_seq: String,
    locus: Locus,
    score: ScoreInfo,
    kind: RecoveryKind,
    cluster_key: String,
    gene_key: String,
    flank_anchor: Option<FlankAnchor>,
    gap_anchor: Option<GapAnchor>,
    imx_anchor: Option<ImxAnchor>,
    module: Option<ModuleRel>,
    final_tag: Option<String>,
}

impl Recovery {
    fn current_tag(&self) -> &str {
        self.final_tag.as_deref().unwrap_or(&self.tag)
    }
}

// ===========================================================================
// Helpers
// ===========================================================================

#[inline]
fn is_gap(c: u8) -> bool {
    c == b'-' || c == b'.'
}

fn count_residues(seq: &str) -> usize {
    seq.bytes().filter(|&c| !is_gap(c)).count()
}

fn data_cols_set(seq: &str) -> BTreeSet<usize> {
    seq.bytes()
        .enumerate()
        .filter_map(|(i, c)| if !is_gap(c) { Some(i) } else { None })
        .collect()
}

fn seq_residue_bounds(seq: &str) -> (i64, i64) {
    let mut first: i64 = -1;
    let mut last: i64 = -1;
    for (i, c) in seq.bytes().enumerate() {
        if !is_gap(c) {
            if first == -1 {
                first = i as i64;
            }
            last = i as i64;
        }
    }
    (first, last)
}

fn interval_overlap_bp(s_a: usize, e_a: usize, s_b: usize, e_b: usize) -> i64 {
    let lo = s_a.max(s_b) as i64;
    let hi = e_a.min(e_b) as i64;
    (hi - lo + 1).max(0)
}

fn tag_from_header(hdr: &str) -> &str {
    let mut it = hdr.split('|');
    for _ in 0..3 {
        if it.next().is_none() {
            return hdr;
        }
    }
    match it.next() {
        Some(t) => t,
        None => hdr,
    }
}

// AA index for BLOSUM62 / aa_index_arr lookups (20 = unknown sentinel).
#[inline]
fn aa_index(c: u8) -> usize {
    match c {
        b'A' => 0, b'R' => 1, b'N' => 2, b'D' => 3, b'C' => 4,
        b'Q' => 5, b'E' => 6, b'G' => 7, b'H' => 8, b'I' => 9,
        b'L' => 10, b'K' => 11, b'M' => 12, b'F' => 13, b'P' => 14,
        b'S' => 15, b'T' => 16, b'W' => 17, b'Y' => 18, b'V' => 19,
        b'a' => 0, b'r' => 1, b'n' => 2, b'd' => 3, b'c' => 4,
        b'q' => 5, b'e' => 6, b'g' => 7, b'h' => 8, b'i' => 9,
        b'l' => 10, b'k' => 11, b'm' => 12, b'f' => 13, b'p' => 14,
        b's' => 15, b't' => 16, b'w' => 17, b'y' => 18, b'v' => 19,
        _ => 20,
    }
}

#[rustfmt::skip]
const BLOSUM62: [[i8; 20]; 20] = [
  //  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    [ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0], // A
    [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3], // R
    [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3], // N
    [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3], // D
    [ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1], // C
    [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2], // Q
    [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2], // E
    [ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3], // G
    [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3], // H
    [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3], // I
    [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1], // L
    [-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2], // K
    [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1], // M
    [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1], // F
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2], // P
    [ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2], // S
    [ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0], // T
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3], // W
    [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1], // Y
    [ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4], // V
];

fn blosum62_score(a: u8, b: u8) -> i32 {
    let ai = aa_index(a);
    let bi = aa_index(b);
    if ai >= 20 || bi >= 20 {
        return -4;
    }
    BLOSUM62[ai][bi] as i32
}

/// Fraction of pairwise BLOSUM-positive aligned residues between two
/// equal-length aligned sequences; mirrors `core._blosum_identity`.
fn blosum_identity(seq_a: &str, seq_b: &str) -> f64 {
    if seq_a.is_empty() || seq_b.is_empty() {
        return 0.0;
    }
    let a = seq_a.as_bytes();
    let b = seq_b.as_bytes();
    let n = a.len().min(b.len());
    let mut scored = 0usize;
    let mut matched = 0usize;
    for i in 0..n {
        let ca = a[i];
        let cb = b[i];
        if is_gap(ca) || is_gap(cb) {
            continue;
        }
        let ai = aa_index(ca);
        let bi = aa_index(cb);
        if ai >= 20 || bi >= 20 {
            continue;
        }
        scored += 1;
        if BLOSUM62[ai][bi] > 0 {
            matched += 1;
        }
    }
    if scored == 0 {
        0.0
    } else {
        matched as f64 / scored as f64
    }
}

fn fmt_region(loc: &Locus) -> String {
    format!(
        "{}:{}-{}({})",
        loc.scaffold, loc.bp_start, loc.bp_end, loc.strand
    )
}

// ===========================================================================
// PSSM scoring (window cache + per-candidate scoring) -- mirrors
// `core._score_window` / `_score_candidate` but the column scores stay in a
// dense [f64; 20] array (no HashMap allocation per column) and the candidate
// loop runs in tight bytes.
// ===========================================================================

/// Per-column PSSM entry: ([score per AA index], col_max).
type PssmCol = ([f64; 20], f64);

/// Window-keyed cache: (win_start, win_end) -> (data_cols, per-column pssm, full_max).
struct WindowCache {
    cache: HashMap<(usize, usize), (Vec<usize>, Vec<Option<PssmCol>>, f64)>,
    aln_len: usize,
}

impl WindowCache {
    fn new(aln_len: usize) -> Self {
        Self {
            cache: HashMap::new(),
            aln_len,
        }
    }

    fn get_or_build<'a>(
        &'a mut self,
        ref_seqs: &[&[u8]],
        win_start: usize,
        win_end: usize,
    ) -> &'a (Vec<usize>, Vec<Option<PssmCol>>, f64) {
        let key = (win_start, win_end);
        if !self.cache.contains_key(&key) {
            let built = build_pssm_for_window(ref_seqs, win_start, win_end, self.aln_len);
            self.cache.insert(key, built);
        }
        self.cache.get(&key).unwrap()
    }
}

/// Combined ref-data-cols + PSSM build, returning a column-indexed
/// `Vec<Option<PssmCol>>` of length `aln_len` so candidate scoring is a
/// branch-free linear scan.
fn build_pssm_for_window(
    ref_seqs: &[&[u8]],
    win_start: usize,
    win_end: usize,
    aln_len: usize,
) -> (Vec<usize>, Vec<Option<PssmCol>>, f64) {
    let mut pssm: Vec<Option<PssmCol>> = vec![None; aln_len];
    let mut data_cols: Vec<usize> = Vec::new();
    let mut full_max = 0.0f64;

    let n_refs = ref_seqs.len();
    if n_refs == 0 || win_end <= win_start {
        return (data_cols, pssm, full_max);
    }
    let n_refs_f = n_refs as f64;
    let min_ref_occ = 0.30f64;

    for col in win_start..win_end.min(aln_len) {
        let mut occ_count = 0usize;
        let mut counts = [0i32; 20];
        let mut total = 0i32;
        for r in ref_seqs {
            if col < r.len() {
                let c = r[col];
                if !is_gap(c) {
                    occ_count += 1;
                    let ai = aa_index(c);
                    if ai < 20 {
                        counts[ai] += 1;
                        total += 1;
                    }
                }
            }
        }
        if (occ_count as f64) / n_refs_f < min_ref_occ {
            continue;
        }
        data_cols.push(col);
        if total == 0 {
            continue;
        }
        let total_f = total as f64;
        let mut scores = [0f64; 20];
        let mut col_max = f64::NEG_INFINITY;
        for q in 0..20usize {
            let qb = match q {
                0 => b'A', 1 => b'R', 2 => b'N', 3 => b'D', 4 => b'C',
                5 => b'Q', 6 => b'E', 7 => b'G', 8 => b'H', 9 => b'I',
                10 => b'L', 11 => b'K', 12 => b'M', 13 => b'F', 14 => b'P',
                15 => b'S', 16 => b'T', 17 => b'W', 18 => b'Y', 19 => b'V',
                _ => unreachable!(),
            };
            let mut s = 0f64;
            for (i, &cnt) in counts.iter().enumerate() {
                if cnt > 0 {
                    let rb = match i {
                        0 => b'A', 1 => b'R', 2 => b'N', 3 => b'D', 4 => b'C',
                        5 => b'Q', 6 => b'E', 7 => b'G', 8 => b'H', 9 => b'I',
                        10 => b'L', 11 => b'K', 12 => b'M', 13 => b'F', 14 => b'P',
                        15 => b'S', 16 => b'T', 17 => b'W', 18 => b'Y', 19 => b'V',
                        _ => unreachable!(),
                    };
                    s += blosum62_score(qb, rb) as f64 * cnt as f64;
                }
            }
            let normalized = s / total_f;
            scores[q] = normalized;
            if normalized > col_max {
                col_max = normalized;
            }
        }
        if col_max > 0.0 {
            pssm[col] = Some((scores, col_max));
            full_max += col_max;
        }
    }

    (data_cols, pssm, full_max)
}

/// Window scoring result -- mirrors `_WindowScore`.
#[derive(Clone, Debug)]
struct WindowScore {
    gap_residues: usize,
    blosum_pct: f64,
    covered_pct: f64,
    scored_cols: usize,
    pssm_total_cols: usize,
    win_start: usize,
    win_end: usize,
}

/// Score one candidate sequence against a cached window PSSM.
fn score_window(
    cand_seq: &[u8],
    ref_seqs: &[&[u8]],
    win_start: usize,
    win_end: usize,
    cache: &mut WindowCache,
) -> WindowScore {
    let (data_cols, pssm, full_max) = cache.get_or_build(ref_seqs, win_start, win_end);

    let mut gap_residues = 0usize;
    for &col in data_cols {
        if col < cand_seq.len() && !is_gap(cand_seq[col]) {
            gap_residues += 1;
        }
    }

    let mut cand_total = 0f64;
    let mut max_total = 0f64;
    let mut scored_cols = 0usize;
    let mut pssm_total_cols = 0usize;
    for (col, entry) in pssm.iter().enumerate() {
        let Some((scores, col_max)) = entry else { continue };
        pssm_total_cols += 1;
        if col >= cand_seq.len() {
            continue;
        }
        let c = cand_seq[col];
        if is_gap(c) {
            continue;
        }
        let ai = aa_index(c);
        let v = if ai < 20 { scores[ai] } else { -4.0 };
        cand_total += v;
        max_total += *col_max;
        scored_cols += 1;
    }
    let blosum_pct = if *full_max > 0.0 {
        (cand_total / *full_max) * 100.0
    } else {
        0.0
    };
    let covered_pct = if max_total > 0.0 {
        (cand_total / max_total) * 100.0
    } else {
        0.0
    };

    WindowScore {
        gap_residues,
        blosum_pct,
        covered_pct,
        scored_cols,
        pssm_total_cols,
        win_start,
        win_end,
    }
}

// ===========================================================================
// classify_pair + dedup_group -- mirrors recovery.classify_pair / _dedup_group
// ===========================================================================

#[derive(Clone, Debug)]
struct PairRelation {
    kind: PairKind,
    genomic_dist: usize,
}

fn classify_pair(
    cand: &Recovery,
    acc: &Recovery,
    cand_is_imx_scan: bool,
) -> PairRelation {
    let cand_cols = &cand.locus.msa_cols;
    let acc_cols = &acc.locus.msa_cols;
    let col_overlap = cand_cols.intersection(acc_cols).count();

    let cand_msa_min = cand_cols.iter().next().copied().unwrap_or(0);
    let cand_msa_max = cand_cols.iter().next_back().copied().unwrap_or(0);
    let acc_msa_min = acc_cols.iter().next().copied().unwrap_or(0);
    let acc_msa_max = acc_cols.iter().next_back().copied().unwrap_or(0);

    let cand_span = (cand_msa_max + 1).saturating_sub(cand_msa_min).max(1);
    let acc_span = (acc_msa_max + 1).saturating_sub(acc_msa_min).max(1);
    let overlap_min = cand_msa_min.max(acc_msa_min);
    let overlap_max = cand_msa_max.min(acc_msa_max);
    let msa_overlap_len = if overlap_max >= overlap_min {
        (overlap_max - overlap_min) + 1
    } else {
        0
    };
    let shorter_span = cand_span.min(acc_span);
    let msa_overlap_frac = if shorter_span > 0 {
        msa_overlap_len as f64 / shorter_span as f64
    } else {
        0.0
    };

    let same_locus =
        cand.locus.scaffold == acc.locus.scaffold && cand.locus.strand == acc.locus.strand;
    let genomic_overlap = same_locus
        && interval_overlap_bp(
            cand.locus.bp_start, cand.locus.bp_end,
            acc.locus.bp_start, acc.locus.bp_end,
        ) > 0;
    let genomic_dist = if same_locus {
        (cand.locus.bp_start as i64 - acc.locus.bp_start as i64).unsigned_abs() as usize
    } else {
        0
    };

    let smaller = cand.score.gap_residues.min(acc.score.gap_residues);
    let larger = cand.score.gap_residues.max(acc.score.gap_residues);
    let size_ratio = if larger > 0 { smaller as f64 / larger as f64 } else { 0.0 };
    let col_overlap_frac_smaller = if smaller > 0 {
        col_overlap as f64 / smaller as f64
    } else {
        0.0
    };
    let col_overlap_frac_larger = if larger > 0 {
        col_overlap as f64 / larger as f64
    } else {
        0.0
    };
    let score_ratio_ok = acc.score.covered_blosum <= 0.0
        || cand.score.covered_blosum >= acc.score.covered_blosum * MODULE_SCORE_RATIO;

    let msa_module_skip =
        shorter_span > 0 && msa_overlap_frac >= MODULE_OVERLAP_FRAC;

    if col_overlap > DEDUP_MAX_OVERLAP {
        if genomic_overlap {
            return PairRelation { kind: PairKind::SlotDuplicate, genomic_dist };
        }
        let module_geometry = same_locus
            && larger > 0
            && size_ratio >= MIN_MODULE_SIZE_RATIO
            && col_overlap_frac_smaller >= MIN_MODULE_COL_OVERLAP
            && col_overlap_frac_larger >= MIN_MODULE_COL_OVERLAP_LARGER
            && genomic_dist >= MIN_INTRON_BP
            && score_ratio_ok;
        if module_geometry {
            return PairRelation { kind: PairKind::Module, genomic_dist };
        }
        return PairRelation { kind: PairKind::SlotDuplicate, genomic_dist };
    }

    if !same_locus {
        return PairRelation { kind: PairKind::Independent, genomic_dist };
    }

    if !cand_is_imx_scan && !msa_module_skip {
        let mono_break = if cand.locus.strand == "+" {
            (cand.locus.bp_start > acc.locus.bp_start && cand_msa_min < acc_msa_min)
                || (cand.locus.bp_start < acc.locus.bp_start && cand_msa_min > acc_msa_min)
        } else {
            (cand.locus.bp_start > acc.locus.bp_start && cand_msa_min > acc_msa_min)
                || (cand.locus.bp_start < acc.locus.bp_start && cand_msa_min < acc_msa_min)
        };
        if mono_break {
            return PairRelation { kind: PairKind::MonotonicityBreak, genomic_dist };
        }
    }

    if genomic_overlap {
        return PairRelation { kind: PairKind::GenomicOverlap, genomic_dist };
    }

    PairRelation { kind: PairKind::Independent, genomic_dist }
}

#[derive(Default, Clone, Debug)]
struct DedupCounts {
    n_in: usize,
    accepted: usize,
    module: usize,
    rej_mono: usize,
    rej_genomic: usize,
    rej_overlap: usize,
}

/// Greedy per-group dedup; mirrors `_dedup_group`.  Mutates incoming
/// recoveries (sets `module` field on module partners).  Returns the
/// indices of accepted entries plus counts.
fn dedup_group(
    sorted_cands: &mut [(String, String, Recovery)],
    scope: &str,
    group_label: &str,
    module_type_for: impl Fn(&Recovery) -> &'static str,
    rejection_log: &mut Vec<String>,
    debug_level: i32,
) -> (Vec<usize>, DedupCounts) {
    let mut accepted_idx: Vec<usize> = Vec::new();
    let mut claimed_cols: BTreeSet<usize> = BTreeSet::new();
    let mut counts = DedupCounts::default();
    counts.n_in = sorted_cands.len();

    for i in 0..sorted_cands.len() {
        let (tag, region, data_cols, is_imx_scan, full_blosum, covered_blosum, gap_residues_score) = {
            let r = &sorted_cands[i].2;
            (
                r.tag.clone(),
                fmt_region(&r.locus),
                r.locus.msa_cols.clone(),
                r.imx_anchor.is_some(),
                r.score.full_blosum,
                r.score.covered_blosum,
                r.score.gap_residues,
            )
        };
        let overlap = data_cols.intersection(&claimed_cols).count();
        let overlap_thr = (DEDUP_MAX_OVERLAP as f64).min(data_cols.len() as f64 * 0.5);

        let mut mono_violator: Option<usize> = None;
        let mut genomic_violator: Option<usize> = None;
        let mut module_partner: Option<(usize, PairRelation)> = None;

        for &acc_i in &accepted_idx {
            let rel = classify_pair(&sorted_cands[i].2, &sorted_cands[acc_i].2, is_imx_scan);
            match rel.kind {
                PairKind::MonotonicityBreak if mono_violator.is_none() => {
                    mono_violator = Some(acc_i);
                }
                PairKind::GenomicOverlap if genomic_violator.is_none() => {
                    genomic_violator = Some(acc_i);
                }
                PairKind::Module if module_partner.is_none() => {
                    module_partner = Some((acc_i, rel));
                }
                _ => {}
            }
        }

        if debug_level > 0 {
            rejection_log.push(format!(
                "DEBUG_DEDUP {} data_cols={} overlap={} thr={:.0} claimed={} blosum={:.1}%/{:.1}% gap_residues={} scope={} {} region={}",
                tag,
                data_cols.len(),
                overlap,
                overlap_thr,
                claimed_cols.len(),
                full_blosum,
                covered_blosum,
                gap_residues_score,
                scope,
                group_label,
                region,
            ));
        }

        if let Some(violator) = mono_violator {
            let partner_tag = sorted_cands[violator].2.tag.clone();
            rejection_log.push(format!(
                "REJECTED {} reason=monotonicity_break partner={} scope={} region={}",
                tag, partner_tag, scope, region,
            ));
            counts.rej_mono += 1;
            continue;
        }
        if let Some(violator) = genomic_violator {
            let partner_tag = sorted_cands[violator].2.tag.clone();
            rejection_log.push(format!(
                "REJECTED {} reason=genomic_overlap partner={} scope={} region={}",
                tag, partner_tag, scope, region,
            ));
            counts.rej_genomic += 1;
            continue;
        }

        let mut overlap_rej = false;
        if let Some((acc_i, rel)) = module_partner {
            let module_type = module_type_for(&sorted_cands[i].2);
            let acc_tag = sorted_cands[acc_i].2.tag.clone();
            let acc_cols = sorted_cands[acc_i].2.locus.msa_cols.clone();
            sorted_cands[i].2.module = Some(ModuleRel {
                partner_tag: acc_tag.clone(),
                module_type: module_type.to_string(),
                slot_cols: data_cols.clone(),
            });
            sorted_cands[acc_i].2.module = Some(ModuleRel {
                partner_tag: tag.clone(),
                module_type: module_type.to_string(),
                slot_cols: acc_cols,
            });
            rejection_log.push(format!(
                "MODULE_DETECTED {} partner={} type={} dist={}bp scope={} region={}",
                tag, acc_tag, module_type, rel.genomic_dist, scope, region,
            ));
            counts.module += 1;
        } else if (overlap as f64) > overlap_thr {
            rejection_log.push(format!(
                "REJECTED {} reason=slot_duplicate overlap={} scope={} region={}",
                tag, overlap, scope, region,
            ));
            counts.rej_overlap += 1;
            overlap_rej = true;
        }
        if overlap_rej {
            continue;
        }

        accepted_idx.push(i);
        for c in &data_cols {
            claimed_cols.insert(*c);
        }
        counts.accepted += 1;
    }

    (accepted_idx, counts)
}

fn emit_group_summary(
    rejection_log: &mut Vec<String>,
    scope: &str,
    group_label: &str,
    counts: &DedupCounts,
) {
    if counts.n_in > 1
        || counts.module > 0
        || counts.rej_mono > 0
        || counts.rej_genomic > 0
        || counts.rej_overlap > 0
    {
        rejection_log.push(format!(
            "GROUP_SUMMARY in={} accepted={} module={} rej_mono={} rej_genomic={} rej_overlap={} scope={} {}",
            counts.n_in, counts.accepted, counts.module,
            counts.rej_mono, counts.rej_genomic, counts.rej_overlap,
            scope, group_label,
        ));
    }
}

// ===========================================================================
// Cross-cluster genomic dedup -- mirrors _cross_cluster_genomic_dedup
// ===========================================================================

fn cross_cluster_genomic_dedup(
    pending_by_cluster: &mut HashMap<String, Vec<(String, String, Recovery)>>,
    rejection_log: &mut Vec<String>,
) {
    // Flatten (cluster_key, idx_in_cluster, rec_clone) for pairwise comparison.
    let mut all: Vec<(String, usize, Recovery)> = Vec::new();
    for (ck, cands) in pending_by_cluster.iter() {
        for (i, (_, _, r)) in cands.iter().enumerate() {
            all.push((ck.clone(), i, r.clone()));
        }
    }

    fn rank(r: &Recovery) -> (f64, usize) {
        (r.score.full_blosum, count_residues(&r.aa_seq))
    }

    // Cluster keys are "N" (base row) or "N_M" (isoform M of base N).
    fn base_of(ck: &str) -> &str {
        match ck.rsplit_once('_') {
            Some((b, _)) => b,
            None => ck,
        }
    }
    fn is_base_key(ck: &str) -> bool {
        !ck.contains('_')
    }

    let mut losers: HashSet<(String, usize)> = HashSet::new();
    let mut loser_winner: HashMap<(String, usize), Recovery> = HashMap::new();
    let n = all.len();
    for i in 0..n {
        if losers.contains(&(all[i].0.clone(), all[i].1)) {
            continue;
        }
        for j in (i + 1)..n {
            if losers.contains(&(all[j].0.clone(), all[j].1)) {
                continue;
            }
            if all[i].2.locus.scaffold != all[j].2.locus.scaffold {
                continue;
            }
            if all[i].2.locus.strand != all[j].2.locus.strand {
                continue;
            }
            let overlap = interval_overlap_bp(
                all[i].2.locus.bp_start, all[i].2.locus.bp_end,
                all[j].2.locus.bp_start, all[j].2.locus.bp_end,
            );
            if overlap <= 0 {
                continue;
            }
            let len_i = (all[i].2.locus.bp_end + 1).saturating_sub(all[i].2.locus.bp_start);
            let len_j = (all[j].2.locus.bp_end + 1).saturating_sub(all[j].2.locus.bp_start);
            let min_len = len_i.min(len_j);
            if (overlap as f64) < CROSS_CLUSTER_OVERLAP_FRAC * (min_len as f64) {
                continue;
            }
            let ck_i = all[i].2.cluster_key.as_str();
            let ck_j = all[j].2.cluster_key.as_str();
            // Same locus claimed by a base row and one of its own isoform
            // rows: always keep the base copy. A base-keyed gap is fanned out
            // to the cluster's isoforms by propagate_iso_gaps (which sources
            // only from base keys), whereas an isoform-keyed gap is stranded
            // in that one row and can never reach the base or its siblings.
            // The base therefore strictly dominates, so this outranks score.
            let (winner_i, loser_i) = if base_of(ck_i) == base_of(ck_j)
                && is_base_key(ck_i) != is_base_key(ck_j)
            {
                if is_base_key(ck_i) { (i, j) } else { (j, i) }
            } else if rank(&all[i].2) >= rank(&all[j].2) {
                (i, j)
            } else {
                (j, i)
            };
            let loser_key = (all[loser_i].0.clone(), all[loser_i].1);
            losers.insert(loser_key.clone());
            loser_winner.insert(loser_key, all[winner_i].2.clone());
            if loser_i == i {
                break;
            }
        }
    }

    // Emit rejection lines + filter out losers from pending_by_cluster.
    for (ck, cands) in pending_by_cluster.iter_mut() {
        let mut keep: Vec<(String, String, Recovery)> = Vec::with_capacity(cands.len());
        for (i, entry) in cands.drain(..).enumerate() {
            if losers.contains(&(ck.clone(), i)) {
                let region = fmt_region(&entry.2.locus);
                let winner = loser_winner.get(&(ck.clone(), i)).unwrap();
                rejection_log.push(format!(
                    "REJECTED {} reason=duplicate_genomic_locus_across_clusters partner={} cluster={} partner_cluster={} full_blosum={:.1}% partner_full_blosum={:.1}% region={}",
                    entry.2.tag,
                    winner.tag,
                    entry.2.cluster_key,
                    winner.cluster_key,
                    entry.2.score.full_blosum,
                    winner.score.full_blosum,
                    region,
                ));
            } else {
                keep.push(entry);
            }
        }
        *cands = keep;
    }
}

// ===========================================================================
// IMX late kick + slot-sibling enumeration -- mirrors _imx_late_kick /
// _collect_imx_slot_siblings
// ===========================================================================

fn collect_imx_slot_siblings(
    candidate_cluster_key: &str,
    left_flank: &str,
    right_flank: &str,
    imx_slot_node: &str,
    clusters: &HashMap<String, Vec<String>>,
) -> HashSet<String> {
    let mut siblings: HashSet<String> = HashSet::new();
    if clusters.is_empty() {
        return siblings;
    }
    let base_key = match candidate_cluster_key.rsplit_once('_') {
        Some((b, _)) => b,
        None => candidate_cluster_key,
    };
    for (ck, tokens) in clusters.iter() {
        let other_base = match ck.rsplit_once('_') {
            Some((b, _)) => b,
            None => ck.as_str(),
        };
        if other_base != base_key {
            continue;
        }
        let lp = tokens.iter().position(|t| t == left_flank);
        let rp = tokens.iter().position(|t| t == right_flank);
        if let (Some(lp), Some(rp)) = (lp, rp) {
            let (lo, hi) = if lp < rp { (lp, rp) } else { (rp, lp) };
            for t in &tokens[lo + 1..hi] {
                if t != imx_slot_node {
                    siblings.insert(t.clone());
                }
            }
        }
    }
    siblings
}

fn imx_late_kick(
    surviving_gaps: Vec<Recovery>,
    seq_by_token: &HashMap<String, String>,
    clusters: &HashMap<String, Vec<String>>,
    rejection_log: &mut Vec<String>,
    gap_counts: &mut HashMap<&'static str, usize>,
) -> (Vec<Recovery>, HashSet<String>, HashSet<String>) {
    let mut kicked_headers: HashSet<String> = HashSet::new();
    let mut kicked_tags: HashSet<String> = HashSet::new();
    let mut kept: Vec<Recovery> = Vec::with_capacity(surviving_gaps.len());

    for rec in surviving_gaps {
        let Some(imx) = rec.imx_anchor.clone() else {
            kept.push(rec);
            continue;
        };
        if rec.module.is_some() {
            kept.push(rec);
            continue;
        }
        let sibs = collect_imx_slot_siblings(
            &rec.cluster_key,
            &imx.left_flank,
            &imx.right_flank,
            &imx.slot_node,
            clusters,
        );
        let sib_seqs: Vec<(String, String)> = sibs
            .iter()
            .filter_map(|t| seq_by_token.get(t).map(|s| (t.clone(), s.clone())))
            .collect();
        if sib_seqs.is_empty() {
            kept.push(rec);
            continue;
        }
        let sib_lens: Vec<usize> = sib_seqs.iter().map(|(_, s)| count_residues(s)).collect();
        let min_len = *sib_lens.iter().min().unwrap();
        let max_len = *sib_lens.iter().max().unwrap();
        let r_len = count_residues(&rec.aa_seq);
        let region = fmt_region(&rec.locus);

        if (r_len as f64) < (min_len as f64) * 0.8 || (r_len as f64) > (max_len as f64) * 1.2 {
            *gap_counts.entry("kick_imx_length").or_default() += 1;
            kicked_tags.insert(rec.current_tag().to_string());
            kicked_headers.insert(rec.header.clone());
            rejection_log.push(format!(
                "REJECTED {} reason=imx_length cand_len={} range=[{},{}] slot={} cluster={} region={}",
                rec.tag, r_len, min_len, max_len, imx.slot_node, rec.cluster_key, region,
            ));
            continue;
        }

        let mut best_sib = String::new();
        let mut best_sim = 0f64;
        for (t, s) in &sib_seqs {
            let sim = blosum_identity(&rec.aa_seq, s);
            if sim > best_sim {
                best_sim = sim;
                best_sib = t.clone();
            }
        }
        if best_sim < 0.25 {
            *gap_counts.entry("kick_imx_similarity").or_default() += 1;
            kicked_tags.insert(rec.current_tag().to_string());
            kicked_headers.insert(rec.header.clone());
            rejection_log.push(format!(
                "REJECTED {} reason=imx_similarity best_sim={:.1}% thr=25.0% against={} slot={} cluster={} region={}",
                rec.tag,
                best_sim * 100.0,
                best_sib,
                imx.slot_node,
                rec.cluster_key,
                region,
            ));
            continue;
        }

        kept.push(rec);
    }

    (kept, kicked_headers, kicked_tags)
}

// ===========================================================================
// Native supersession -- mirrors _supersede_natives
// ===========================================================================

fn supersede_natives(
    surviving_gaps: &[Recovery],
    gff_nodes: &HashMap<String, GffEntryInput>,
    clusters: &HashMap<String, Vec<String>>,
) -> (HashSet<String>, HashMap<String, String>, HashSet<String>, Vec<String>) {
    let mut superseded_names: HashSet<String> = HashSet::new();
    let mut supersession_pairs: HashMap<String, String> = HashMap::new();
    let mut refind_rejects: HashSet<String> = HashSet::new();
    let mut audit: Vec<String> = Vec::new();

    if clusters.is_empty() {
        return (superseded_names, supersession_pairs, refind_rejects, audit);
    }

    let mut native_tokens: BTreeSet<String> = BTreeSet::new();
    for toks in clusters.values() {
        for t in toks {
            native_tokens.insert(t.clone());
        }
    }

    let imx_loci: Vec<(String, usize, usize, String, String)> = surviving_gaps
        .iter()
        .filter(|r| r.imx_anchor.is_some())
        .map(|r| {
            (
                r.locus.scaffold.clone(),
                r.locus.bp_start, r.locus.bp_end,
                r.locus.strand.clone(),
                r.current_tag().to_string(),
            )
        })
        .collect();
    let plain_loci: Vec<(String, usize, usize, String, String)> = surviving_gaps
        .iter()
        .filter(|r| r.imx_anchor.is_none() && r.current_tag().starts_with("GAP_"))
        .map(|r| {
            (
                r.locus.scaffold.clone(),
                r.locus.bp_start, r.locus.bp_end,
                r.locus.strand.clone(),
                r.current_tag().to_string(),
            )
        })
        .collect();

    let final_recovery_prefix = |t: &str| t.starts_with("FLANK_") || t.starts_with("GAP_");

    fn covers(
        g_scaf: &str, g_start: usize, g_end: usize, g_strand: &str,
        n_scaf: &str, n_start: usize, n_end: usize, n_strand: &str, n_len: usize,
    ) -> Option<(bool, i64)> {
        if g_scaf != n_scaf || g_strand != n_strand {
            return None;
        }
        let contained = n_start >= g_start && n_end <= g_end;
        let overlap = interval_overlap_bp(g_start, g_end, n_start, n_end);
        if contained || (overlap as f64) / (n_len as f64) >= NATIVE_OVERLAP_FRAC {
            Some((contained, overlap))
        } else {
            None
        }
    }

    for n_tok in &native_tokens {
        if final_recovery_prefix(n_tok) {
            continue;
        }
        let Some(n_entry) = gff_nodes.get(n_tok) else { continue };
        let n_len = (n_entry.end + 1).saturating_sub(n_entry.start);
        if n_len == 0 {
            continue;
        }
        // IMX supersession: stop on first hit.
        for g in &imx_loci {
            if let Some((contained, overlap)) =
                covers(&g.0, g.1, g.2, &g.3, &n_entry.scaffold, n_entry.start, n_entry.end, &n_entry.strand, n_len)
            {
                superseded_names.insert(n_tok.clone());
                supersession_pairs.insert(n_tok.clone(), g.4.clone());
                audit.push(format!(
                    "SUPERSEDED_NATIVE {} superseded_by={} native={}:{}-{}({}) gap={}:{}-{}({}) overlap_bp={}/{} contained={}",
                    n_tok, g.4,
                    n_entry.scaffold, n_entry.start, n_entry.end, n_entry.strand,
                    g.0, g.1, g.2, g.3,
                    overlap, n_len, contained,
                ));
                break;
            }
        }
        // Plain gap refind rejection: scan all.
        for g in &plain_loci {
            if let Some((contained, overlap)) =
                covers(&g.0, g.1, g.2, &g.3, &n_entry.scaffold, n_entry.start, n_entry.end, &n_entry.strand, n_len)
            {
                refind_rejects.insert(g.4.clone());
                audit.push(format!(
                    "REJECTED_REFIND {} duplicates_native={} native={}:{}-{}({}) gap={}:{}-{}({}) overlap_bp={}/{} contained={}",
                    g.4, n_tok,
                    n_entry.scaffold, n_entry.start, n_entry.end, n_entry.strand,
                    g.0, g.1, g.2, g.3,
                    overlap, n_len, contained,
                ));
            }
        }
    }

    (superseded_names, supersession_pairs, refind_rejects, audit)
}

// ===========================================================================
// split_calc ports: build_ref_columns, calculate_split, is_module_pair,
// detect_module_groups
// ===========================================================================

#[derive(Default, Clone, Debug)]
struct RefColumn {
    counts: HashMap<u8, i32>, // includes gaps
}

impl RefColumn {
    fn add(&mut self, b: u8) {
        *self.counts.entry(b).or_default() += 1;
    }
    fn get(&self, b: u8) -> i32 {
        *self.counts.get(&b).unwrap_or(&0)
    }
    fn contains(&self, b: u8) -> bool {
        self.counts.contains_key(&b)
    }
}

fn find_index_pair_bytes(seq: &[u8], gap: u8) -> (usize, usize) {
    let s = seq.iter().position(|&c| c != gap).unwrap_or(0);
    let e = seq.iter().rposition(|&c| c != gap).map(|i| i + 1).unwrap_or(0);
    (s, e)
}

fn build_ref_columns(ref_seqs: &[String]) -> Vec<RefColumn> {
    let aln_len = ref_seqs.iter().map(|s| s.len()).max().unwrap_or(0);
    let mut cols: Vec<RefColumn> = vec![RefColumn::default(); aln_len];
    for seq in ref_seqs {
        let bytes = seq.as_bytes();
        let (s, e) = find_index_pair_bytes(bytes, b'-');
        for i in s..e {
            if i < cols.len() {
                cols[i].add(bytes[i]);
            }
        }
    }
    cols
}

fn _score_letter(letter: u8, ref_col: Option<&RefColumn>) -> i32 {
    let mut score = match ref_col {
        Some(rc) if rc.contains(letter) => rc.get(letter),
        _ => -1,
    };
    if letter == b'-' {
        score -= 1;
    }
    score
}

fn calculate_split(
    seq_a: &str,
    seq_b: &str,
    overlap_start: usize,
    overlap_end: usize,
    ref_columns: &[RefColumn],
) -> usize {
    if overlap_end <= overlap_start {
        return overlap_start;
    }
    let l = overlap_end - overlap_start;
    let mut contrib_a = vec![0i32; l];
    let mut contrib_b = vec![0i32; l];
    let a_bytes = seq_a.as_bytes();
    let b_bytes = seq_b.as_bytes();
    for j in 0..l {
        let col = overlap_start + j;
        let ref_col = if col < ref_columns.len() {
            Some(&ref_columns[col])
        } else {
            None
        };
        let ca = if col < a_bytes.len() { a_bytes[col] } else { b'-' };
        let cb = if col < b_bytes.len() { b_bytes[col] } else { b'-' };
        contrib_a[j] = _score_letter(ca, ref_col);
        contrib_b[j] = _score_letter(cb, ref_col);
    }
    let mut pref_a = vec![0i32; l + 1];
    let mut pref_b = vec![0i32; l + 1];
    for i in 0..l {
        pref_a[i + 1] = pref_a[i] + contrib_a[i];
        pref_b[i + 1] = pref_b[i] + contrib_b[i];
    }
    let total_b = pref_b[l];
    let mut highest_score = -1i32;
    let mut highest_pos = overlap_start;
    for i in overlap_start..=overlap_end {
        let j = i - overlap_start;
        let score = pref_a[j] + (total_b - pref_b[j]);
        if score >= highest_score {
            highest_score = score;
            highest_pos = i;
        }
    }
    highest_pos
}

fn is_module_pair_aligned(
    seq_a: &str,
    seq_b: &str,
    parent_a: &str,
    gpos_a: usize,
    parent_b: &str,
    gpos_b: usize,
) -> bool {
    if parent_a != parent_b {
        return false;
    }
    let dist = (gpos_a as i64 - gpos_b as i64).unsigned_abs() as usize;
    if dist < MODULE_MIN_INTRON_BP_SPLIT {
        return false;
    }
    let a_bytes = seq_a.as_bytes();
    let b_bytes = seq_b.as_bytes();
    let (a_start, a_end) = find_index_pair_bytes(a_bytes, b'-');
    let (b_start, b_end) = find_index_pair_bytes(b_bytes, b'-');
    let coords_a: BTreeSet<usize> = (a_start..a_end)
        .filter(|&i| i < a_bytes.len() && a_bytes[i] != b'-')
        .collect();
    let coords_b: BTreeSet<usize> = (b_start..b_end)
        .filter(|&i| i < b_bytes.len() && b_bytes[i] != b'-')
        .collect();
    if coords_a.is_empty() || coords_b.is_empty() {
        return false;
    }
    let smaller = coords_a.len().min(coords_b.len());
    let larger = coords_a.len().max(coords_b.len());
    if (smaller as f64) / (larger as f64) < MODULE_MIN_SIZE_RATIO_SPLIT {
        return false;
    }
    let overlap = coords_a.intersection(&coords_b).count();
    (overlap as f64) / (smaller as f64) >= MODULE_MIN_OVERLAP_SPLIT
}

fn detect_module_groups(
    members: &[String],
    records_map: &HashMap<String, String>,
    header_to_gff: &dyn Fn(&str) -> Option<(String, usize, String)>,
) -> Vec<Vec<String>> {
    let mut parent: HashMap<String, String> = HashMap::with_capacity(members.len());
    for h in members {
        parent.insert(h.clone(), h.clone());
    }
    fn find(parent: &mut HashMap<String, String>, h: &str) -> String {
        let mut cur = h.to_string();
        loop {
            let p = parent.get(&cur).cloned().unwrap_or_else(|| cur.clone());
            if p == cur {
                return cur;
            }
            // Path compression: parent[cur] = parent[parent[cur]]
            let gp = parent.get(&p).cloned().unwrap_or_else(|| p.clone());
            parent.insert(cur.clone(), gp.clone());
            cur = gp;
        }
    }
    fn union(parent: &mut HashMap<String, String>, a: &str, b: &str) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            parent.insert(ra, rb);
        }
    }

    for i in 0..members.len() {
        let hi = &members[i];
        let Some(gi) = header_to_gff(hi) else { continue };
        let Some(seq_i) = records_map.get(hi) else { continue };
        for j in (i + 1)..members.len() {
            let hj = &members[j];
            let Some(gj) = header_to_gff(hj) else { continue };
            let Some(seq_j) = records_map.get(hj) else { continue };
            if is_module_pair_aligned(seq_i, seq_j, &gi.0, gi.1, &gj.0, gj.1) {
                union(&mut parent, hi, hj);
            }
        }
    }

    let mut groups: HashMap<String, Vec<String>> = HashMap::new();
    for h in members {
        let r = find(&mut parent, h);
        groups.entry(r).or_default().push(h.clone());
    }
    let mut out: Vec<Vec<String>> = Vec::new();
    for (_, group) in groups {
        if group.len() < 2 {
            continue;
        }
        let is_minus = group.iter().any(|h| {
            header_to_gff(h).map(|(_, _, s)| s == "-").unwrap_or(false)
        });
        let mut sorted = group;
        sorted.sort_by(|a, b| {
            let ga = header_to_gff(a).map(|(_, gp, _)| gp).unwrap_or(0);
            let gb = header_to_gff(b).map(|(_, gp, _)| gp).unwrap_or(0);
            if is_minus { gb.cmp(&ga) } else { ga.cmp(&gb) }
        });
        out.push(sorted);
    }
    out
}

// is_same_kmer: hamming distance over equal-length slices.
fn is_same_kmer(a: &str, b: &str) -> bool {
    if a.len() != b.len() {
        return false;
    }
    a.as_bytes()
        .iter()
        .zip(b.as_bytes())
        .all(|(x, y)| x == y)
}

// ===========================================================================
// apply_cull_trim / apply_drop
// ===========================================================================

fn apply_cull_trim(rec: &mut Recovery, drops: &HashSet<usize>) {
    if drops.is_empty() {
        return;
    }
    let orig_aa_len = rec.aa_seq.len();
    let mut left_trim = 0usize;
    while left_trim < orig_aa_len && drops.contains(&left_trim) {
        left_trim += 1;
    }
    if left_trim == orig_aa_len {
        return;
    }
    let mut right_trim = 0usize;
    while right_trim < orig_aa_len && drops.contains(&(orig_aa_len - 1 - right_trim)) {
        right_trim += 1;
    }
    if left_trim == 0 && right_trim == 0 {
        return;
    }
    let new_locus = rec.locus.trim_aa(left_trim, right_trim);
    if new_locus.bp_start > new_locus.bp_end {
        return;
    }
    rec.locus = new_locus;
}

fn apply_drop(
    hdr: &str,
    drop_cols: &BTreeSet<usize>,
    records_map: &mut HashMap<String, String>,
    nt_cull_map: &mut HashMap<String, HashSet<usize>>,
    orig_aa_len_by_header: &HashMap<String, usize>,
) -> Result<usize, String> {
    if drop_cols.is_empty() {
        return Ok(0);
    }
    let Some(seq) = records_map.get(hdr).cloned() else {
        return Ok(0);
    };
    let seq_bytes = seq.as_bytes();
    let residue_cols: Vec<usize> = drop_cols
        .iter()
        .filter(|&&c| c < seq_bytes.len() && !is_gap(seq_bytes[c]))
        .copied()
        .collect();

    let mut chars: Vec<u8> = seq_bytes.to_vec();
    for &c in drop_cols {
        if c < chars.len() {
            chars[c] = b'-';
        }
    }
    records_map.insert(
        hdr.to_string(),
        String::from_utf8(chars).unwrap_or(seq.clone()),
    );

    if residue_cols.is_empty() {
        return Ok(0);
    }

    let orig_len = *orig_aa_len_by_header.get(hdr).unwrap_or(&0);
    if orig_len == 0 {
        return Err(format!(
            "apply_drop: no orig AA len for {}; NT drops cannot sync",
            hdr
        ));
    }
    let existing = nt_cull_map.get(hdr).cloned().unwrap_or_default();
    let kept_sorted: Vec<usize> = (0..orig_len).filter(|i| !existing.contains(i)).collect();
    let residue_col_set: BTreeSet<usize> = residue_cols.into_iter().collect();
    let mut aligned_idx = 0usize;
    let mut new_drops: Vec<usize> = Vec::new();
    for (col, &ch) in seq_bytes.iter().enumerate() {
        if is_gap(ch) {
            continue;
        }
        if residue_col_set.contains(&col) && aligned_idx < kept_sorted.len() {
            new_drops.push(kept_sorted[aligned_idx]);
        }
        aligned_idx += 1;
    }
    let count = new_drops.len();
    if !new_drops.is_empty() {
        let entry = nt_cull_map.entry(hdr.to_string()).or_default();
        for d in new_drops {
            entry.insert(d);
        }
    }
    Ok(count)
}

// ===========================================================================
// run_final_split
// ===========================================================================

#[allow(clippy::too_many_arguments)]
fn run_final_split(
    gene_tag: &str,
    records: Vec<(String, String)>,
    cluster_by_header: &HashMap<String, String>,
    orig_aa_len_by_header: &HashMap<String, usize>,
    nt_cull_map: &mut HashMap<String, HashSet<usize>>,
    gff_nodes: &HashMap<String, GffEntryInput>,
    skip_exon_split: bool,
    native_cluster_members: &HashMap<String, Vec<String>>,
    debug: bool,
) -> (
    Vec<(String, String)>,
    HashSet<String>,
    HashMap<String, (usize, usize)>,
    Vec<String>,
) {
    let mut split_kicked: HashSet<String> = HashSet::new();
    let mut split_calc_log: Vec<String> = Vec::new();
    let mut split_trim: HashMap<String, (usize, usize)> = HashMap::new();

    if skip_exon_split {
        if debug {
            split_calc_log.push(format!("gene={} SKIPPED (--skip-exon-split)", gene_tag));
        }
        return (records, split_kicked, split_trim, split_calc_log);
    }

    if cluster_by_header.is_empty() || records.is_empty() {
        return (records, split_kicked, split_trim, split_calc_log);
    }

    let mut records_map: HashMap<String, String> =
        records.iter().map(|(h, s)| (h.clone(), s.clone())).collect();
    let ref_only_seqs: Vec<String> = records
        .iter()
        .filter(|(h, _)| h.ends_with('.'))
        .map(|(_, s)| s.clone())
        .collect();
    let ref_columns = build_ref_columns(&ref_only_seqs);

    // cluster_members: cluster_key -> [headers including natives + recoveries]
    let mut cluster_members: HashMap<String, Vec<String>> = HashMap::new();
    for (hdr, ck) in cluster_by_header {
        if records_map.contains_key(hdr) {
            cluster_members.entry(ck.clone()).or_default().push(hdr.clone());
        }
    }
    for (ck, hdrs) in native_cluster_members {
        let entry = cluster_members.entry(ck.clone()).or_default();
        let existing: HashSet<String> = entry.iter().cloned().collect();
        for hdr in hdrs {
            if records_map.contains_key(hdr) && !existing.contains(hdr) {
                entry.push(hdr.clone());
            }
        }
    }

    let header_to_gff = |hdr: &str| -> Option<(String, usize, String)> {
        let tag = tag_from_header(hdr);
        let entry = gff_nodes.get(tag)?;
        Some((entry.scaffold.clone(), entry.start, entry.strand.clone()))
    };

    let mut gene_cluster_sections: Vec<String> = Vec::new();

    let mut cluster_keys: Vec<&String> = cluster_members.keys().collect();
    cluster_keys.sort();
    for ck in cluster_keys {
        let cluster_all = cluster_members.get(ck).cloned().unwrap_or_default();
        let groups = detect_module_groups(&cluster_all, &records_map, &header_to_gff);

        let mut module_non_reps: HashSet<String> = HashSet::new();
        for group in &groups {
            for sib in group.iter().skip(1) {
                module_non_reps.insert(sib.clone());
            }
        }

        let mut members: Vec<String> = cluster_all
            .iter()
            .filter(|m| !module_non_reps.contains(*m))
            .cloned()
            .collect();
        if members.len() < 2 && groups.is_empty() {
            continue;
        }

        members.sort_by(|a, b| {
            let (first_a, _) = seq_residue_bounds(records_map.get(a).map(String::as_str).unwrap_or(""));
            let (first_b, _) = seq_residue_bounds(records_map.get(b).map(String::as_str).unwrap_or(""));
            (first_a, a).cmp(&(first_b, b))
        });

        let mut cluster_block: Vec<String> = Vec::new();
        let n = members.len();
        let mut i = 0;
        'outer: while i < n {
            let a_hdr = members[i].clone();
            if split_kicked.contains(&a_hdr) {
                i += 1;
                continue;
            }
            for j in (i + 1)..n {
                let b_hdr = members[j].clone();
                if split_kicked.contains(&b_hdr) {
                    continue;
                }
                let seq_a = records_map.get(&a_hdr).cloned().unwrap_or_default();
                let seq_b = records_map.get(&b_hdr).cloned().unwrap_or_default();
                let data_cols_a = data_cols_set(&seq_a);
                let data_cols_b = data_cols_set(&seq_b);
                if data_cols_a.intersection(&data_cols_b).next().is_none() {
                    continue;
                }
                let (a_first_i, a_last_i) = seq_residue_bounds(&seq_a);
                let (b_first_i, b_last_i) = seq_residue_bounds(&seq_b);
                if a_first_i < 0 || b_first_i < 0 {
                    continue;
                }
                let a_first = a_first_i as usize;
                let a_last = a_last_i as usize;
                let b_first = b_first_i as usize;
                let b_last = b_last_i as usize;
                let overlap_start = a_first.max(b_first);
                let overlap_end = a_last.min(b_last) + 1;
                if overlap_end <= overlap_start {
                    continue;
                }

                let a_tag = tag_from_header(&a_hdr).to_string();
                let b_tag = tag_from_header(&b_hdr).to_string();
                let gff_a = header_to_gff(&a_hdr);
                let gff_b = header_to_gff(&b_hdr);
                if let (Some(ga), Some(gb)) = (gff_a, gff_b) {
                    if is_module_pair_aligned(&seq_a, &seq_b, &ga.0, ga.1, &gb.0, gb.1) {
                        if debug {
                            cluster_block.push(format!(
                                "gene={} cluster={} MODULE_SKIP {} vs {}",
                                gene_tag, ck, a_tag, b_tag,
                            ));
                        }
                        continue;
                    }
                }

                let kmer_a: &str = &seq_a[overlap_start..overlap_end.min(seq_a.len())];
                let kmer_b: &str = &seq_b[overlap_start..overlap_end.min(seq_b.len())];
                if is_same_kmer(kmer_a, kmer_b) {
                    continue;
                }

                let split_k = calculate_split(&seq_a, &seq_b, overlap_start, overlap_end, &ref_columns);

                let is_recovery_prefix = |t: &str| {
                    t.starts_with("GAP_") || t.starts_with("FLANK_")
                        || t.starts_with("_TMPGAP_") || t.starts_with("_TMPFLANK_")
                };
                let a_native = !is_recovery_prefix(&a_tag);
                let b_native = !is_recovery_prefix(&b_tag);

                let a_drop_cols: BTreeSet<usize> = (split_k..=a_last).collect();
                let b_drop_cols: BTreeSet<usize> = (b_first..split_k).collect();

                let a_dropped = if a_native {
                    0
                } else {
                    match apply_drop(&a_hdr, &a_drop_cols, &mut records_map, nt_cull_map, orig_aa_len_by_header) {
                        Ok(n) => n,
                        Err(_) => 0,
                    }
                };
                let b_dropped = if b_native {
                    0
                } else {
                    match apply_drop(&b_hdr, &b_drop_cols, &mut records_map, nt_cull_map, orig_aa_len_by_header) {
                        Ok(n) => n,
                        Err(_) => 0,
                    }
                };
                if a_dropped > 0 {
                    let e = split_trim.entry(a_tag.clone()).or_insert((0, 0));
                    e.1 += a_dropped;
                }
                if b_dropped > 0 {
                    let e = split_trim.entry(b_tag.clone()).or_insert((0, 0));
                    e.0 += b_dropped;
                }

                let a_remaining = count_residues(records_map.get(&a_hdr).map(String::as_str).unwrap_or(""));
                let b_remaining = count_residues(records_map.get(&b_hdr).map(String::as_str).unwrap_or(""));

                if debug {
                    cluster_block.push(format!(
                        "gene={} cluster={} PAIR {} vs {} split_at={} trimmed={}/{} remaining={}/{}",
                        gene_tag, ck, a_tag, b_tag, split_k, a_dropped, b_dropped,
                        a_remaining, b_remaining,
                    ));
                }

                let mut a_kicked = false;
                if a_remaining < MIN_AA_AFTER_SPLIT && !a_native {
                    split_kicked.insert(a_hdr.clone());
                    a_kicked = true;
                    if debug {
                        cluster_block.push(format!(
                            "gene={} cluster={} KICK {} reason=<{}_aa_after_split remaining={}",
                            gene_tag, ck, a_tag, MIN_AA_AFTER_SPLIT, a_remaining,
                        ));
                    }
                }
                if b_remaining < MIN_AA_AFTER_SPLIT && !b_native {
                    split_kicked.insert(b_hdr.clone());
                    if debug {
                        cluster_block.push(format!(
                            "gene={} cluster={} KICK {} reason=<{}_aa_after_split remaining={}",
                            gene_tag, ck, b_tag, MIN_AA_AFTER_SPLIT, b_remaining,
                        ));
                    }
                }
                if a_kicked {
                    i += 1;
                    continue 'outer;
                }
            }
            i += 1;
        }

        // Module sibling trim against the rep.
        for group in &groups {
            let Some(rep_hdr) = group.first() else { continue };
            if !records_map.contains_key(rep_hdr) {
                continue;
            }
            let rep_seq = records_map.get(rep_hdr).cloned().unwrap_or_default();
            let (rep_first_i, rep_last_i) = seq_residue_bounds(&rep_seq);
            if rep_first_i < 0 {
                continue;
            }
            let rep_start = rep_first_i as usize;
            let rep_end = (rep_last_i as usize) + 1;
            let rep_tag = tag_from_header(rep_hdr).to_string();
            let rep_kicked = split_kicked.contains(rep_hdr);

            for sib_hdr in group.iter().skip(1) {
                if split_kicked.contains(sib_hdr) {
                    continue;
                }
                let Some(sib_seq) = records_map.get(sib_hdr).cloned() else { continue };
                let (sib_first_i, sib_last_i) = seq_residue_bounds(&sib_seq);
                if sib_first_i < 0 {
                    continue;
                }
                let sib_first = sib_first_i as usize;
                let sib_end = (sib_last_i as usize) + 1;
                let mut trim_cols: BTreeSet<usize> = BTreeSet::new();
                for c in sib_first..(sib_end.min(rep_start)) {
                    trim_cols.insert(c);
                }
                for c in sib_first.max(rep_end)..sib_end {
                    trim_cols.insert(c);
                }
                if trim_cols.is_empty() {
                    continue;
                }

                let before_sib = sib_seq.clone();
                let dropped = apply_drop(sib_hdr, &trim_cols, &mut records_map, nt_cull_map, orig_aa_len_by_header)
                    .unwrap_or(0);
                let remaining = count_residues(records_map.get(sib_hdr).map(String::as_str).unwrap_or(""));

                let before_bytes = before_sib.as_bytes();
                let left_trim = trim_cols
                    .iter()
                    .filter(|&&c| c < rep_start && c < before_bytes.len() && !is_gap(before_bytes[c]))
                    .count();
                let right_trim = dropped.saturating_sub(left_trim);
                let tag = tag_from_header(sib_hdr).to_string();
                if left_trim > 0 {
                    let e = split_trim.entry(tag.clone()).or_insert((0, 0));
                    e.0 += left_trim;
                }
                if right_trim > 0 {
                    let e = split_trim.entry(tag.clone()).or_insert((0, 0));
                    e.1 += right_trim;
                }
                if debug {
                    let rep_kick_note = if rep_kicked { " rep_kicked=True" } else { "" };
                    cluster_block.push(format!(
                        "gene={} cluster={} MODULE_TRIM rep={} sib={} dropped={} left={} right={} window=[{},{}) remaining={}{}",
                        gene_tag, ck, rep_tag, tag, dropped, left_trim, right_trim,
                        rep_start, rep_end, remaining, rep_kick_note,
                    ));
                }
                if remaining < MIN_AA_AFTER_SPLIT {
                    split_kicked.insert(sib_hdr.clone());
                    if debug {
                        cluster_block.push(format!(
                            "gene={} cluster={} KICK {} reason=<{}_aa_after_module_trim remaining={}",
                            gene_tag, ck, tag, MIN_AA_AFTER_SPLIT, remaining,
                        ));
                    }
                }
            }
        }

        if !cluster_block.is_empty() {
            gene_cluster_sections.extend(cluster_block);
        }
    }

    if !gene_cluster_sections.is_empty() {
        split_calc_log.extend(gene_cluster_sections);
    }

    let mut new_records: Vec<(String, String)> = Vec::with_capacity(records.len());
    for (h, s) in &records {
        if split_kicked.contains(h) {
            continue;
        }
        new_records.push((h.clone(), records_map.get(h).cloned().unwrap_or(s.clone())));
    }

    // Drop tags whose trim sums to (0, 0).
    let split_trim: HashMap<String, (usize, usize)> = split_trim
        .into_iter()
        .filter(|(_, (l, r))| *l > 0 || *r > 0)
        .collect();

    (new_records, split_kicked, split_trim, split_calc_log)
}

// ===========================================================================
// genomic sort key (for final fasta ordering)
// ===========================================================================

fn genomic_sort_key(
    header: &str,
    gff_nodes: &HashMap<String, GffEntryInput>,
) -> (String, usize, usize, i64, String) {
    let fields: Vec<&str> = header.split('|').collect();
    let node_field = fields.get(3).copied().unwrap_or("");
    let entry = gff_nodes.get(node_field);
    let (parent, start, end) = match entry {
        Some(e) => (e.scaffold.clone(), e.start, e.end),
        None => (String::new(), 0, 0),
    };
    let primary = node_field
        .split("&&").next().unwrap_or("")
        .split('_').next().unwrap_or("");
    let node_int: i64 = primary.parse().unwrap_or(0);
    (parent, start, end, node_int, node_field.to_string())
}

/// Restore `*` (stop codons) on the native rows of the freshly hmm_aligned
/// records.
///
/// Natives are folded into the hmm_align reference template, which
/// `aligner.rs` masks `*`->`X` so hmmbuild will accept them.  `--mapali` then
/// copies those rows into the output verbatim (only re-gapping to absorb
/// candidate insert columns), so the emitted native rows carry `X` where they
/// had a stop.  `--mapali` preserves residue identity and order, so we walk
/// each native's original stop offsets and flip just those positions back.
///
/// Must run *before* cull, while each native row's non-gap count still equals
/// its ungapped source length.  Only original `*` positions are touched, so
/// genuine `X` residues are left as-is.
fn restore_native_stops(
    aligned: &mut [(String, String)],
    natives_aa: &[(String, String)],
) {
    // header -> residue indices (in ungapped order) that were stop codons.
    let stop_positions: HashMap<&str, Vec<usize>> = natives_aa
        .iter()
        .filter_map(|(h, s)| {
            let stops: Vec<usize> = s
                .bytes()
                .filter(|&b| b != b'-' && b != b'.')
                .enumerate()
                .filter_map(|(i, b)| (b == b'*').then_some(i))
                .collect();
            (!stops.is_empty()).then(|| (h.as_str(), stops))
        })
        .collect();

    if stop_positions.is_empty() {
        return;
    }

    for (header, seq) in aligned.iter_mut() {
        let Some(stops) = stop_positions.get(header.as_str()) else { continue };
        let mut stop_iter = stops.iter().copied().peekable();
        let mut residue_idx = 0usize;
        // SAFETY: we only overwrite individual ASCII `X` bytes with ASCII `*`,
        // which keeps the string valid UTF-8 (protein rows are ASCII).
        let bytes = unsafe { seq.as_mut_vec() };
        for b in bytes.iter_mut() {
            if *b == b'-' {
                continue;
            }
            if stop_iter.peek() == Some(&residue_idx) {
                *b = b'*';
                stop_iter.next();
            }
            residue_idx += 1;
        }
    }
}

// ===========================================================================
// Top-level: exonfinder_process_gene
// ===========================================================================

#[pyfunction]
#[pyo3(signature = (
    gene_key,
    aa_file,
    refs_aa,
    natives_aa,
    flank_inputs,
    gap_inputs,
    nt_seqs,
    gff_nodes,
    clusters_for_gene,
    taxa = None,
    tmpdir = None,
    debug_level = 0,
    skip_exon_split = false,
    disable_column_cull = false,
    log_dir = None,
))]
pub fn exonfinder_process_gene(
    py: Python<'_>,
    gene_key: String,
    aa_file: String,
    refs_aa: Vec<(String, String)>,
    natives_aa: Vec<(String, String)>,
    flank_inputs: Vec<FlankInput>,
    gap_inputs: Vec<GapInput>,
    nt_seqs: HashMap<String, String>,
    gff_nodes: HashMap<String, GffEntryInput>,
    clusters_for_gene: HashMap<String, Vec<String>>,
    taxa: Option<String>,
    tmpdir: Option<String>,
    debug_level: i32,
    skip_exon_split: bool,
    disable_column_cull: bool,
    log_dir: Option<String>,
) -> PyResult<Py<PyDict>> {
    let t_total = Instant::now();
    let mut rejection_log: Vec<String> = Vec::new();
    let mut recovered_log: Vec<String> = Vec::new();
    let mut worker_gff_lines: Vec<String> = Vec::new();

    let gene_name = aa_file.split('.').next().unwrap_or(&aa_file).to_string();

    // -------------------------------------------------------------------
    // 1. Reference template + candidate list for hmm_align.
    //
    // The template (hmmbuild input *and* --mapali) is the resolve output
    // MSA: the orthoset refs plus the already-aligned natives.  --mapali maps
    // an existing alignment onto the model *without* realigning it, so the
    // natives come back frozen in their resolve columns and only the new
    // flank/gap stubs are aligned fresh.  Building the profile from
    // refs+natives also gives those stubs a data-informed model for free.
    //
    // aligner.rs masks `*`->`X` on the template so hmmbuild accepts it; the
    // native stops are restored right after the align call (see
    // restore_native_stops) so the emitted AA matches the resolve input.
    // -------------------------------------------------------------------
    if flank_inputs.is_empty() && gap_inputs.is_empty() {
        // No candidates: short-circuit, return the input fasta unchanged
        // (refs + already-aligned natives).
        return empty_result(py, &gene_key, &refs_aa, &natives_aa);
    }

    let mut mapali: Vec<(String, String)> =
        Vec::with_capacity(refs_aa.len() + natives_aa.len());
    mapali.extend(refs_aa.iter().cloned());
    mapali.extend(natives_aa.iter().cloned());

    let mut cands_aa: Vec<(String, String)> =
        Vec::with_capacity(flank_inputs.len() + gap_inputs.len());
    // Flank stubs: one per flank input.
    for fi in &flank_inputs {
        cands_aa.push((fi.header.clone(), fi.aa_seq.clone()));
    }
    // Gap stubs: deduplicated by header (Python keeps the first aa_seq).
    let mut gap_header_seen: HashSet<String> = HashSet::new();
    for gi in &gap_inputs {
        if gap_header_seen.insert(gi.header.clone()) {
            cands_aa.push((gi.header.clone(), gi.aa_seq.clone()));
        }
    }

    // -------------------------------------------------------------------
    // 2. hmm_align (re-uses existing pyfunction; pass our py token).
    // -------------------------------------------------------------------
    let t_hmm = Instant::now();
    let mut aligned = hmm_align(
        py,
        cands_aa,
        mapali,
        tmpdir.clone(),
        Some(gene_name.clone()),
        taxa.clone(),
    )?;
    // Natives were `*`->`X` masked for hmmbuild; put their stop codons back
    // before cull so the emitted AA matches the resolve input residues.
    restore_native_stops(&mut aligned, &natives_aa);
    let t_hmm = t_hmm.elapsed().as_secs_f64();

    // -------------------------------------------------------------------
    // 3. cull_columns (unless disabled).
    // -------------------------------------------------------------------
    let mut cull_nt_seqs: HashMap<String, String> = nt_seqs.clone();
    for fi in &flank_inputs {
        cull_nt_seqs.insert(fi.header.clone(), fi.nt_seq.clone());
    }
    let mut gap_nt_seen: HashSet<String> = HashSet::new();
    for gi in &gap_inputs {
        if gap_nt_seen.insert(gi.header.clone()) {
            cull_nt_seqs.insert(gi.header.clone(), gi.nt_seq.clone());
        }
    }

    let t_cull = Instant::now();
    let (mut culled_records, nt_cull_map) = if disable_column_cull {
        (aligned, HashMap::new())
    } else {
        let (cr, nc, _gff_cull, _gff_intron, _exon_split) = cull_columns(
            aligned,
            ".",
            0.33,
            10,
            12,
            25,
            true,
            0.5,
            10,
            1.0,
            Some(cull_nt_seqs.clone()),
            debug_level,
            log_dir.clone(),
        )?;
        (cr, nc)
    };
    let t_cull = t_cull.elapsed().as_secs_f64();

    // -------------------------------------------------------------------
    // 4. Low-complexity filter (call inner directly).
    // -------------------------------------------------------------------
    let t_flank_total = Instant::now();
    let culled_headers: HashSet<String> = culled_records
        .iter()
        .filter(|(h, _)| !h.ends_with('.'))
        .map(|(h, _)| h.clone())
        .collect();
    let mut lc_kicked: HashSet<String> = HashSet::new();
    for h in &culled_headers {
        if let Some(nt) = cull_nt_seqs.get(h) {
            if crate::is_low_complexity_nt(
                nt,
                60, 15, 2.5, 3.75, 30, 0.5, 100, 4,
            ) {
                lc_kicked.insert(h.clone());
            }
        }
    }
    for h in &lc_kicked {
        rejection_log.push(format!("REJECTED {} reason=low_complexity_nt", h));
        cull_nt_seqs.remove(h);
    }
    if !lc_kicked.is_empty() {
        culled_records.retain(|(h, _)| !lc_kicked.contains(h));
    }

    // -------------------------------------------------------------------
    // 5. Build Recovery skeletons.
    // -------------------------------------------------------------------
    let aln_len = culled_records.first().map(|(_, s)| s.len()).unwrap_or(0);

    // header -> per-header gap-input list
    let mut gap_inputs_by_header: HashMap<String, Vec<GapInput>> = HashMap::new();
    for gi in gap_inputs {
        gap_inputs_by_header
            .entry(gi.header.clone())
            .or_default()
            .push(gi);
    }

    // Flank Recovery objects keyed by header (one each).
    let mut flank_rec_by_header: HashMap<String, Recovery> = HashMap::new();
    for fi in flank_inputs {
        let rec = Recovery {
            tag: fi.tag.clone(),
            header: fi.header.clone(),
            aa_seq: fi.aa_seq.clone(),
            nt_seq: fi.nt_seq.clone(),
            locus: Locus {
                scaffold: fi.scaffold.clone(),
                strand: fi.strand.clone(),
                bp_start: fi.hit_start,
                bp_end: fi.hit_end,
                msa_cols: BTreeSet::new(),
                frame: fi.frame,
            },
            score: ScoreInfo {
                full_blosum: 0.0,
                covered_blosum: 0.0,
                gap_residues: 0,
                win_start: fi.gap_start,
                win_end: fi.gap_end,
            },
            kind: RecoveryKind::Flank,
            cluster_key: fi.cluster_key.clone(),
            gene_key: fi.gene_key.clone(),
            flank_anchor: Some(FlankAnchor {
                adj_node: fi.node_name.clone(),
                is_leading: fi.is_leading,
            }),
            gap_anchor: None,
            imx_anchor: None,
            module: None,
            final_tag: None,
        };
        flank_rec_by_header.insert(fi.header.clone(), rec);
    }

    // Apply cull trims to flank recoveries.
    let nt_drops_view: HashMap<String, HashSet<usize>> = nt_cull_map.clone();
    for (header, _) in &culled_records {
        if header.ends_with('.') {
            continue;
        }
        if let Some(rec) = flank_rec_by_header.get_mut(header) {
            if let Some(drops) = nt_drops_view.get(header) {
                apply_cull_trim(rec, drops);
            }
        }
    }

    // -------------------------------------------------------------------
    // 6. Flank scoring + dedup.
    // -------------------------------------------------------------------
    let mut surviving_flanks: Vec<Recovery> = Vec::new();
    if !flank_rec_by_header.is_empty() {
        let ref_seq_strings: Vec<String> = culled_records
            .iter()
            .filter(|(h, _)| h.ends_with('.'))
            .map(|(_, s)| s.clone())
            .collect();
        let ref_seqs_bytes: Vec<&[u8]> =
            ref_seq_strings.iter().map(|s| s.as_bytes()).collect();

        // token -> seq (excluding the flank stubs themselves)
        let mut node_to_seq: HashMap<String, String> = HashMap::new();
        let flank_header_set: HashSet<String> =
            flank_rec_by_header.keys().cloned().collect();
        for (h, s) in &culled_records {
            if h.ends_with('.') || flank_header_set.contains(h) {
                continue;
            }
            let tag = tag_from_header(h);
            if tag != h {
                node_to_seq.insert(tag.to_string(), s.clone());
            }
        }

        let mut filtered_culled: Vec<(String, String)> = Vec::with_capacity(culled_records.len());
        let mut pending_flanks: Vec<(String, String, Recovery)> = Vec::new();
        let mut window_cache = WindowCache::new(aln_len);

        for (header, seq) in &culled_records {
            if !flank_rec_by_header.contains_key(header) {
                filtered_culled.push((header.clone(), seq.clone()));
                continue;
            }
            let mut rec = flank_rec_by_header.remove(header).unwrap();
            let anchor = rec.flank_anchor.clone().unwrap();
            let is_leading = anchor.is_leading;
            let adj_node = anchor.adj_node.clone();

            let Some(adj_seq) = node_to_seq.get(&adj_node) else { continue };
            let (adj_first_i, adj_last_i) = seq_residue_bounds(adj_seq);
            if adj_first_i < 0 {
                continue;
            }
            let adj_first = adj_first_i as usize;
            let adj_last = adj_last_i as usize;

            let (win_start, win_end) = if is_leading {
                (0usize, adj_first)
            } else {
                (adj_last + 1, aln_len)
            };
            if win_end <= win_start {
                continue;
            }

            let score = score_window(
                seq.as_bytes(),
                &ref_seqs_bytes,
                win_start, win_end,
                &mut window_cache,
            );
            if score.gap_residues < MIN_GAP_RESIDUES {
                continue;
            }

            let region = fmt_region(&rec.locus);
            let direction = if is_leading { "leading" } else { "trailing" };
            let suffix = format!("direction={} adj_node={}", direction, adj_node);

            if debug_level > 0 {
                rejection_log.push(format!(
                    "DEBUG_SCORE {} blosum={:.1}%/{:.1}% scored={}/{} gap_residues={} {} region={}",
                    rec.tag,
                    score.blosum_pct, score.covered_pct,
                    score.scored_cols, score.pssm_total_cols,
                    score.gap_residues, suffix, region,
                ));
            }

            if score.covered_pct < MIN_COVERED_PCT {
                rejection_log.push(format!(
                    "REJECTED {} reason=low_blosum blosum={:.1}%/{:.1}% thr={}% scored={}/{} gap_residues={} aa_len={} {} region={}",
                    rec.tag,
                    score.blosum_pct, score.covered_pct,
                    MIN_COVERED_PCT,
                    score.scored_cols, score.pssm_total_cols,
                    score.gap_residues,
                    count_residues(seq),
                    suffix, region,
                ));
                continue;
            }

            rec.locus.msa_cols = data_cols_set(seq);
            rec.score = ScoreInfo {
                full_blosum: score.blosum_pct,
                covered_blosum: score.covered_pct,
                gap_residues: score.gap_residues,
                win_start: score.win_start,
                win_end: score.win_end,
            };
            pending_flanks.push((header.clone(), seq.clone(), rec));
        }

        // Group by (adj_node, is_leading) and dedup.
        let mut groups: HashMap<(String, bool), Vec<(String, String, Recovery)>> = HashMap::new();
        for entry in pending_flanks {
            let key = (
                entry.2.flank_anchor.as_ref().unwrap().adj_node.clone(),
                entry.2.flank_anchor.as_ref().unwrap().is_leading,
            );
            groups.entry(key).or_default().push(entry);
        }

        let mut group_keys: Vec<(String, bool)> = groups.keys().cloned().collect();
        group_keys.sort();

        for key in group_keys {
            let mut cands = groups.remove(&key).unwrap();
            cands.sort_by(|a, b| {
                b.2.score.full_blosum.partial_cmp(&a.2.score.full_blosum)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let group_label = format!(
                "adj_node={} direction={}",
                key.0,
                if key.1 { "leading" } else { "trailing" }
            );
            let (accepted_idx, counts) = dedup_group(
                &mut cands,
                "flank",
                &group_label,
                |r| {
                    if r.flank_anchor.as_ref().map(|a| a.is_leading).unwrap_or(false) {
                        "N-terminal"
                    } else {
                        "C-terminal"
                    }
                },
                &mut rejection_log,
                debug_level,
            );
            for i in &accepted_idx {
                let (h, s, r) = &cands[*i];
                filtered_culled.push((h.clone(), s.clone()));
                surviving_flanks.push(r.clone());
            }
            emit_group_summary(&mut rejection_log, "flank", &group_label, &counts);
        }

        culled_records = filtered_culled;
    }
    let t_flank = t_flank_total.elapsed().as_secs_f64();

    // -------------------------------------------------------------------
    // 7. Gap scoring + per-cluster dedup + cross-cluster dedup + IMX late kick.
    // -------------------------------------------------------------------
    let t_gap_total = Instant::now();
    let mut surviving_gaps: Vec<Recovery> = Vec::new();
    let mut imx_late_kicks: HashSet<String> = HashSet::new();

    let mut gap_recoveries_by_header: HashMap<String, Vec<Recovery>> = HashMap::new();
    if !gap_inputs_by_header.is_empty() {
        let culled_header_set: HashSet<String> = culled_records
            .iter()
            .filter(|(h, _)| !h.ends_with('.'))
            .map(|(h, _)| h.clone())
            .collect();
        for header in culled_header_set.iter() {
            let Some(raw_entries) = gap_inputs_by_header.get(header) else { continue };
            for gi in raw_entries {
                let imx_anchor = if !gi.imx_slot_node.is_empty() {
                    Some(ImxAnchor {
                        slot_node: gi.imx_slot_node.clone(),
                        base_cluster_key: gi.cluster_key.clone(),
                        left_flank: gi.node_a_name.clone(),
                        right_flank: gi.node_b_name.clone(),
                    })
                } else {
                    None
                };
                let mut rec = Recovery {
                    tag: gi.tag.clone(),
                    header: gi.header.clone(),
                    aa_seq: gi.aa_seq.clone(),
                    nt_seq: gi.nt_seq.clone(),
                    locus: Locus {
                        scaffold: gi.scaffold.clone(),
                        strand: gi.strand.clone(),
                        bp_start: gi.hit_start,
                        bp_end: gi.hit_end,
                        msa_cols: BTreeSet::new(),
                        frame: gi.frame,
                    },
                    score: ScoreInfo {
                        full_blosum: 0.0,
                        covered_blosum: 0.0,
                        gap_residues: 0,
                        win_start: gi.gap_start,
                        win_end: gi.gap_end,
                    },
                    kind: RecoveryKind::Gap,
                    cluster_key: gi.cluster_key.clone(),
                    gene_key: gi.gene_key.clone(),
                    flank_anchor: None,
                    gap_anchor: Some(GapAnchor {
                        node_a: gi.node_a_name.clone(),
                        node_b: gi.node_b_name.clone(),
                    }),
                    imx_anchor,
                    module: None,
                    final_tag: None,
                };
                if let Some(drops) = nt_cull_map.get(header) {
                    apply_cull_trim(&mut rec, drops);
                }
                gap_recoveries_by_header.entry(header.clone()).or_default().push(rec);
            }
        }
    }

    if !gap_recoveries_by_header.is_empty() {
        let ref_seq_strings: Vec<String> = culled_records
            .iter()
            .filter(|(h, _)| h.ends_with('.'))
            .map(|(_, s)| s.clone())
            .collect();
        let ref_seqs_bytes: Vec<&[u8]> =
            ref_seq_strings.iter().map(|s| s.as_bytes()).collect();

        let mut gap_node_to_seq: HashMap<String, String> = HashMap::new();
        let gap_header_set: HashSet<String> =
            gap_recoveries_by_header.keys().cloned().collect();
        for (h, s) in &culled_records {
            if h.ends_with('.') || gap_header_set.contains(h) {
                continue;
            }
            let tag = tag_from_header(h);
            if tag != h {
                gap_node_to_seq.insert(tag.to_string(), s.clone());
            }
        }

        let mut filtered_culled: Vec<(String, String)> = Vec::with_capacity(culled_records.len());
        let mut pending_by_cluster: HashMap<String, Vec<(String, String, Recovery)>> = HashMap::new();
        let mut window_cache = WindowCache::new(aln_len);
        let mut bounds_cache: HashMap<String, (i64, i64)> = HashMap::new();
        let mut gap_counts: HashMap<&'static str, usize> = HashMap::new();

        for (header, seq) in &culled_records {
            let Some(recs) = gap_recoveries_by_header.remove(header) else {
                filtered_culled.push((header.clone(), seq.clone()));
                continue;
            };
            let mut best_rec: Option<Recovery> = None;
            let mut best_blosum: f64 = -1.0;
            let mut all_rejected: Vec<String> = Vec::new();

            for mut rec in recs {
                let node_a_name = rec.gap_anchor.as_ref().unwrap().node_a.clone();
                let node_b_name = rec.gap_anchor.as_ref().unwrap().node_b.clone();
                let seq_a = gap_node_to_seq.get(&node_a_name).cloned();
                let seq_b = gap_node_to_seq.get(&node_b_name).cloned();

                let strand = rec.locus.strand.clone();
                let cluster = rec.cluster_key.clone();
                let gap_start = rec.score.win_start;
                let gap_end = rec.score.win_end;

                let (win_start, win_end) = match (seq_a, seq_b) {
                    (Some(sa), Some(sb)) => {
                        let a_bounds = *bounds_cache
                            .entry(node_a_name.clone())
                            .or_insert_with(|| seq_residue_bounds(&sa));
                        let b_bounds = *bounds_cache
                            .entry(node_b_name.clone())
                            .or_insert_with(|| seq_residue_bounds(&sb));
                        let (a_first, a_last) = a_bounds;
                        let (b_first, b_last) = b_bounds;
                        if a_first == -1 || b_first == -1 {
                            (gap_start, gap_end.min(aln_len))
                        } else if strand == "-" {
                            ((b_last as usize) + 1, a_first as usize)
                        } else {
                            ((a_last as usize) + 1, b_first as usize)
                        }
                    }
                    _ => (gap_start, gap_end.min(aln_len)),
                };

                let score = score_window(
                    seq.as_bytes(),
                    &ref_seqs_bytes,
                    win_start, win_end,
                    &mut window_cache,
                );
                if score.gap_residues < MIN_GAP_RESIDUES {
                    *gap_counts.entry("kick_gap_residues").or_default() += 1;
                    continue;
                }
                let region = fmt_region(&rec.locus);
                let suffix = format!("cluster={}", cluster);
                if debug_level > 0 {
                    rejection_log.push(format!(
                        "DEBUG_SCORE {} blosum={:.1}%/{:.1}% scored={}/{} gap_residues={} {} region={}",
                        rec.tag,
                        score.blosum_pct, score.covered_pct,
                        score.scored_cols, score.pssm_total_cols,
                        score.gap_residues, suffix, region,
                    ));
                }
                if score.covered_pct < MIN_COVERED_PCT {
                    *gap_counts.entry("kick_low_blosum").or_default() += 1;
                    all_rejected.push(format!(
                        "REJECTED {} reason=low_blosum blosum={:.1}%/{:.1}% thr={}% scored={}/{} gap_residues={} aa_len={} {} region={}",
                        rec.tag,
                        score.blosum_pct, score.covered_pct,
                        MIN_COVERED_PCT,
                        score.scored_cols, score.pssm_total_cols,
                        score.gap_residues,
                        count_residues(seq),
                        suffix, region,
                    ));
                    continue;
                }
                if score.blosum_pct > best_blosum {
                    *gap_counts.entry("pass_blosum").or_default() += 1;
                    best_blosum = score.blosum_pct;
                    rec.locus.msa_cols = data_cols_set(seq);
                    rec.score = ScoreInfo {
                        full_blosum: score.blosum_pct,
                        covered_blosum: score.covered_pct,
                        gap_residues: score.gap_residues,
                        win_start: score.win_start,
                        win_end: score.win_end,
                    };
                    best_rec = Some(rec);
                }
            }

            if let Some(br) = best_rec {
                pending_by_cluster.entry(br.cluster_key.clone()).or_default().push(
                    (header.clone(), seq.clone(), br),
                );
            } else {
                rejection_log.extend(all_rejected);
            }
        }

        // Cross-cluster genomic dedup.
        cross_cluster_genomic_dedup(&mut pending_by_cluster, &mut rejection_log);

        let mut cluster_keys: Vec<String> = pending_by_cluster.keys().cloned().collect();
        cluster_keys.sort();
        for ck in cluster_keys {
            let mut cands = pending_by_cluster.remove(&ck).unwrap();
            cands.sort_by(|a, b| {
                b.2.score.full_blosum.partial_cmp(&a.2.score.full_blosum)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let group_label = format!("cluster={}", ck);
            let (accepted_idx, counts) = dedup_group(
                &mut cands,
                "gap",
                &group_label,
                |_r| "IMX",
                &mut rejection_log,
                debug_level,
            );
            for i in &accepted_idx {
                let (h, s, r) = &cands[*i];
                filtered_culled.push((h.clone(), s.clone()));
                surviving_gaps.push(r.clone());
            }
            emit_group_summary(&mut rejection_log, "gap", &group_label, &counts);
        }
        culled_records = filtered_culled;

        // IMX late kick.
        let seq_by_token: HashMap<String, String> = culled_records
            .iter()
            .filter(|(h, _)| !h.ends_with('.'))
            .filter_map(|(h, s)| {
                let t = tag_from_header(h);
                if t != h { Some((t.to_string(), s.clone())) } else { None }
            })
            .collect();
        let (kept, kicked_headers_imx, new_kicks) = imx_late_kick(
            surviving_gaps,
            &seq_by_token,
            &clusters_for_gene,
            &mut rejection_log,
            &mut gap_counts,
        );
        surviving_gaps = kept;
        imx_late_kicks.extend(new_kicks);
        if !kicked_headers_imx.is_empty() {
            culled_records.retain(|(h, _)| !kicked_headers_imx.contains(h));
        }
    }
    let t_gap = t_gap_total.elapsed().as_secs_f64();

    // -------------------------------------------------------------------
    // 8. Finalize: orig_aa_len, kicked sets, supersession.
    // -------------------------------------------------------------------
    let t_finalize_start = Instant::now();
    let mut cluster_by_header: HashMap<String, String> = HashMap::new();
    let mut orig_aa_len_by_header: HashMap<String, usize> = HashMap::new();
    for fr in &surviving_flanks {
        if !fr.cluster_key.is_empty() {
            cluster_by_header.insert(fr.header.clone(), fr.cluster_key.clone());
            orig_aa_len_by_header.insert(fr.header.clone(), fr.aa_seq.len());
        }
    }
    for gr in &surviving_gaps {
        if !gr.cluster_key.is_empty() {
            cluster_by_header.insert(gr.header.clone(), gr.cluster_key.clone());
            orig_aa_len_by_header.insert(gr.header.clone(), gr.aa_seq.len());
        }
    }
    for (header, nt_seq) in cull_nt_seqs.iter() {
        if header.ends_with('.') {
            continue;
        }
        orig_aa_len_by_header
            .entry(header.clone())
            .or_insert(nt_seq.len() / 3);
    }

    let mut kicked_headers: HashSet<String> = HashSet::new();
    let mut kicked_nodes: HashSet<String> = HashSet::new();
    kicked_nodes.extend(imx_late_kicks.iter().cloned());

    let culled_header_set: HashSet<String> = culled_records.iter().map(|(h, _)| h.clone()).collect();
    // Natives + flanks/gaps that disappeared during cull/scoring.
    let existing_headers: HashSet<String> =
        natives_aa.iter().map(|(h, _)| h.clone()).collect();
    for h in &existing_headers {
        if !culled_header_set.contains(h) {
            let node_name = tag_from_header(h);
            if !node_name.is_empty() {
                kicked_nodes.insert(node_name.to_string());
            }
        }
    }

    // -------------------------------------------------------------------
    // 9. Assign final tags + emit gff/recovered_log lines.
    // -------------------------------------------------------------------
    let mut flank_renames: HashMap<String, String> = HashMap::new();
    let mut gap_renames: HashMap<String, String> = HashMap::new();
    assign_tags(
        &mut surviving_flanks,
        RecoveryKind::Flank,
        &gene_key,
        &mut flank_renames,
        &mut worker_gff_lines,
        &mut recovered_log,
    );
    assign_tags(
        &mut surviving_gaps,
        RecoveryKind::Gap,
        &gene_key,
        &mut gap_renames,
        &mut worker_gff_lines,
        &mut recovered_log,
    );
    let mut rename_map: HashMap<String, String> = HashMap::new();
    rename_map.extend(flank_renames);
    rename_map.extend(gap_renames);

    // Mutable: reorder refs first then candidates, uppercase, apply renames.
    let mut refs_out: Vec<(String, String)> = Vec::new();
    let mut cands_out: Vec<(String, String)> = Vec::new();
    for (h, s) in &culled_records {
        let upper = s.to_ascii_uppercase();
        if h.ends_with('.') {
            refs_out.push((h.clone(), upper));
        } else {
            cands_out.push((h.clone(), upper));
        }
    }
    culled_records.clear();
    culled_records.extend(refs_out);
    culled_records.extend(cands_out);

    let renamed_records: Vec<(String, String)>;
    let renamed_cluster_by_header: HashMap<String, String>;
    let renamed_orig_aa_len: HashMap<String, usize>;
    let mut renamed_nt_cull_map: HashMap<String, HashSet<usize>>;
    let renamed_cull_nt_seqs: HashMap<String, String>;
    if rename_map.is_empty() {
        renamed_records = culled_records.clone();
        renamed_cluster_by_header = cluster_by_header.clone();
        renamed_orig_aa_len = orig_aa_len_by_header.clone();
        renamed_nt_cull_map = nt_cull_map.clone();
        renamed_cull_nt_seqs = cull_nt_seqs.clone();
    } else {
        let rh = |hdr: &str| -> String {
            let parts: Vec<&str> = hdr.split('|').collect();
            if parts.len() <= 3 {
                return hdr.to_string();
            }
            if let Some(new_tag) = rename_map.get(parts[3]) {
                let mut v: Vec<String> = parts.iter().map(|s| s.to_string()).collect();
                v[3] = new_tag.clone();
                v.join("|")
            } else {
                hdr.to_string()
            }
        };
        renamed_records = culled_records.iter().map(|(h, s)| (rh(h), s.clone())).collect();
        renamed_cluster_by_header = cluster_by_header
            .iter()
            .map(|(h, v)| (rh(h), v.clone()))
            .collect();
        renamed_orig_aa_len = orig_aa_len_by_header
            .iter()
            .map(|(h, v)| (rh(h), *v))
            .collect();
        renamed_nt_cull_map = nt_cull_map.iter().map(|(h, v)| (rh(h), v.clone())).collect();
        renamed_cull_nt_seqs = cull_nt_seqs
            .iter()
            .map(|(h, v)| (rh(h), v.clone()))
            .collect();
    }

    // Supersession (pre-split, IMX only) + plain-gap refind rejection.
    let (
        superseded_names,
        supersession_pairs,
        refind_rejects,
        supersession_audit,
    ) = supersede_natives(&surviving_gaps, &gff_nodes, &clusters_for_gene);

    // Drop rejected re-find gaps from the surviving set and worker GFF.
    let mut renamed_records = renamed_records;
    if !refind_rejects.is_empty() {
        surviving_gaps = surviving_gaps
            .into_iter()
            .filter(|r| !refind_rejects.contains(r.current_tag()))
            .collect();
        worker_gff_lines.retain(|ln| {
            !refind_rejects
                .iter()
                .any(|t| ln.contains(&format!("Name={};", t)))
        });
    }
    let drop_tags: HashSet<String> = superseded_names
        .iter()
        .cloned()
        .chain(refind_rejects.iter().cloned())
        .collect();
    let mut renamed_cluster_by_header = renamed_cluster_by_header;
    if !drop_tags.is_empty() {
        renamed_records = renamed_records
            .into_iter()
            .filter(|(h, _)| {
                if h.ends_with('.') {
                    return true;
                }
                let t = tag_from_header(h);
                !drop_tags.contains(t)
            })
            .collect();
        renamed_cluster_by_header.retain(|h, _| {
            let t = tag_from_header(h);
            !drop_tags.contains(t)
        });
        kicked_nodes.extend(superseded_names.iter().cloned());
    }

    // -------------------------------------------------------------------
    // 10. Build native_cluster_members for run_final_split.
    // -------------------------------------------------------------------
    let mut native_cluster_members: HashMap<String, Vec<String>> = HashMap::new();
    if !renamed_cluster_by_header.is_empty() {
        let mut token_to_header: HashMap<String, String> = HashMap::new();
        for (h, _) in &renamed_records {
            if h.ends_with('.') {
                continue;
            }
            let t = tag_from_header(h);
            if t != h.as_str() {
                token_to_header.insert(t.to_string(), h.clone());
            }
        }
        let recovery_cks: HashSet<String> = renamed_cluster_by_header.values().cloned().collect();
        for ck in &recovery_cks {
            let mut seen_h: HashSet<String> = HashSet::new();
            if let Some(toks) = clusters_for_gene.get(ck) {
                let entry = native_cluster_members.entry(ck.clone()).or_default();
                for tok in toks {
                    if let Some(hdr) = token_to_header.get(tok) {
                        if !seen_h.contains(hdr) {
                            entry.push(hdr.clone());
                            seen_h.insert(hdr.clone());
                        }
                    }
                }
            }
        }
    }

    // -------------------------------------------------------------------
    // 11. Extend gff_nodes with final-tag entries so run_final_split's
    //     header_to_gff lookup resolves renamed recovery headers.
    // -------------------------------------------------------------------
    let mut gff_nodes_extended = gff_nodes.clone();
    for r in surviving_flanks.iter().chain(surviving_gaps.iter()) {
        if let Some(ft) = &r.final_tag {
            gff_nodes_extended.insert(
                ft.clone(),
                GffEntryInput {
                    scaffold: r.locus.scaffold.clone(),
                    start: r.locus.bp_start,
                    end: r.locus.bp_end,
                    strand: r.locus.strand.clone(),
                },
            );
        }
    }

    // -------------------------------------------------------------------
    // 12. run_final_split.
    // -------------------------------------------------------------------
    let (post_split_records, split_kicked_headers, gene_split_trim, gene_split_calc_log) =
        run_final_split(
            &gene_key,
            renamed_records,
            &renamed_cluster_by_header,
            &renamed_orig_aa_len,
            &mut renamed_nt_cull_map,
            &gff_nodes_extended,
            skip_exon_split,
            &native_cluster_members,
            debug_level >= 2,
        );

    // Apply split trims to Recovery loci.  Build an index map up-front, then
    // mutate by index so we never need raw pointers or unsafe aliasing.
    if !gene_split_trim.is_empty() {
        let mut tag_to_idx: HashMap<String, (bool, usize)> = HashMap::new();
        for (i, r) in surviving_flanks.iter().enumerate() {
            tag_to_idx.insert(r.current_tag().to_string(), (true, i));
        }
        for (i, r) in surviving_gaps.iter().enumerate() {
            tag_to_idx.insert(r.current_tag().to_string(), (false, i));
        }
        for (node_name, (l, r_)) in &gene_split_trim {
            if *l == 0 && *r_ == 0 {
                continue;
            }
            if let Some(&(is_flank, i)) = tag_to_idx.get(node_name) {
                if is_flank {
                    let rec = &mut surviving_flanks[i];
                    rec.locus = rec.locus.trim_aa(*l, *r_);
                } else {
                    let rec = &mut surviving_gaps[i];
                    rec.locus = rec.locus.trim_aa(*l, *r_);
                }
            }
        }
    }

    let mut split_kicked_tags: HashSet<String> = HashSet::new();
    for hdr in &split_kicked_headers {
        let t = tag_from_header(hdr);
        if !t.is_empty() {
            split_kicked_tags.insert(t.to_string());
        }
    }
    kicked_headers.extend(split_kicked_headers.iter().cloned());
    kicked_nodes.extend(split_kicked_tags.iter().cloned());

    // -------------------------------------------------------------------
    // 12. Sort cands + build NT records (codon-dropped per nt_cull_map).
    // -------------------------------------------------------------------
    let mut refs_final: Vec<(String, String)> = Vec::new();
    let mut cands_final: Vec<(String, String)> = Vec::new();
    for (h, s) in &post_split_records {
        if h.ends_with('.') {
            refs_final.push((h.clone(), s.clone()));
        } else {
            cands_final.push((h.clone(), s.clone()));
        }
    }
    cands_final.sort_by(|a, b| {
        let ka = genomic_sort_key(&a.0, &gff_nodes_extended);
        let kb = genomic_sort_key(&b.0, &gff_nodes_extended);
        ka.cmp(&kb)
    });
    let mut final_records: Vec<(String, String)> = Vec::with_capacity(refs_final.len() + cands_final.len());
    final_records.extend(refs_final);
    final_records.extend(cands_final);

    // NT records: codon-drop per renamed_nt_cull_map.
    let t_nt = Instant::now();
    let mut nt_records: Vec<(String, String)> = Vec::new();
    for (header, _) in &final_records {
        if header.ends_with('.') {
            continue;
        }
        let Some(nt_seq) = renamed_cull_nt_seqs.get(header) else { continue };
        let drops = renamed_nt_cull_map.get(header);
        let nt_out = if let Some(drops) = drops {
            let mut drop_sorted: Vec<usize> = drops.iter().copied().collect();
            drop_sorted.sort_unstable_by(|a, b| b.cmp(a));
            let mut codons: Vec<String> = (0..nt_seq.len()).step_by(3)
                .map(|i| nt_seq[i..(i + 3).min(nt_seq.len())].to_string())
                .collect();
            for idx in drop_sorted {
                if idx < codons.len() {
                    codons.remove(idx);
                }
            }
            codons.concat().to_ascii_uppercase()
        } else {
            nt_seq.to_ascii_uppercase()
        };
        nt_records.push((header.clone(), nt_out));
    }
    let t_nt = t_nt.elapsed().as_secs_f64();

    // -------------------------------------------------------------------
    // 13. Build node_occupancy (only when gene has gaps).
    // -------------------------------------------------------------------
    let gene_has_gaps = final_records
        .iter()
        .any(|(h, _)| !h.ends_with('.') && tag_from_header(h).starts_with("GAP_"));
    let mut node_occupancy: HashMap<String, BTreeSet<usize>> = HashMap::new();
    if gene_has_gaps {
        for (h, s) in &final_records {
            if h.ends_with('.') {
                continue;
            }
            let t = tag_from_header(h);
            if t != h.as_str() {
                node_occupancy.insert(t.to_string(), data_cols_set(s));
            }
        }
    }

    let t_finalize = t_finalize_start.elapsed().as_secs_f64();
    let _ = t_total;

    // -------------------------------------------------------------------
    // 14. Build Python return dict.
    // -------------------------------------------------------------------
    build_result_dict(
        py,
        &gene_key,
        kicked_headers,
        kicked_nodes,
        surviving_flanks,
        surviving_gaps,
        rejection_log,
        (t_hmm, t_cull, t_nt, t_flank, t_gap, t_finalize),
        gene_split_trim,
        gene_split_calc_log,
        split_kicked_tags,
        worker_gff_lines,
        recovered_log,
        node_occupancy,
        supersession_pairs,
        supersession_audit,
        final_records,
        nt_records,
    )
}

fn assign_tags(
    survivors: &mut Vec<Recovery>,
    kind: RecoveryKind,
    gene_key: &str,
    renames: &mut HashMap<String, String>,
    worker_gff_lines: &mut Vec<String>,
    recovered_log: &mut Vec<String>,
) {
    survivors.sort_by(|a, b| {
        (
            &a.locus.scaffold,
            &a.locus.strand,
            a.locus.bp_start,
            a.locus.bp_end,
            a.locus.frame,
        )
            .cmp(&(
                &b.locus.scaffold,
                &b.locus.strand,
                b.locus.bp_start,
                b.locus.bp_end,
                b.locus.frame,
            ))
    });
    let gff_source = if kind == RecoveryKind::Flank { "FlankScan" } else { "ExonFinder" };
    let note_prefix = if kind == RecoveryKind::Flank { "flank_recovered" } else { "gap_recovered" };

    for (idx, r) in survivors.iter_mut().enumerate() {
        let tmp_tag = r.tag.clone();
        let gene_base = r.gene_key.split('.').next().unwrap_or(&r.gene_key).to_string();
        let final_tag = format!("{}_{}_{}", kind.as_str(), gene_base, idx);
        r.final_tag = Some(final_tag.clone());
        renames.insert(tmp_tag, final_tag.clone());

        let culled_aa_len = r.locus.bp_end.saturating_sub(r.locus.bp_start) / 3;
        let (attrs, anchor_line, header_suffix) = match kind {
            RecoveryKind::Flank => {
                let fa = r.flank_anchor.as_ref().unwrap();
                let direction = if fa.is_leading { "leading" } else { "trailing" };
                (
                    format!(
                        "ID={};Name={};Parent={};Note={},frame={},aa_len={}",
                        final_tag, final_tag, gene_base,
                        note_prefix, r.locus.frame, culled_aa_len,
                    ),
                    format!("  adj_node={}", fa.adj_node),
                    format!(" ({})", direction),
                )
            }
            RecoveryKind::Gap => {
                let ga = r.gap_anchor.as_ref().unwrap();
                (
                    format!(
                        "ID={};Name={};Parent={};Note={},frame={},aa_len={},between={}+{},cluster={}",
                        final_tag, final_tag, gene_base,
                        note_prefix, r.locus.frame, culled_aa_len,
                        ga.node_a, ga.node_b, r.cluster_key,
                    ),
                    format!("  between={}+{}", ga.node_a, ga.node_b),
                    String::new(),
                )
            }
        };

        worker_gff_lines.push(format!(
            "{}\t{}\texon\t{}\t{}\t.\t{}\t.\t{}",
            r.locus.scaffold, gff_source,
            r.locus.bp_start, r.locus.bp_end,
            r.locus.strand, attrs,
        ));

        let region = fmt_region(&r.locus);
        let module_label = match r.module.as_ref() {
            Some(c) => format!(" MODULE={}", c.module_type),
            None => String::new(),
        };
        recovered_log.push(format!(
            "RECOVERED {}{} gene={} cluster={} win_cols={}-{} frame={}{}\n  region={}\n{}\n  full_blosum={:.1}% covered_blosum={:.1}% gap_residues={}\n  aa={}\n  nt={}",
            final_tag, header_suffix, gene_key, r.cluster_key,
            r.score.win_start, r.score.win_end, r.locus.frame, module_label,
            region, anchor_line,
            r.score.full_blosum, r.score.covered_blosum, r.score.gap_residues,
            r.aa_seq, r.nt_seq,
        ));
    }
}

// ===========================================================================
// Python dict construction
// ===========================================================================

fn empty_result(
    py: Python<'_>,
    gene_key: &str,
    refs_aa: &[(String, String)],
    natives_aa: &[(String, String)],
) -> PyResult<Py<PyDict>> {
    let mut records: Vec<(String, String)> =
        Vec::with_capacity(refs_aa.len() + natives_aa.len());
    records.extend(refs_aa.iter().cloned());
    records.extend(natives_aa.iter().cloned());
    build_result_dict(
        py,
        gene_key,
        HashSet::new(),
        HashSet::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        HashMap::new(),
        Vec::new(),
        HashSet::new(),
        Vec::new(),
        Vec::new(),
        HashMap::new(),
        HashMap::new(),
        Vec::new(),
        records,
        Vec::new(),
    )
}

fn recovery_to_dict(py: Python<'_>, r: &Recovery) -> PyResult<Py<PyDict>> {
    let d = PyDict::new(py);
    d.set_item("tag", &r.tag)?;
    d.set_item("header", &r.header)?;
    d.set_item("aa_seq", &r.aa_seq)?;
    d.set_item("nt_seq", &r.nt_seq)?;

    let locus = PyDict::new(py);
    locus.set_item("scaffold", &r.locus.scaffold)?;
    locus.set_item("strand", &r.locus.strand)?;
    locus.set_item("bp_start", r.locus.bp_start)?;
    locus.set_item("bp_end", r.locus.bp_end)?;
    let cols_list = PyList::new(py, r.locus.msa_cols.iter().copied())?;
    locus.set_item("msa_cols", cols_list)?;
    locus.set_item("frame", r.locus.frame)?;
    d.set_item("locus", locus)?;

    let score = PyDict::new(py);
    score.set_item("full_blosum", r.score.full_blosum)?;
    score.set_item("covered_blosum", r.score.covered_blosum)?;
    score.set_item("gap_residues", r.score.gap_residues)?;
    score.set_item("win_start", r.score.win_start)?;
    score.set_item("win_end", r.score.win_end)?;
    d.set_item("score", score)?;

    d.set_item("kind", r.kind.as_str())?;
    d.set_item("cluster_key", &r.cluster_key)?;
    d.set_item("gene_key", &r.gene_key)?;

    if let Some(fa) = &r.flank_anchor {
        let fd = PyDict::new(py);
        fd.set_item("adj_node", &fa.adj_node)?;
        fd.set_item("is_leading", fa.is_leading)?;
        d.set_item("flank_anchor", fd)?;
    } else {
        d.set_item("flank_anchor", py.None())?;
    }
    if let Some(ga) = &r.gap_anchor {
        let gd = PyDict::new(py);
        gd.set_item("node_a", &ga.node_a)?;
        gd.set_item("node_b", &ga.node_b)?;
        d.set_item("gap_anchor", gd)?;
    } else {
        d.set_item("gap_anchor", py.None())?;
    }
    if let Some(mx) = &r.imx_anchor {
        let md = PyDict::new(py);
        md.set_item("slot_node", &mx.slot_node)?;
        md.set_item("base_cluster_key", &mx.base_cluster_key)?;
        md.set_item("left_flank", &mx.left_flank)?;
        md.set_item("right_flank", &mx.right_flank)?;
        d.set_item("imx_anchor", md)?;
    } else {
        d.set_item("imx_anchor", py.None())?;
    }
    if let Some(c) = &r.module {
        let cd = PyDict::new(py);
        cd.set_item("partner_tag", &c.partner_tag)?;
        cd.set_item("module_type", &c.module_type)?;
        let cols = PyList::new(py, c.slot_cols.iter().copied())?;
        cd.set_item("slot_cols", cols)?;
        d.set_item("module", cd)?;
    } else {
        d.set_item("module", py.None())?;
    }
    match &r.final_tag {
        Some(ft) => d.set_item("final_tag", ft)?,
        None => d.set_item("final_tag", py.None())?,
    }
    Ok(d.into())
}

#[allow(clippy::too_many_arguments)]
fn build_result_dict(
    py: Python<'_>,
    gene_key: &str,
    kicked_headers: HashSet<String>,
    kicked_nodes: HashSet<String>,
    surviving_flanks: Vec<Recovery>,
    surviving_gaps: Vec<Recovery>,
    rejection_log: Vec<String>,
    timings: (f64, f64, f64, f64, f64, f64),
    split_trim: HashMap<String, (usize, usize)>,
    split_calc_log: Vec<String>,
    split_kicked_tags: HashSet<String>,
    gff_lines: Vec<String>,
    recovered_log: Vec<String>,
    node_occupancy: HashMap<String, BTreeSet<usize>>,
    supersession_pairs: HashMap<String, String>,
    supersession_audit: Vec<String>,
    post_split_aa_records: Vec<(String, String)>,
    post_split_nt_records: Vec<(String, String)>,
) -> PyResult<Py<PyDict>> {
    let out = PyDict::new(py);

    out.set_item("gene_key", gene_key)?;
    out.set_item("kicked_headers", PyList::new(py, kicked_headers.iter())?)?;
    out.set_item("kicked_nodes", PyList::new(py, kicked_nodes.iter())?)?;

    let flanks_list = PyList::empty(py);
    for r in &surviving_flanks {
        flanks_list.append(recovery_to_dict(py, r)?)?;
    }
    out.set_item("surviving_flanks", flanks_list)?;
    let gaps_list = PyList::empty(py);
    for r in &surviving_gaps {
        gaps_list.append(recovery_to_dict(py, r)?)?;
    }
    out.set_item("surviving_gaps", gaps_list)?;

    out.set_item("rejection_log", PyList::new(py, rejection_log.iter())?)?;
    out.set_item("recovered_log", PyList::new(py, recovered_log.iter())?)?;
    out.set_item("gff_lines", PyList::new(py, gff_lines.iter())?)?;
    out.set_item("split_calc_log", PyList::new(py, split_calc_log.iter())?)?;

    out.set_item(
        "timings",
        (timings.0, timings.1, timings.2, timings.3, timings.4, timings.5),
    )?;

    let trim = PyDict::new(py);
    for (k, (l, r_)) in &split_trim {
        trim.set_item(k, (*l, *r_))?;
    }
    out.set_item("split_trim", trim)?;

    out.set_item("split_kicked_tags", PyList::new(py, split_kicked_tags.iter())?)?;

    let occ = PyDict::new(py);
    for (k, cols) in &node_occupancy {
        occ.set_item(k, PyList::new(py, cols.iter().copied())?)?;
    }
    out.set_item("node_occupancy", occ)?;

    let sp = PyDict::new(py);
    for (k, v) in &supersession_pairs {
        sp.set_item(k, v)?;
    }
    out.set_item("supersession_pairs", sp)?;
    out.set_item("supersession_audit", PyList::new(py, supersession_audit.iter())?)?;

    out.set_item("post_split_aa_records", post_split_aa_records)?;
    out.set_item("post_split_nt_records", post_split_nt_records)?;

    Ok(out.into())
}
