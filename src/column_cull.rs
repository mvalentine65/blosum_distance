use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

// ---------------------------------------------------------------------------
// Constants for retained-intron detection
// ---------------------------------------------------------------------------
const MIN_INTRON_NT: usize = 30;

// ---------------------------------------------------------------------------
// Constants for stop-codon BLOSUM trim
// ---------------------------------------------------------------------------

/// Maximum distance (non-gap AA) from a data edge for a stop to be
/// eligible for BLOSUM fragment scoring.  Stops deeper than this on
/// both sides are internal pseudogenic stops -- skip them.
const STOP_TRIM_MAX_DEPTH_AA: usize = 40;

/// The smaller fragment must score below this to be considered junk.
const STOP_TRIM_BAD_THRESHOLD: f64 = 0.50;

/// The larger fragment must score above this to confirm it is real.
const STOP_TRIM_GOOD_THRESHOLD: f64 = 0.60;

/// Minimum fraction of refs that must have data at a column for it
/// to enter the PSSM (matches exonfinder.py _compute_ref_data_cols).
const STOP_TRIM_MIN_REF_OCC: f64 = 0.30;

/// BLOSUM62 scoring matrix (upper-case AA only).
/// Returns score for (query, reference) pair.
fn blosum62(a: u8, b: u8) -> f64 {
    let idx = |c: u8| -> usize {
        match c {
            b'A' => 0, b'R' => 1, b'N' => 2, b'D' => 3, b'C' => 4,
            b'Q' => 5, b'E' => 6, b'G' => 7, b'H' => 8, b'I' => 9,
            b'L' => 10, b'K' => 11, b'M' => 12, b'F' => 13, b'P' => 14,
            b'S' => 15, b'T' => 16, b'W' => 17, b'Y' => 18, b'V' => 19,
            _ => 20,
        }
    };
    #[rustfmt::skip]
    const MAT: [[i8; 20]; 20] = [
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
    let ai = idx(a);
    let bi = idx(b);
    if ai >= 20 || bi >= 20 {
        return -4.0;
    }
    MAT[ai][bi] as f64
}

/// Build a BLOSUM62 PSSM from reference sequences.
///
/// For each column in `data_cols`, computes a weighted score for every
/// possible query amino acid based on the reference residue frequencies.
/// Returns {col: (scores_by_aa, max_possible_score)}.
fn build_blosum_pssm(
    ref_seqs: &[&[u8]],
    n_refs: usize,
    aln_len: usize,
) -> HashMap<usize, (HashMap<u8, f64>, f64)> {
    let mut pssm: HashMap<usize, (HashMap<u8, f64>, f64)> = HashMap::new();
    let aas: [u8; 20] = *b"ARNDCQEGHILKMFPSTWYV";

    for col in 0..aln_len {
        // 30% ref occupancy gate.
        let mut counts: HashMap<u8, f64> = HashMap::new();
        let mut total = 0usize;
        for rseq in ref_seqs {
            if col >= rseq.len() { continue; }
            let c = rseq[col];
            if c == b'-' || c == b'*' || c == b'.' { continue; }
            let cu = c.to_ascii_uppercase();
            *counts.entry(cu).or_default() += 1.0;
            total += 1;
        }
        if total == 0 { continue; }
        if (total as f64 / n_refs as f64) < STOP_TRIM_MIN_REF_OCC { continue; }

        let total_f = total as f64;
        let mut scores: HashMap<u8, f64> = HashMap::new();
        for &qa in &aas {
            let mut s = 0.0f64;
            for (&ref_aa, &cnt) in &counts {
                s += blosum62(qa, ref_aa) * cnt;
            }
            scores.insert(qa, s / total_f);
        }
        let col_max = scores.values().cloned().fold(f64::NEG_INFINITY, f64::max);
        if col_max > 0.0 {
            pssm.insert(col, (scores, col_max));
        }
    }
    pssm
}

/// Score a fragment of a sequence against a precomputed BLOSUM PSSM.
///
/// Returns (ratio, scored_cols) where ratio = candidate_total / max_total.
fn score_fragment(
    seq: &[u8],
    pssm: &HashMap<usize, (HashMap<u8, f64>, f64)>,
    start: usize,
    end: usize,
) -> (f64, usize) {
    let mut cand_total = 0.0f64;
    let mut max_total = 0.0f64;
    let mut scored = 0usize;
    for i in start..end {
        if i >= seq.len() { break; }
        let c = seq[i];
        if c == b'-' || c == b'*' { continue; }
        if let Some((scores, col_max)) = pssm.get(&i) {
            cand_total += scores.get(&c.to_ascii_uppercase()).copied().unwrap_or(-4.0);
            max_total += col_max;
            scored += 1;
        }
    }
    if max_total > 0.0 {
        (cand_total / max_total, scored)
    } else {
        (0.0, 0)
    }
}

/// Trim the smaller fragment around a stop codon if BLOSUM PSSM scoring
/// shows it is junk while the larger fragment is real.
///
/// Operates on sequences in post-gap-cull column space.  Returns a new
/// alignment with trimmed fragments masked to gap.
fn mask_stop_blosum(
    records: &[(String, Vec<u8>)],
    ref_suffix: &str,
    debug: i32,
    log_dir: &Option<String>,
) -> (Vec<(String, Vec<u8>)>, HashMap<String, HashSet<usize>>) {
    let aln_len = if records.is_empty() { 0 } else { records[0].1.len() };

    let ref_seqs: Vec<&[u8]> = records
        .iter()
        .filter(|(h, _)| h.ends_with(ref_suffix))
        .map(|(_, s)| s.as_slice())
        .collect();

    let n_refs = ref_seqs.len();
    let mut masked_cols: HashMap<String, HashSet<usize>> = HashMap::new();

    if n_refs == 0 || aln_len == 0 {
        return (records.to_vec(), masked_cols);
    }

    // Build PSSM once for all sequences.
    let pssm = build_blosum_pssm(&ref_seqs, n_refs, aln_len);

    let mut stop_log: Vec<(String, usize, usize, usize, f64, f64, String)> = Vec::new();
    let mut result = Vec::with_capacity(records.len());

    for (header, seq) in records {
        if header.ends_with(ref_suffix) {
            result.push((header.clone(), seq.clone()));
            continue;
        }

        let stops: Vec<usize> = seq.iter().enumerate()
            .filter(|(_, &c)| c == b'*')
            .map(|(i, _)| i)
            .collect();

        if stops.is_empty() {
            result.push((header.clone(), seq.clone()));
            continue;
        }

        let data_start = match seq.iter().position(|&c| c != b'-') {
            Some(p) => p,
            None => {
                result.push((header.clone(), seq.clone()));
                continue;
            }
        };
        let data_end = seq.iter().rposition(|&c| c != b'-').unwrap() + 1;

        let mut cols_to_mask: HashSet<usize> = HashSet::new();

        for &star in &stops {
            // Depth gate.
            let aa_left = (data_start..star)
                .filter(|&i| seq[i] != b'-' && seq[i] != b'*')
                .count();
            let aa_right = (star + 1..data_end)
                .filter(|&i| seq[i] != b'-' && seq[i] != b'*')
                .count();

            if aa_left > STOP_TRIM_MAX_DEPTH_AA
                && aa_right > STOP_TRIM_MAX_DEPTH_AA
            {
                continue;
            }

            // Score both fragments.
            let (left_score, left_scored) = score_fragment(seq, &pssm, data_start, star);
            let (right_score, right_scored) = score_fragment(seq, &pssm, star + 1, data_end);

            // Determine smaller/larger.
            let (smaller_score, larger_score, trim_side) = if aa_right <= aa_left {
                (right_score, left_score, "right")
            } else {
                (left_score, right_score, "left")
            };

            // Apply thresholds.
            if smaller_score >= STOP_TRIM_BAD_THRESHOLD {
                continue;
            }
            if larger_score <= STOP_TRIM_GOOD_THRESHOLD {
                continue;
            }

            // Mask the smaller fragment + the stop itself.
            if trim_side == "right" {
                cols_to_mask.extend(star..data_end);
            } else {
                cols_to_mask.extend(data_start..=star);
            }

            if debug >= 1 {
                stop_log.push((
                    header.clone(),
                    star,
                    aa_left,
                    aa_right,
                    left_score,
                    right_score,
                    trim_side.to_string(),
                ));
            }
        }

        if cols_to_mask.is_empty() {
            result.push((header.clone(), seq.clone()));
        } else {
            let mut masked = seq.clone();
            for &i in &cols_to_mask {
                if i < masked.len() {
                    masked[i] = b'-';
                }
            }
            masked_cols.insert(header.clone(), cols_to_mask);
            result.push((header.clone(), masked));
        }
    }

    // Debug log.
    if !stop_log.is_empty() && debug >= 1 {
        if let Some(ld) = log_dir {
            let log_path = Path::new(ld).join("stop_blosum_trim_debug.csv");
            if let Ok(mut f) = File::create(&log_path) {
                let _ = writeln!(
                    f,
                    "header,stop_col,aa_left,aa_right,left_score,right_score,trimmed_side"
                );
                for (h, stop_col, aa_l, aa_r, ls, rs, side) in &stop_log {
                    let safe = if h.contains(',') {
                        format!("\"{}\"", h)
                    } else {
                        h.clone()
                    };
                    let _ = writeln!(
                        f,
                        "{},{},{},{},{:.3},{:.3},{}",
                        safe, stop_col, aa_l, aa_r, ls, rs, side
                    );
                }
            }
        }
    }

    (result, masked_cols)
}

fn is_donor(a: u8, b: u8) -> bool {
    (a == b'G' && b == b'T') || (a == b'G' && b == b'C')
}

fn is_acceptor(a: u8, b: u8) -> bool {
    a == b'A' && b == b'G'
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Single-pass normalize: replace '.' with '-', strip ' ', uppercase. One allocation.
fn normalize_records(records: &[(String, String)]) -> Vec<(String, Vec<u8>)> {
    records
        .iter()
        .map(|(h, s)| {
            let mut norm = Vec::with_capacity(s.len());
            for &b in s.as_bytes() {
                match b {
                    b'.' => norm.push(b'-'),
                    b' ' => {}
                    _ => norm.push(b.to_ascii_uppercase()),
                }
            }
            (h.clone(), norm)
        })
        .collect()
}

/// Apply a boolean column mask in-place: set masked positions to b'-'.
fn mask_columns_bool(seq: &mut [u8], mask: &[bool]) {
    for (i, m) in mask.iter().enumerate() {
        if *m && i < seq.len() {
            seq[i] = b'-';
        }
    }
}

fn has_splice_sites(
    nt_seq: &str,
    block_first_residue: usize,
    block_residue_count: usize,
) -> bool {
    let base_start = block_first_residue * 3;
    let base_end = (block_first_residue + block_residue_count) * 3;

    let intron_len = base_end - base_start;
    if intron_len < MIN_INTRON_NT {
        return false;
    }

    let nt = nt_seq.as_bytes();
    for phase in 0..3usize {
        let d_start = base_start + phase;
        let a_end = base_end + phase;
        if d_start + 1 >= nt.len() || a_end > nt.len() {
            continue;
        }
        let da = nt[d_start].to_ascii_uppercase();
        let db = nt[d_start + 1].to_ascii_uppercase();
        let aa = nt[a_end - 2].to_ascii_uppercase();
        let ab = nt[a_end - 1].to_ascii_uppercase();
        if is_donor(da, db) && is_acceptor(aa, ab) {
            return true;
        }
    }
    false
}

/// Detect and mask retained same-frame introns in candidate sequences.
fn mask_retained_introns(
    normalized: &[(String, Vec<u8>)],
    ref_suffix: &str,
    nt_seqs: &Option<HashMap<String, String>>,
    min_ref_supported_gap: usize,
    debug: i32,
    log_dir: &Option<String>,
) -> (Vec<(String, Vec<u8>)>, HashMap<String, HashSet<usize>>) {
    let min_intron_aa: usize = 10;
    let ref_gap_threshold: f64 = 0.95;
    let max_bridge_residues: usize = 8;

    if normalized.is_empty() {
        return (Vec::new(), HashMap::new());
    }

    let aln_len = normalized[0].1.len();
    let ref_seqs: Vec<&[u8]> = normalized
        .iter()
        .filter(|(h, _)| h.ends_with(ref_suffix))
        .map(|(_, s)| s.as_slice())
        .collect();

    if ref_seqs.is_empty() {
        return (normalized.to_vec(), HashMap::new());
    }

    let n_refs = ref_seqs.len();

    // Pre-compute reference-gap columns as Vec<bool>.
    let mut ref_gap_col = vec![false; aln_len];
    for i in 0..aln_len {
        let gap_count = ref_seqs.iter().filter(|s| s[i] == b'-').count();
        ref_gap_col[i] = (gap_count as f64 / n_refs as f64) >= ref_gap_threshold;
    }

    let mut intron_log: Vec<(String, usize, usize, usize, String, String, usize, usize, usize)> = Vec::new();
    let mut result = Vec::with_capacity(normalized.len());
    let mut intron_columns: HashMap<String, HashSet<usize>> = HashMap::new();

    for (header, seq) in normalized {
        if header.ends_with(ref_suffix) {
            result.push((header.clone(), seq.clone()));
            continue;
        }

        let sb = seq.as_slice();

        // Find contiguous insertion blocks.
        let mut blocks: Vec<(usize, usize)> = Vec::new();
        let mut block_start: Option<usize> = None;
        for i in 0..aln_len {
            let in_insertion = ref_gap_col[i] && sb[i] != b'-';
            if in_insertion {
                if block_start.is_none() {
                    block_start = Some(i);
                }
            } else if let Some(bs) = block_start {
                blocks.push((bs, i));
                block_start = None;
            }
        }
        if let Some(bs) = block_start {
            blocks.push((bs, aln_len));
        }

        // Merge nearby blocks.
        if blocks.len() > 1 {
            let mut merged: Vec<(usize, usize)> = vec![blocks[0]];
            for &(bs, be) in &blocks[1..] {
                let (_, prev_be) = *merged.last().unwrap();
                let gap_residues = sb[prev_be..bs].iter().filter(|&&c| c != b'-').count();
                if gap_residues <= max_bridge_residues {
                    merged.last_mut().unwrap().1 = be;
                } else {
                    merged.push((bs, be));
                }
            }
            blocks = merged;
        }

        let mut cols_to_mask: HashSet<usize> = HashSet::new();
        let candidate_nt = nt_seqs.as_ref().and_then(|m| m.get(header));

        for &(bs, be) in &blocks {
            let block_residues: Vec<u8> =
                sb[bs..be].iter().filter(|&&c| c != b'-').copied().collect();
            if block_residues.len() < min_intron_aa {
                continue;
            }

            // The ref-gap run backing this block must be wide enough to
            // plausibly represent an intron.  A tiny ref-gap (e.g. a
            // 3-column alignment indel) should not trigger a split even
            // if it contains a stop codon.  This mirrors the
            // min_ref_supported_gap check used by the sparse tail trim.
            let ref_gap_cols_in_block = (bs..be).filter(|&i| ref_gap_col[i]).count();
            if ref_gap_cols_in_block < min_ref_supported_gap {
                continue;
            }

            let has_stop = block_residues.contains(&b'*');

            let has_splice = if let Some(nt) = candidate_nt {
                let block_first_residue = sb[..bs].iter().filter(|&&c| c != b'-').count();
                has_splice_sites(nt, block_first_residue, block_residues.len())
            } else {
                false
            };

            if !has_stop && !has_splice {
                continue;
            }

            // Check that most residues sit at ref-gap columns.
            let ref_gap_residues = (bs..be)
                .filter(|&i| sb[i] != b'-' && ref_gap_col[i])
                .count();
            let total_block_residues = (bs..be).filter(|&i| sb[i] != b'-').count();
            if total_block_residues > 0
                && (ref_gap_residues as f64 / total_block_residues as f64) < 0.5
            {
                continue;
            }

            cols_to_mask.extend(bs..be);
            if debug >= 1 {
                let mut reason = Vec::new();
                if has_stop {
                    reason.push("stop");
                }
                if has_splice {
                    reason.push("splice");
                }
                let block_str: String = block_residues.iter().map(|&c| c as char).collect();
                intron_log.push((
                    header.clone(),
                    bs,
                    be,
                    block_residues.len(),
                    block_str,
                    reason.join("+"),
                    ref_gap_cols_in_block,
                    ref_gap_residues,
                    total_block_residues,
                ));
            }
        }

        if cols_to_mask.is_empty() {
            result.push((header.clone(), seq.clone()));
        } else {
            let mut masked = seq.clone();
            for &i in &cols_to_mask {
                if i < masked.len() {
                    masked[i] = b'-';
                }
            }
            intron_columns.insert(header.clone(), cols_to_mask);
            result.push((header.clone(), masked));
        }
    }

    // Debug log
    if !intron_log.is_empty() && debug >= 1 {
        if let Some(ld) = log_dir {
            let log_path = Path::new(ld).join("retained_intron_debug.csv");
            if let Ok(mut f) = File::create(&log_path)
            {
                let _ = writeln!(
                    f,
                    "header,block_start,block_end,residue_count,block_content,reason,ref_gap_cols,ref_gap_residues,total_block_residues"
                );
                for (h, bs, be, rc, content, reason, rgc, rgr, tbr) in &intron_log {
                    let safe = if h.contains(',') {
                        format!("\"{}\"", h)
                    } else {
                        h.clone()
                    };
                    let _ = writeln!(f, "{},{},{},{},{},{},{},{},{}", safe, bs, be, rc, content, reason, rgc, rgr, tbr);
                }
            }
        }
    }

    (result, intron_columns)
}

fn scan_gap_run(seq: &[u8], ref_has_data: &[bool], start: usize) -> (usize, usize) {
    let n = seq.len();
    let mut j = start;
    let mut ref_supported = 0usize;
    while j < n && seq[j] == b'-' {
        if ref_has_data[j] {
            ref_supported += 1;
        }
        j += 1;
    }
    (j, ref_supported)
}

fn block_is_real(
    total_res: usize,
    best_contig: usize,
    anchor_run: usize,
    tail_max_data: usize,
) -> bool {
    best_contig >= anchor_run || total_res > tail_max_data
}

fn pick_final_block_start(
    seq: &[u8],
    ref_has_data: &[bool],
    min_ref_supported_gap: usize,
    anchor_run: usize,
    tail_max_data: usize,
) -> usize {
    let n = seq.len();
    let mut block_start = 0usize;
    let mut total_res = 0usize;
    let mut contig = 0usize;
    let mut best_contig = 0usize;

    let mut best_nonempty_start: Option<usize> = None;
    let mut best_nonempty_res: usize = 0;

    let mut i = 0usize;
    while i < n {
        if seq[i] == b'-' {
            contig = 0;
            let (j, ref_supported) = scan_gap_run(seq, ref_has_data, i);

            if ref_supported >= min_ref_supported_gap {
                if total_res > 0 {
                    if block_is_real(total_res, best_contig, anchor_run, tail_max_data) {
                        return block_start;
                    }
                    if total_res > best_nonempty_res {
                        best_nonempty_res = total_res;
                        best_nonempty_start = Some(block_start);
                    }
                }
                block_start = j;
                total_res = 0;
                best_contig = 0;
            }

            i = j;
            continue;
        }

        total_res += 1;
        contig += 1;
        if contig > best_contig {
            best_contig = contig;
        }
        i += 1;
    }

    if total_res > 0 {
        if block_is_real(total_res, best_contig, anchor_run, tail_max_data) {
            return block_start;
        }
        if total_res > best_nonempty_res {
            best_nonempty_start = Some(block_start);
        }
    }

    best_nonempty_start.unwrap_or(0)
}

fn trim_sparse_tail_end(
    seq: &[u8],
    ref_has_data: &[bool],
    start: usize,
    min_ref_supported_gap: usize,
    _anchor_run: usize,
    tail_max_data: usize,
) -> usize {
    let n = seq.len();
    let mut end = n;
    let mut total_res = 0usize;
    let mut contig = 0usize;
    let mut best_contig = 0usize;

    let mut i = if n > 0 { n - 1 } else { return 0 };
    loop {
        if i >= end {
            if end == 0 || end <= start {
                break;
            }
            i = end - 1;
            continue;
        }
        if i < start {
            break;
        }

        if seq[i] == b'-' {
            contig = 0;
            let mut j = i;
            let mut ref_supported = 0usize;
            while j >= start && seq[j] == b'-' {
                if ref_has_data[j] {
                    ref_supported += 1;
                }
                if j == 0 {
                    break;
                }
                j -= 1;
            }
            // Handle the boundary: if seq[j] != '-' or j < start
            let gap_start = if seq[j] == b'-' { j } else { j + 1 };

            if ref_supported >= min_ref_supported_gap {
                if total_res <= tail_max_data {
                    end = gap_start;
                    total_res = 0;
                    best_contig = 0;
                    contig = 0;
                    if gap_start == 0 || j < start {
                        break;
                    }
                    i = j;
                    continue;
                }
                break;
            }

            if j < start || (j == 0 && seq[0] == b'-') {
                break;
            }
            i = j;
            continue;
        }

        total_res += 1;
        contig += 1;
        if contig > best_contig {
            best_contig = contig;
        }
        if i == 0 || i <= start {
            break;
        }
        i -= 1;
    }

    end
}

fn codon_drops_from_mask(orig_seq: &[u8], masked_columns: &[bool]) -> HashSet<usize> {
    let mut drops = HashSet::new();
    let mut res_idx = 0usize;
    for (col, &aa) in orig_seq.iter().enumerate() {
        if aa != b'-' {
            if col < masked_columns.len() && masked_columns[col] {
                drops.insert(res_idx);
            }
            res_idx += 1;
        }
    }
    drops
}

fn count_edge_residue_trims_from_removed(
    orig_seq: &[u8],
    removed_columns: &[bool],
) -> Option<(usize, usize)> {
    let kept_indices: Vec<usize> = (0..orig_seq.len())
        .filter(|&i| !(i < removed_columns.len() && removed_columns[i]) && orig_seq[i] != b'-')
        .collect();
    if kept_indices.is_empty() {
        return None;
    }
    let first_kept = kept_indices[0];
    let last_kept = *kept_indices.last().unwrap();

    let left_trim = (0..first_kept)
        .filter(|&i| orig_seq[i] != b'-' && i < removed_columns.len() && removed_columns[i])
        .count();
    let right_trim = (last_kept + 1..orig_seq.len())
        .filter(|&i| orig_seq[i] != b'-' && i < removed_columns.len() && removed_columns[i])
        .count();
    Some((left_trim, right_trim))
}

/// Group contiguous column indices into (start, end_exclusive) blocks.
fn group_contiguous_blocks(cols: &HashSet<usize>) -> Vec<(usize, usize)> {
    if cols.is_empty() {
        return Vec::new();
    }
    let mut sorted_cols: Vec<usize> = cols.iter().copied().collect();
    sorted_cols.sort_unstable();

    let mut blocks: Vec<(usize, usize)> = Vec::new();
    let mut bs = sorted_cols[0];
    let mut be = bs;
    for &c in &sorted_cols[1..] {
        if c == be + 1 {
            be = c;
        } else {
            blocks.push((bs, be + 1));
            bs = c;
            be = c;
        }
    }
    blocks.push((bs, be + 1));
    blocks
}

fn classify_intron_splits(orig_seq: &[u8], intron_cols: &HashSet<usize>) -> Vec<(usize, usize)> {
    if intron_cols.is_empty() {
        let total_res = orig_seq.iter().filter(|&&c| c != b'-').count();
        return if total_res > 0 {
            vec![(0, total_res)]
        } else {
            Vec::new()
        };
    }

    let res_positions: Vec<usize> = (0..orig_seq.len()).filter(|&i| orig_seq[i] != b'-').collect();
    if res_positions.is_empty() {
        return Vec::new();
    }
    let total_res = res_positions.len();

    // res_prefix[i] = number of non-gap chars in orig_seq[..i]
    let mut res_prefix = vec![0usize; orig_seq.len() + 1];
    for i in 0..orig_seq.len() {
        res_prefix[i + 1] = res_prefix[i] + if orig_seq[i] != b'-' { 1 } else { 0 };
    }

    // Reuse shared block grouping.
    let blocks = group_contiguous_blocks(intron_cols);

    // Compute all intron residues for subtraction.
    let mut all_intron_residues: HashSet<usize> = HashSet::new();
    for &(block_start, block_end) in &blocks {
        for i in block_start..block_end {
            if orig_seq[i] != b'-' {
                all_intron_residues.insert(res_prefix[i]);
            }
        }
    }

    // Convert blocks to residue-space intervals and classify.
    let mut intron_intervals: Vec<(usize, usize)> = Vec::new();
    for &(block_start, block_end) in &blocks {
        let intron_residues: Vec<usize> = (block_start..block_end)
            .filter(|&i| orig_seq[i] != b'-')
            .collect();
        if intron_residues.is_empty() {
            continue;
        }
        let first_res_col = intron_residues[0];
        let res_start = res_prefix[first_res_col];
        let res_end = res_start + intron_residues.len();

        let has_left = res_start > 0;
        let has_right = res_end < total_res;
        if has_left && has_right {
            intron_intervals.push((res_start, res_end));
        }
    }

    if intron_intervals.is_empty() {
        let kept: HashSet<usize> = (0..total_res)
            .filter(|r| !all_intron_residues.contains(r))
            .collect();
        if kept.is_empty() {
            return Vec::new();
        }
        let first = *kept.iter().min().unwrap();
        let last = *kept.iter().max().unwrap();
        return vec![(first, last - first + 1)];
    }

    let mut kept_residues: Vec<usize> = (0..total_res)
        .filter(|r| !all_intron_residues.contains(r))
        .collect();
    kept_residues.sort_unstable();
    if kept_residues.is_empty() {
        return Vec::new();
    }

    let seq_start = kept_residues[0];
    let seq_end = *kept_residues.last().unwrap() + 1;

    intron_intervals.sort_unstable();
    let mut exons: Vec<(usize, usize)> = Vec::new();
    let mut cursor = seq_start;
    for &(intron_start, intron_end) in &intron_intervals {
        if intron_start > cursor {
            exons.push((cursor, intron_start - cursor));
        }
        cursor = intron_end;
    }
    if cursor < seq_end {
        exons.push((cursor, seq_end - cursor));
    }

    exons
}

// ---------------------------------------------------------------------------
// Main cull_columns function
// ---------------------------------------------------------------------------

#[pyfunction]
#[pyo3(signature = (
    records,
    ref_suffix = ".",
    max_allowed_gaps_in_ref = 0.33,
    min_ref_supported_gap = 10,
    anchor_run = 12,
    tail_max_data = 25,
    include_edge = true,
    edge_ref_data_fraction = 0.5,
    min_seq_len = 10,
    gap_cull_threshold = 1.0,
    nt_seqs = None,
    debug = 0,
    log_dir = None,
))]
pub fn cull_columns(
    records: Vec<(String, String)>,
    ref_suffix: &str,
    max_allowed_gaps_in_ref: f64,
    min_ref_supported_gap: usize,
    anchor_run: usize,
    tail_max_data: usize,
    include_edge: bool,
    edge_ref_data_fraction: f64,
    min_seq_len: usize,
    gap_cull_threshold: f64,
    nt_seqs: Option<HashMap<String, String>>,
    debug: i32,
    log_dir: Option<String>,
) -> PyResult<(
    Vec<(String, String)>,
    HashMap<String, HashSet<usize>>,
    HashMap<String, (usize, usize)>,
    HashMap<String, Vec<(usize, usize)>>,
    HashMap<String, Vec<(usize, usize)>>,
)> {
    if records.is_empty() {
        return Ok((Vec::new(), HashMap::new(), HashMap::new(), HashMap::new(), HashMap::new()));
    }

    let mut gap_cull_log: Vec<(String, usize, usize)> = Vec::new();

    let normalized = normalize_records(&records);
    let orig_seq_map: HashMap<&str, &[u8]> = normalized
        .iter()
        .map(|(h, s)| (h.as_str(), s.as_slice()))
        .collect();
    let orig_aln_len = normalized[0].1.len();
    let mut aln_len = orig_aln_len;

    for (_, seq) in &normalized {
        if seq.len() != aln_len {
            return Err(PyRuntimeError::new_err(
                "ERROR: alignment records not the same length",
            ));
        }
    }


    // ---------------------------------------------------------------
    // Cumulative removal tracker: one Vec<bool> per non-ref record,
    // always in original column space.  Every cull marks into this.
    // ---------------------------------------------------------------
    let mut removed_columns: HashMap<String, Vec<bool>> = HashMap::new();
    let mut intron_removed_original: HashMap<String, HashSet<usize>> = HashMap::new();
    for (h, _) in &normalized {
        if !h.ends_with(ref_suffix) {
            removed_columns.insert(h.clone(), vec![false; orig_aln_len]);
        }
    }

    // ---------------------------------------------------------------
    // Step 1: Ref edge cull + data boundary + sparse tail trim
    // ---------------------------------------------------------------

    let ref_seqs_pre: Vec<&[u8]> = normalized
        .iter()
        .filter(|(h, _)| h.ends_with(ref_suffix))
        .map(|(_, s)| s.as_slice())
        .collect();

    let mut edge_mask = vec![false; aln_len];

    let ref_non_gap: Option<Vec<usize>> = if !ref_seqs_pre.is_empty() {
        let mut rng = vec![0usize; aln_len];
        for seq in &ref_seqs_pre {
            for (i, &ch) in seq.iter().enumerate() {
                if ch != b'-' {
                    rng[i] += 1;
                }
            }
        }
        Some(rng)
    } else {
        None
    };

    if include_edge {
        if let Some(rng) = &ref_non_gap {
            let n_refs = ref_seqs_pre.len() as f64;
            let mut left = 0usize;
            while left < aln_len {
                if (rng[left] as f64 / n_refs) >= edge_ref_data_fraction {
                    break;
                }
                left += 1;
            }
            let mut right = aln_len.saturating_sub(1);
            while right > 0 {
                if (rng[right] as f64 / n_refs) >= edge_ref_data_fraction {
                    break;
                }
                if right == 0 {
                    break;
                }
                right -= 1;
            }
            if left <= right {
                for i in 0..left {
                    edge_mask[i] = true;
                }
                for i in (right + 1)..aln_len {
                    edge_mask[i] = true;
                }
            }
        }
    }

    // Record edge mask.
    for (h, _) in &normalized {
        if h.ends_with(ref_suffix) { continue; }
        if let Some(rv) = removed_columns.get_mut(h) {
            for i in 0..orig_aln_len {
                if edge_mask[i] { rv[i] = true; }
            }
        }
    }

    let ref_has_data: Vec<bool> = if let Some(rng) = &ref_non_gap {
        let n_refs = ref_seqs_pre.len() as f64;
        (0..aln_len)
            .map(|i| (1.0 - (rng[i] as f64 / n_refs)) < max_allowed_gaps_in_ref)
            .collect()
    } else {
        vec![false; aln_len]
    };

    // Per-sequence: edge mask + data boundary + sparse tail.
    let mut structural_masks: Vec<Vec<bool>> = Vec::new();
    for (header, seq) in &normalized {
        let mut combined_mask = edge_mask.clone();

        if !header.ends_with(ref_suffix) && !ref_seqs_pre.is_empty() {
            let first_data = seq.iter().enumerate().position(|(i, &c)| c != b'-' && !edge_mask[i]);
            let (start, end_data) = if let Some(fd) = first_data {
                let last_data = seq.iter().enumerate().rposition(|(i, &c)| c != b'-' && !edge_mask[i]).unwrap();
                (fd, last_data + 1)
            } else {
                (0usize, 0usize)
            };

            for i in 0..start {
                combined_mask[i] = true;
            }
            for i in end_data..aln_len {
                combined_mask[i] = true;
            }

            let eb: Vec<u8> = seq.iter().enumerate().map(|(i, &c)| {
                if combined_mask[i] { b'-' } else { c }
            }).collect();

            let left_trim_point = pick_final_block_start(
                &eb, &ref_has_data, min_ref_supported_gap, anchor_run, tail_max_data,
            );
            let right_trim_point = trim_sparse_tail_end(
                &eb, &ref_has_data, left_trim_point, min_ref_supported_gap, anchor_run, tail_max_data,
            );

            let left_aa_trimmed = (start..left_trim_point).filter(|&i| eb[i] != b'-').count();
            let right_aa_trimmed = (right_trim_point..end_data).filter(|&i| eb[i] != b'-').count();
            if left_aa_trimmed > 0 || right_aa_trimmed > 0 {
                gap_cull_log.push((header.clone(), left_aa_trimmed, right_aa_trimmed));
            }

            for i in start..left_trim_point {
                combined_mask[i] = true;
            }
            for i in right_trim_point..aln_len {
                combined_mask[i] = true;
            }

            // Record into removed_columns.
            if let Some(rv) = removed_columns.get_mut(header) {
                for i in 0..orig_aln_len {
                    if combined_mask[i] { rv[i] = true; }
                }
            }
        }

        structural_masks.push(combined_mask);
    }

    // Apply structural masks.
    let after_structural: Vec<(String, Vec<u8>)> = normalized
        .iter()
        .zip(structural_masks.iter())
        .map(|((h, seq), mask)| {
            let mut masked = seq.clone();
            mask_columns_bool(&mut masked, mask);
            (h.clone(), masked)
        })
        .collect();

    // ---------------------------------------------------------------
    // Step 2: Gap cull (remove all-gap columns from structural trims)
    // ---------------------------------------------------------------

    let gap_cull_mask: Vec<bool> = {
        let mut col_non_gap = vec![0usize; aln_len];
        for (_, seq) in &after_structural {
            for (i, &ch) in seq.iter().enumerate() {
                if ch != b'-' {
                    col_non_gap[i] += 1;
                }
            }
        }
        (0..aln_len).map(|i| col_non_gap[i] == 0).collect()
    };

    // Record gap cull (original space).
    for (h, _) in &normalized {
        if h.ends_with(ref_suffix) { continue; }
        if let Some(rv) = removed_columns.get_mut(h) {
            for i in 0..orig_aln_len {
                if gap_cull_mask[i] { rv[i] = true; }
            }
        }
    }

    // Build column index map: post-gap-cull col -> original col.
    let current_to_orig: Vec<usize> = (0..orig_aln_len)
        .filter(|&i| !gap_cull_mask[i])
        .collect();

    let after_gap_cull: Vec<(String, Vec<u8>)> = if current_to_orig.len() < orig_aln_len {
        after_structural
            .iter()
            .map(|(h, s)| {
                let new_seq: Vec<u8> = current_to_orig.iter().map(|&i| s[i]).collect();
                (h.clone(), new_seq)
            })
            .collect()
    } else {
        after_structural
    };

    aln_len = current_to_orig.len();

    // ---------------------------------------------------------------
    // Step 3: Retained intron detection (on structurally clean aln)
    // ---------------------------------------------------------------

    let (after_intron, intron_masked_cols) =
        mask_retained_introns(&after_gap_cull, ref_suffix, &nt_seqs, min_ref_supported_gap, debug, &log_dir);

    // Record intron masked cols (post-gap-cull space -> original).
    for (header, cols) in &intron_masked_cols {
        for &col in cols {
            if col < current_to_orig.len() {
                let orig_col = current_to_orig[col];
                if let Some(rv) = removed_columns.get_mut(header) {
                    if orig_col < rv.len() { rv[orig_col] = true; }
                }
                intron_removed_original.entry(header.clone()).or_default().insert(orig_col);
            }
        }
    }

    // ---------------------------------------------------------------
    // Step 4: Stop codon BLOSUM trim
    // ---------------------------------------------------------------

    let (after_stop_trim, stop_masked_cols) =
        mask_stop_blosum(&after_intron, ref_suffix, debug, &log_dir);

    // Record stop masked cols (post-gap-cull space -> original).
    for (header, cols) in &stop_masked_cols {
        for &col in cols {
            if col < current_to_orig.len() {
                let orig_col = current_to_orig[col];
                if let Some(rv) = removed_columns.get_mut(header) {
                    if orig_col < rv.len() { rv[orig_col] = true; }
                }
                intron_removed_original.entry(header.clone()).or_default().insert(orig_col);
            }
        }
    }

    // ---------------------------------------------------------------
    // Step 5: Min length filter + global empty column removal
    // ---------------------------------------------------------------

    let mut masked_records: Vec<(String, Vec<u8>)> = Vec::new();
    for (idx, (header, _)) in records.iter().enumerate() {
        if idx >= after_stop_trim.len() { continue; }
        let seq = &after_stop_trim[idx].1;

        if !header.ends_with(ref_suffix) && min_seq_len > 0 {
            let non_gap = seq.iter().filter(|&&c| c != b'-').count();
            if non_gap < min_seq_len {
                continue;
            }
        }

        masked_records.push((header.clone(), seq.clone()));
    }

    if masked_records.is_empty() {
        return Ok((Vec::new(), HashMap::new(), HashMap::new(), HashMap::new(), HashMap::new()));
    }

    // Global empty column removal.
    let masked_len = masked_records[0].1.len();
    let mut global_empty_mask = vec![false; masked_len];
    let mut has_global_empty = false;
    {
        let mut masked_non_gap = vec![0usize; masked_len];
        for (_, seq) in &masked_records {
            for (i, &ch) in seq.iter().enumerate() {
                if ch != b'-' {
                    masked_non_gap[i] += 1;
                }
            }
        }
        for i in 0..masked_len {
            if masked_non_gap[i] == 0 {
                global_empty_mask[i] = true;
                has_global_empty = true;
            }
        }
    }

    // Record global empty (map through current_to_orig).
    if has_global_empty {
        for (i, &is_empty) in global_empty_mask.iter().enumerate() {
            if !is_empty { continue; }
            if i < current_to_orig.len() {
                let orig_col = current_to_orig[i];
                for (h, _) in &normalized {
                    if h.ends_with(ref_suffix) { continue; }
                    if let Some(rv) = removed_columns.get_mut(h) {
                        if orig_col < rv.len() { rv[orig_col] = true; }
                    }
                }
            }
        }
    }

    let trimmed_records: Vec<(String, String)> = if has_global_empty {
        let keep_indices: Vec<usize> = (0..masked_len)
            .filter(|&i| !global_empty_mask[i])
            .collect();
        masked_records.iter().map(|(h, s)| {
            let new_seq: String = keep_indices.iter().map(|&i| s[i] as char).collect();
            (h.clone(), new_seq)
        }).collect()
    } else {
        masked_records.iter().map(|(h, s)| {
            let new_seq: String = s.iter().map(|&c| c as char).collect();
            (h.clone(), new_seq)
        }).collect()
    };

    // ---------------------------------------------------------------
    // Build output maps from removed_columns
    // ---------------------------------------------------------------

    let mut nt_trim_map: HashMap<String, HashSet<usize>> = HashMap::new();
    let mut gff_trim_map: HashMap<String, (usize, usize)> = HashMap::new();

    for (header, removed) in &removed_columns {
        if !masked_records.iter().any(|(h, _)| h == header) { continue; }

        let orig_seq = orig_seq_map.get(header.as_str()).copied().unwrap_or(&[] as &[u8]);
        let drops = codon_drops_from_mask(orig_seq, removed);
        if !drops.is_empty() {
            nt_trim_map.insert(header.clone(), drops);
        }

        // gff_trim_map: edge trims only (removed minus intron cols).
        let mut edge_removed = removed.clone();
        if let Some(ic) = intron_removed_original.get(header.as_str()) {
            for &col in ic {
                if col < orig_aln_len { edge_removed[col] = false; }
            }
        }
        if let Some(trims) = count_edge_residue_trims_from_removed(orig_seq, &edge_removed) {
            if trims.0 > 0 || trims.1 > 0 {
                gff_trim_map.insert(header.clone(), trims);
            }
        }
    }

    // Build gff_intron_map and exon_split_map from intron_removed_original
    // (already in original column space).
    let mut gff_intron_map: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    let mut exon_split_map: HashMap<String, Vec<(usize, usize)>> = HashMap::new();

    for (header, orig_cols) in &intron_removed_original {
        let orig_seq = match orig_seq_map.get(header.as_str()) {
            Some(s) => *s,
            None => continue,
        };

        let blocks = group_contiguous_blocks(orig_cols);
        let mut splits: Vec<(usize, usize)> = Vec::new();
        for &(block_start, block_end) in &blocks {
            let res_offset = orig_seq[..block_start].iter().filter(|&&c| c != b'-').count();
            let intron_res = (block_start..block_end).filter(|&i| orig_seq[i] != b'-').count();
            if intron_res > 0 { splits.push((res_offset, intron_res)); }
        }
        if !splits.is_empty() { gff_intron_map.insert(header.clone(), splits); }

        let exons = classify_intron_splits(orig_seq, orig_cols);
        if exons.len() > 1 { exon_split_map.insert(header.clone(), exons); }
    }
    if !gap_cull_log.is_empty() && debug >= 1 {
        if let Some(ld) = &log_dir {
            let log_path = Path::new(ld).join("gap_cull_debug.csv");
            if let Ok(mut f) = File::create(&log_path)
            {
                let _ = writeln!(f, "header,left_aa_trimmed,right_aa_trimmed");
                for (h, left, right) in &gap_cull_log {
                    let safe = if h.contains(',') {
                        format!("\"{}\"", h)
                    } else {
                        h.clone()
                    };
                    let _ = writeln!(f, "{},{},{}", safe, left, right);
                }
            }
        }
    }

    Ok((
        trimmed_records,
        nt_trim_map,
        gff_trim_map,
        gff_intron_map,
        exon_split_map,
    ))
}

// ---------------------------------------------------------------------------
// apply_gff_culls
// ---------------------------------------------------------------------------

fn extract_gff_attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    let prefix = format!("{}=", key);
    let start = attrs.find(&prefix)?;
    let val_start = start + prefix.len();
    let val_end = attrs[val_start..]
        .find(|c: char| c == ';' || c == ' ')
        .map(|i| val_start + i)
        .unwrap_or(attrs.len());
    Some(&attrs[val_start..val_end])
}

#[pyfunction]
#[pyo3(signature = (source_gff_path, output_gff_path, culls, intron_splits = None))]
pub fn apply_gff_culls(
    source_gff_path: &str,
    output_gff_path: &str,
    culls: HashMap<(String, String, Option<String>), Option<(usize, usize)>>,
    intron_splits: Option<
        HashMap<(String, String, Option<String>), Vec<(usize, usize)>>,
    >,
) -> PyResult<(bool, HashSet<String>)> {
    let intron_splits = intron_splits.unwrap_or_default();

    if culls.is_empty() && intron_splits.is_empty() {
        return Ok((false, HashSet::new()));
    }

    let file = File::open(source_gff_path)
        .map_err(|e| PyRuntimeError::new_err(format!("Cannot open GFF: {}", e)))?;
    let reader = BufReader::new(file);

    let mut gff_lines: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

        if line.starts_with('#') {
            gff_lines.push(line);
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let attributes = fields[8];
        let name = match extract_gff_attr(attributes, "Name") {
            Some(n) => n.to_string(),
            None => {
                // No modification needed, push original line directly.
                gff_lines.push(line);
                continue;
            }
        };
        let gene = match extract_gff_attr(attributes, "Parent") {
            Some(g) => g.to_string(),
            None => {
                gff_lines.push(line);
                continue;
            }
        };
        let frame: Option<String> = extract_gff_attr(attributes, "Note")
            .map(|n| n.split(',').next().unwrap_or("").to_string());

        let trim_key = (gene.clone(), name.clone(), frame.clone());

        // Kicked entirely?
        if let Some(val) = culls.get(&trim_key) {
            if val.is_none() {
                continue;
            }
        }

        let mut start: i64 = fields[3]
            .parse()
            .map_err(|_| PyRuntimeError::new_err("Bad GFF start coordinate"))?;
        let mut end: i64 = fields[4]
            .parse()
            .map_err(|_| PyRuntimeError::new_err("Bad GFF end coordinate"))?;
        let strand = fields[6];

        // Apply edge trims.
        let cull_vals = culls.get(&trim_key).and_then(|v| *v);
        if let Some((left_trim, right_trim)) = cull_vals {
            let trim_left_bp = (left_trim * 3) as i64;
            let trim_right_bp = (right_trim * 3) as i64;
            if strand == "-" {
                start += trim_right_bp;
                end -= trim_left_bp;
            } else {
                start += trim_left_bp;
                end -= trim_right_bp;
            }
            if start > end {
                continue;
            }
        }

        // Apply intron splits.
        if let Some(splits) = intron_splits.get(&trim_key) {
            let left_trim_residues = cull_vals.map(|v| v.0).unwrap_or(0) as i64;
            let mut sorted_splits = splits.clone();
            sorted_splits.sort_by_key(|&(o, _)| o);

            let mut segments: Vec<(i64, i64)> = Vec::new();

            if strand == "-" {
                let mut cursor = end;
                for &(res_offset, intron_res) in &sorted_splits {
                    let mut adjusted_offset = res_offset as i64 - left_trim_residues;
                    let mut ir = intron_res as i64;
                    if adjusted_offset < 0 {
                        ir += adjusted_offset;
                        adjusted_offset = 0;
                    }
                    if ir <= 0 {
                        continue;
                    }
                    let intron_top = end - adjusted_offset * 3;
                    let intron_bottom = intron_top - ir * 3 + 1;
                    if cursor > intron_top {
                        segments.push((intron_top + 1, cursor));
                    }
                    cursor = intron_bottom - 1;
                }
                if cursor >= start {
                    segments.push((start, cursor));
                }
            } else {
                let mut cursor = start;
                for &(res_offset, intron_res) in &sorted_splits {
                    let mut adjusted_offset = res_offset as i64 - left_trim_residues;
                    let mut ir = intron_res as i64;
                    if adjusted_offset < 0 {
                        ir += adjusted_offset;
                        adjusted_offset = 0;
                    }
                    if ir <= 0 {
                        continue;
                    }
                    let intron_start_coord = start + adjusted_offset * 3;
                    let intron_end_coord = intron_start_coord + ir * 3 - 1;
                    if intron_start_coord > cursor {
                        segments.push((cursor, intron_start_coord - 1));
                    }
                    cursor = intron_end_coord + 1;
                }
                if cursor <= end {
                    segments.push((cursor, end));
                }
            }

            let valid_segments: Vec<(i64, i64)> =
                segments.into_iter().filter(|&(s, e)| s <= e).collect();
            let multi = valid_segments.len() > 1;

            let labeled_segments = if strand == "-" && multi {
                valid_segments.into_iter().rev().collect::<Vec<_>>()
            } else {
                valid_segments
            };

            for (seg_i, (seg_start, seg_end)) in labeled_segments.iter().enumerate() {
                let mut row: Vec<String> = fields.iter().map(|s| s.to_string()).collect();
                row[3] = seg_start.to_string();
                row[4] = seg_end.to_string();
                if multi {
                    // Replace Name=XXX with Name=XXX_E{n}
                    let new_attrs = row[8].replace(
                        &format!("Name={}", name),
                        &format!("Name={}_E{}", name, seg_i + 1),
                    );
                    row[8] = new_attrs;
                }
                gff_lines.push(row.join("\t"));
            }
        } else {
            // Avoid re-splitting+joining when no coordinate changes needed.
            let orig_start: i64 = fields[3].parse().unwrap_or(0);
            let orig_end: i64 = fields[4].parse().unwrap_or(0);
            if start == orig_start && end == orig_end {
                gff_lines.push(line);
            } else {
                let mut row: Vec<String> = fields.iter().map(|s| s.to_string()).collect();
                row[3] = start.to_string();
                row[4] = end.to_string();
                gff_lines.push(row.join("\t"));
            }
        }
    }

    // Genomic containment post-pass.
    let mut contained_kicked_names: HashSet<String> = HashSet::new();

    struct ParsedRow {
        line_idx: usize,
        _contig: String,
        start: i64,
        end: i64,
        _strand: String,
        gene: String,
        name: String,
    }

    let mut parsed_rows: Vec<ParsedRow> = Vec::new();
    for (li, line) in gff_lines.iter().enumerate() {
        if line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 {
            continue;
        }
        let nm = match extract_gff_attr(parts[8], "Name") {
            Some(n) => n.to_string(),
            None => continue,
        };
        let pm = match extract_gff_attr(parts[8], "Parent") {
            Some(p) => p.to_string(),
            None => continue,
        };
        let s: i64 = match parts[3].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let e: i64 = match parts[4].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        parsed_rows.push(ParsedRow {
            line_idx: li,
            _contig: parts[0].to_string(),
            start: s,
            end: e,
            _strand: parts[6].to_string(),
            gene: pm,
            name: nm,
        });
    }

    // Group by (gene, contig, strand).
    let mut groups: HashMap<(String, String, String), Vec<usize>> = HashMap::new();
    for (ri, row) in parsed_rows.iter().enumerate() {
        let key = (row.gene.clone(), row._contig.clone(), row._strand.clone());
        groups.entry(key).or_default().push(ri);
    }

    // O(n log n) containment detection instead of O(n^2).
    let mut drop_indices: HashSet<usize> = HashSet::new();
    for members in groups.values() {
        if members.len() < 2 {
            continue;
        }
        // Sort by start ascending, then by end descending (wider intervals first).
        let mut sorted_members: Vec<usize> = members.clone();
        sorted_members.sort_by(|&a, &b| {
            let ra = &parsed_rows[a];
            let rb = &parsed_rows[b];
            ra.start.cmp(&rb.start).then(rb.end.cmp(&ra.end))
        });

        // Walk the sorted list. For each interval, if its end fits within a
        // previously seen interval's end, it is contained. Because we sorted
        // by (start ASC, end DESC), the first interval at each start position
        // is the widest, so a simple running max_end sweep works.
        //
        // Edge case: the original O(n^2) code drops BOTH intervals when they
        // are identical (each contains the other). We replicate that here by
        // tracking duplicate (start, end) pairs.
        let mut max_end = i64::MIN;
        let mut max_end_idx: Option<usize> = None;
        for &mi in &sorted_members {
            let row = &parsed_rows[mi];
            if row.end <= max_end {
                // Contained by a previously seen wider (or equal) interval.
                drop_indices.insert(row.line_idx);
                contained_kicked_names.insert(row.name.clone());
                // If this interval is identical to the one that set max_end,
                // the original code would also drop that one.
                if let Some(prev_mi) = max_end_idx {
                    let prev = &parsed_rows[prev_mi];
                    if row.start == prev.start && row.end == prev.end {
                        drop_indices.insert(prev.line_idx);
                        contained_kicked_names.insert(prev.name.clone());
                    }
                }
            }
            if row.end > max_end {
                max_end = row.end;
                max_end_idx = Some(mi);
            }
        }
    }

    if !drop_indices.is_empty() {
        gff_lines = gff_lines
            .into_iter()
            .enumerate()
            .filter(|(i, _)| !drop_indices.contains(i))
            .map(|(_, l)| l)
            .collect();
    }

    if !gff_lines.is_empty() {
        let mut out = File::create(output_gff_path)
            .map_err(|e| PyRuntimeError::new_err(format!("Cannot create GFF output: {}", e)))?;
        write!(out, "{}", gff_lines.join("\n"))
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok((true, contained_kicked_names))
    } else {
        Ok((false, contained_kicked_names))
    }
}
