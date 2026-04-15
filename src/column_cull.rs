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

fn is_donor(di: &str) -> bool {
    di == "GT" || di == "GC"
}

fn is_acceptor(di: &str) -> bool {
    di == "AG"
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn normalize_records(records: &[(String, String)]) -> Vec<(String, String)> {
    records
        .iter()
        .map(|(h, s)| {
            let norm = s
                .replace('.', "-")
                .replace(' ', "")
                .to_ascii_uppercase();
            (h.clone(), norm)
        })
        .collect()
}

fn mask_columns(seq: &str, columns: &HashSet<usize>) -> String {
    if columns.is_empty() {
        return seq.to_string();
    }
    let mut chars: Vec<u8> = seq.bytes().collect();
    for &i in columns {
        if i < chars.len() {
            chars[i] = b'-';
        }
    }
    String::from_utf8(chars).unwrap()
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
        let donor = std::str::from_utf8(&nt[d_start..d_start + 2])
            .unwrap_or("")
            .to_ascii_uppercase();
        let acceptor = std::str::from_utf8(&nt[a_end - 2..a_end])
            .unwrap_or("")
            .to_ascii_uppercase();
        if is_donor(&donor) && is_acceptor(&acceptor) {
            return true;
        }
    }
    false
}

/// Detect and mask retained same-frame introns in candidate sequences.
fn mask_retained_introns(
    normalized: &[(String, String)],
    ref_suffix: &str,
    nt_seqs: &Option<HashMap<String, String>>,
    debug: i32,
    log_dir: &Option<String>,
) -> (Vec<(String, String)>, HashMap<String, HashSet<usize>>) {
    let min_intron_aa: usize = 3;
    let ref_gap_threshold: f64 = 0.95;
    let max_bridge_residues: usize = 8;

    if normalized.is_empty() {
        return (Vec::new(), HashMap::new());
    }

    let aln_len = normalized[0].1.len();
    let ref_seqs: Vec<&[u8]> = normalized
        .iter()
        .filter(|(h, _)| h.ends_with(ref_suffix))
        .map(|(_, s)| s.as_bytes())
        .collect();

    if ref_seqs.is_empty() {
        return (normalized.to_vec(), HashMap::new());
    }

    let n_refs = ref_seqs.len();

    // Pre-compute reference-gap columns.
    let mut ref_gap_col = vec![false; aln_len];
    for i in 0..aln_len {
        let gap_count = ref_seqs.iter().filter(|s| s[i] == b'-').count();
        ref_gap_col[i] = (gap_count as f64 / n_refs as f64) >= ref_gap_threshold;
    }

    let mut intron_log: Vec<(String, usize, usize, usize, String, String)> = Vec::new();
    let mut result = Vec::with_capacity(normalized.len());
    let mut intron_columns: HashMap<String, HashSet<usize>> = HashMap::new();

    for (header, seq) in normalized {
        if header.ends_with(ref_suffix) {
            result.push((header.clone(), seq.clone()));
            continue;
        }

        let sb = seq.as_bytes();

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
            let block_residues: String =
                sb[bs..be].iter().filter(|&&c| c != b'-').map(|&c| c as char).collect();
            if block_residues.len() < min_intron_aa {
                continue;
            }

            let has_stop = block_residues.contains('*');

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
                intron_log.push((
                    header.clone(),
                    bs,
                    be,
                    block_residues.len(),
                    block_residues.clone(),
                    reason.join("+"),
                ));
            }
        }

        if cols_to_mask.is_empty() {
            result.push((header.clone(), seq.clone()));
        } else {
            let masked = mask_columns(seq, &cols_to_mask);
            intron_columns.insert(header.clone(), cols_to_mask);
            result.push((header.clone(), masked));
        }
    }

    // Debug log
    if !intron_log.is_empty() && debug >= 1 {
        if let Some(ld) = log_dir {
            let log_path = Path::new(ld).join("retained_intron_debug.csv");
            let write_header = !log_path.exists();
            if let Ok(mut f) = std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(&log_path)
            {
                if write_header {
                    let _ = writeln!(
                        f,
                        "header,block_start,block_end,residue_count,block_content,reason"
                    );
                }
                for (h, bs, be, rc, content, reason) in &intron_log {
                    let safe = if h.contains(',') {
                        format!("\"{}\"", h)
                    } else {
                        h.clone()
                    };
                    let _ = writeln!(f, "{},{},{},{},{},{}", safe, bs, be, rc, content, reason);
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

fn codon_drops_from_mask(orig_seq: &str, masked_columns: &HashSet<usize>) -> HashSet<usize> {
    let mut drops = HashSet::new();
    let mut res_idx = 0usize;
    for (col, aa) in orig_seq.bytes().enumerate() {
        if aa != b'-' {
            if masked_columns.contains(&col) {
                drops.insert(res_idx);
            }
            res_idx += 1;
        }
    }
    drops
}

fn count_edge_residue_trims_from_removed(
    orig_seq: &str,
    removed_columns: &HashSet<usize>,
) -> Option<(usize, usize)> {
    let ob = orig_seq.as_bytes();
    let kept_indices: Vec<usize> = (0..ob.len())
        .filter(|&i| !removed_columns.contains(&i) && ob[i] != b'-')
        .collect();
    if kept_indices.is_empty() {
        return None;
    }
    let first_kept = kept_indices[0];
    let last_kept = *kept_indices.last().unwrap();

    let left_trim = (0..first_kept)
        .filter(|&i| ob[i] != b'-' && removed_columns.contains(&i))
        .count();
    let right_trim = (last_kept + 1..ob.len())
        .filter(|&i| ob[i] != b'-' && removed_columns.contains(&i))
        .count();
    Some((left_trim, right_trim))
}

fn classify_intron_splits(orig_seq: &str, intron_cols: &HashSet<usize>) -> Vec<(usize, usize)> {
    let ob = orig_seq.as_bytes();

    if intron_cols.is_empty() {
        let total_res = ob.iter().filter(|&&c| c != b'-').count();
        return if total_res > 0 {
            vec![(0, total_res)]
        } else {
            Vec::new()
        };
    }

    let res_positions: Vec<usize> = (0..ob.len()).filter(|&i| ob[i] != b'-').collect();
    if res_positions.is_empty() {
        return Vec::new();
    }
    let total_res = res_positions.len();

    // res_prefix[i] = number of non-gap chars in orig_seq[..i]
    let mut res_prefix = vec![0usize; ob.len() + 1];
    for i in 0..ob.len() {
        res_prefix[i + 1] = res_prefix[i] + if ob[i] != b'-' { 1 } else { 0 };
    }

    // Group contiguous intron columns into blocks.
    let mut sorted_cols: Vec<usize> = intron_cols.iter().copied().collect();
    sorted_cols.sort();

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

    // Convert blocks to residue-space intervals and classify.
    let mut intron_intervals: Vec<(usize, usize)> = Vec::new();
    for &(block_start, block_end) in &blocks {
        let intron_residues: Vec<usize> = (block_start..block_end)
            .filter(|&i| ob[i] != b'-')
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

    // Compute all intron residues for subtraction.
    let mut all_intron_residues: HashSet<usize> = HashSet::new();
    for &(block_start, block_end) in &blocks {
        for i in block_start..block_end {
            if ob[i] != b'-' {
                all_intron_residues.insert(res_prefix[i]);
            }
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
    kept_residues.sort();
    if kept_residues.is_empty() {
        return Vec::new();
    }

    let seq_start = kept_residues[0];
    let seq_end = *kept_residues.last().unwrap() + 1;

    intron_intervals.sort();
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
    let orig_seq_map: HashMap<&str, &str> = normalized
        .iter()
        .map(|(h, s)| (h.as_str(), s.as_str()))
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

    // Detect and mask retained introns.
    let (after_intron, intron_masked_cols) =
        mask_retained_introns(&normalized, ref_suffix, &nt_seqs, debug, &log_dir);

    // Gap cull: remove all-gap columns.
    let mut gap_cull_columns: HashSet<usize> = HashSet::new();
    if gap_cull_threshold > 0.0 {
        let mut col_non_gap = vec![0usize; aln_len];
        for (_, seq) in &after_intron {
            for (i, ch) in seq.bytes().enumerate() {
                if ch != b'-' {
                    col_non_gap[i] += 1;
                }
            }
        }
        gap_cull_columns = (0..aln_len).filter(|&i| col_non_gap[i] == 0).collect();
    }

    let after_gap_cull: Vec<(String, String)>;
    let orig_index_map: Option<Vec<usize>>;

    if !gap_cull_columns.is_empty() {
        let oimap: Vec<usize> = (0..aln_len)
            .filter(|i| !gap_cull_columns.contains(i))
            .collect();
        after_gap_cull = after_intron
            .iter()
            .map(|(h, s)| {
                let sb = s.as_bytes();
                let new_seq: String = oimap.iter().map(|&i| sb[i] as char).collect();
                (h.clone(), new_seq)
            })
            .collect();
        aln_len = oimap.len();
        orig_index_map = Some(oimap);
    } else {
        after_gap_cull = after_intron;
        orig_index_map = None;
    }

    let current_to_orig: Vec<usize> = match &orig_index_map {
        Some(m) => m.clone(),
        None => (0..aln_len).collect(),
    };

    let ref_seqs: Vec<&[u8]> = after_gap_cull
        .iter()
        .filter(|(h, _)| h.ends_with(ref_suffix))
        .map(|(_, s)| s.as_bytes())
        .collect();

    let mut edge_columns: HashSet<usize> = HashSet::new();

    // Compute ref_non_gap
    let ref_non_gap: Option<Vec<usize>> = if !ref_seqs.is_empty() {
        let mut rng = vec![0usize; aln_len];
        for seq in &ref_seqs {
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

    // Edge trimming.
    if include_edge {
        if let Some(rng) = &ref_non_gap {
            let n_refs = ref_seqs.len() as f64;
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
                edge_columns.extend(0..left);
                edge_columns.extend(right + 1..aln_len);
            }
        }
    }

    // ref_has_data
    let ref_has_data: Vec<bool> = if let Some(rng) = &ref_non_gap {
        let n_refs = ref_seqs.len() as f64;
        (0..aln_len)
            .map(|i| (1.0 - (rng[i] as f64 / n_refs)) < max_allowed_gaps_in_ref)
            .collect()
    } else {
        vec![false; aln_len]
    };

    let mut masked_records: Vec<(String, String)> = Vec::new();
    let mut trim_sets: HashMap<String, HashSet<usize>> = HashMap::new();

    // Process each record.
    for (idx, (header, _)) in records.iter().enumerate() {
        let (_, seq) = &after_gap_cull[idx];
        let sb = seq.as_bytes();

        let edge_masked = mask_columns(seq, &edge_columns);
        let eb = edge_masked.as_bytes();
        let mut to_trim: HashSet<usize> = edge_columns.clone();

        if !header.ends_with(ref_suffix) && !ref_seqs.is_empty() {
            let first_data = eb.iter().position(|&c| c != b'-');
            let (start, end_data) = if let Some(fd) = first_data {
                let last_data = eb.len()
                    - 1
                    - eb.iter().rev().position(|&c| c != b'-').unwrap();
                (fd, last_data + 1)
            } else {
                (0usize, 0usize)
            };

            // Trim outside data region.
            to_trim.extend(0..start);
            to_trim.extend(end_data..eb.len());

            let effective_len = if end_data > start { end_data - start } else { 0 };
            let _limit_len = std::cmp::min((effective_len as f64 * 0.15) as usize, 40);

            let left_trim_point = pick_final_block_start(
                eb,
                &ref_has_data,
                min_ref_supported_gap,
                anchor_run,
                tail_max_data,
            );

            let right_trim_point = trim_sparse_tail_end(
                eb,
                &ref_has_data,
                left_trim_point,
                min_ref_supported_gap,
                anchor_run,
                tail_max_data,
            );

            // Count AA trimmed for debug.
            let left_aa_trimmed = (start..left_trim_point)
                .filter(|&i| eb[i] != b'-')
                .count();
            let right_aa_trimmed = (right_trim_point..end_data)
                .filter(|&i| eb[i] != b'-')
                .count();
            if left_aa_trimmed > 0 || right_aa_trimmed > 0 {
                gap_cull_log.push((header.clone(), left_aa_trimmed, right_aa_trimmed));
            }

            to_trim.extend(start..left_trim_point);
            to_trim.extend(right_trim_point..eb.len());
        }

        let masked_seq = mask_columns(&edge_masked, &to_trim);

        // Min length filter.
        if !header.ends_with(ref_suffix) && min_seq_len > 0 {
            let non_gap = masked_seq.bytes().filter(|&c| c != b'-').count();
            if non_gap < min_seq_len {
                continue;
            }
        }

        masked_records.push((header.clone(), masked_seq));
        trim_sets.insert(header.clone(), to_trim);
    }

    if masked_records.is_empty() {
        return Ok((Vec::new(), HashMap::new(), HashMap::new(), HashMap::new(), HashMap::new()));
    }

    // Global empty column removal.
    let masked_len = masked_records[0].1.len();
    let mut masked_non_gap = vec![0usize; masked_len];
    for (_, seq) in &masked_records {
        for (i, ch) in seq.bytes().enumerate() {
            if ch != b'-' {
                masked_non_gap[i] += 1;
            }
        }
    }
    let global_empty_columns: HashSet<usize> =
        (0..masked_len).filter(|&i| masked_non_gap[i] == 0).collect();

    let trimmed_records: Vec<(String, String)> = if !global_empty_columns.is_empty() {
        let keep_indices: Vec<usize> = (0..masked_len)
            .filter(|i| !global_empty_columns.contains(i))
            .collect();
        masked_records
            .iter()
            .map(|(h, s)| {
                let sb = s.as_bytes();
                let new_seq: String = keep_indices.iter().map(|&i| sb[i] as char).collect();
                (h.clone(), new_seq)
            })
            .collect()
    } else {
        masked_records.clone()
    };

    let global_empty_original: HashSet<usize> = global_empty_columns
        .iter()
        .filter(|&&i| i < current_to_orig.len())
        .map(|&i| current_to_orig[i])
        .collect();

    // Build nt_trim_map, gff_trim_map.
    let mut nt_trim_map: HashMap<String, HashSet<usize>> = HashMap::new();
    let mut gff_trim_map: HashMap<String, (usize, usize)> = HashMap::new();

    for (header, _) in &masked_records {
        if header.ends_with(ref_suffix) {
            continue;
        }

        let to_trim = trim_sets.get(header.as_str()).cloned().unwrap_or_default();
        let mut to_trim_original: HashSet<usize> = to_trim
            .iter()
            .filter(|&&i| i < current_to_orig.len())
            .map(|&i| current_to_orig[i])
            .collect();

        let mut removed_original: HashSet<usize> = gap_cull_columns.clone();
        removed_original.extend(to_trim_original.iter());
        removed_original.extend(global_empty_original.iter());

        let intron_cols = intron_masked_cols.get(header.as_str());
        if let Some(ic) = intron_cols {
            removed_original.extend(ic.iter());
        }

        let orig_seq = orig_seq_map.get(header.as_str()).unwrap_or(&"");
        let drops = codon_drops_from_mask(orig_seq, &removed_original);
        if !drops.is_empty() {
            nt_trim_map.insert(header.clone(), drops);
        }

        let edge_removed: HashSet<usize> = if let Some(ic) = intron_cols {
            removed_original.difference(ic).copied().collect()
        } else {
            removed_original.clone()
        };
        if let Some(trims) = count_edge_residue_trims_from_removed(orig_seq, &edge_removed) {
            if trims.0 > 0 || trims.1 > 0 {
                gff_trim_map.insert(header.clone(), trims);
            }
        }
    }

    // Build gff_intron_map.
    let mut gff_intron_map: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    for (header, cols) in &intron_masked_cols {
        let orig_seq = match orig_seq_map.get(header.as_str()) {
            Some(s) => s.as_bytes(),
            None => continue,
        };

        let mut sorted_cols: Vec<usize> = cols.iter().copied().collect();
        sorted_cols.sort();

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

        let mut splits: Vec<(usize, usize)> = Vec::new();
        for &(block_start, block_end) in &blocks {
            let res_offset = orig_seq[..block_start]
                .iter()
                .filter(|&&c| c != b'-')
                .count();
            let intron_res = (block_start..block_end)
                .filter(|&i| orig_seq[i] != b'-')
                .count();
            if intron_res > 0 {
                splits.push((res_offset, intron_res));
            }
        }
        if !splits.is_empty() {
            gff_intron_map.insert(header.clone(), splits);
        }
    }

    // Build exon_split_map.
    let mut exon_split_map: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    for (header, cols) in &intron_masked_cols {
        let orig_seq = match orig_seq_map.get(header.as_str()) {
            Some(s) => *s,
            None => continue,
        };
        let exons = classify_intron_splits(orig_seq, cols);
        if exons.len() > 1 {
            exon_split_map.insert(header.clone(), exons);
        }
    }

    // Debug CSV
    if !gap_cull_log.is_empty() && debug >= 1 {
        if let Some(ld) = &log_dir {
            let log_path = Path::new(ld).join("gap_cull_debug.csv");
            let write_header = !log_path.exists();
            if let Ok(mut f) = std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(&log_path)
            {
                if write_header {
                    let _ = writeln!(f, "header,left_aa_trimmed,right_aa_trimmed");
                }
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
                gff_lines.push(fields.join("\t"));
                continue;
            }
        };
        let gene = match extract_gff_attr(attributes, "Parent") {
            Some(g) => g.to_string(),
            None => {
                gff_lines.push(fields.join("\t"));
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
            let mut row: Vec<String> = fields.iter().map(|s| s.to_string()).collect();
            row[3] = start.to_string();
            row[4] = end.to_string();
            gff_lines.push(row.join("\t"));
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

    let mut drop_indices: HashSet<usize> = HashSet::new();
    for members in groups.values() {
        if members.len() < 2 {
            continue;
        }
        for &ai in members {
            for &bi in members {
                if ai == bi {
                    continue;
                }
                let a = &parsed_rows[ai];
                let b = &parsed_rows[bi];
                if a.start >= b.start && a.end <= b.end {
                    drop_indices.insert(a.line_idx);
                    contained_kicked_names.insert(a.name.clone());
                }
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
