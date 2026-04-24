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

    let mut intron_log: Vec<(String, usize, usize, usize, String, String)> = Vec::new();
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

/// Mask garbage regions around internal stop codons.
///
/// Uses the stop codon as a seed for a known-bad position, then expands
/// outward in both directions using a sliding window of `recovery_window`
/// scorable residues.  Masking continues until a window is found where
/// at least `recovery_threshold` fraction of positions match the
/// reference character set.  The recovery point is then advanced past
/// any leading mismatches in the recovery window so that only truly
/// matching sequence is kept.
///
/// If the scan reaches the edge of the sequence data without finding a
/// recovery window, the entire side is masked.  The smaller surviving
/// fragment is tossed.
///
/// Columns masked here are merged into `intron_columns` so that
/// downstream NT trimming and GFF coordinate adjustment account for
/// the removal.  Because the orphan eviction guarantees only one
/// contiguous fragment survives, the intron split logic will not
/// produce multiple exons from these columns.
fn mask_stop_sides(
    records: &[(String, Vec<u8>)],
    intron_columns: &mut HashMap<String, HashSet<usize>>,
    ref_suffix: &str,
    recovery_threshold: f64,
    recovery_window: usize,
    debug: i32,
    log_dir: &Option<String>,
) -> Vec<(String, Vec<u8>)> {
    if records.is_empty() {
        return Vec::new();
    }

    let aln_len = records[0].1.len();
    let ref_seqs: Vec<&[u8]> = records
        .iter()
        .filter(|(h, _)| h.ends_with(ref_suffix))
        .map(|(_, s)| s.as_slice())
        .collect();

    if ref_seqs.is_empty() {
        return records.to_vec();
    }

    // Build per-column ref character sets, considering only positions
    // inside each reference's data range (between first and last
    // non-gap character).
    let mut ref_char_sets: Vec<HashSet<u8>> = vec![HashSet::new(); aln_len];
    for rseq in &ref_seqs {
        let first = rseq.iter().position(|&c| c != b'-');
        let last = rseq.iter().rposition(|&c| c != b'-');
        if let (Some(f), Some(l)) = (first, last) {
            for i in f..=l {
                let c = rseq[i];
                if c != b'-' && c != b'*' {
                    ref_char_sets[i].insert(c);
                }
            }
        }
    }

    let mut stop_log: Vec<(String, usize, String, usize)> = Vec::new();
    let mut result = Vec::with_capacity(records.len());

    for (header, seq) in records {
        if header.ends_with(ref_suffix) {
            result.push((header.clone(), seq.clone()));
            continue;
        }

        // Fast path: no stops remaining.
        if !seq.contains(&b'*') {
            result.push((header.clone(), seq.clone()));
            continue;
        }

        // Find data range.
        let data_start = match seq.iter().position(|&c| c != b'-') {
            Some(p) => p,
            None => {
                result.push((header.clone(), seq.clone()));
                continue;
            }
        };
        let data_end = seq.iter().rposition(|&c| c != b'-').unwrap() + 1;

        // Collect stop positions.
        let stops: Vec<usize> = (data_start..data_end)
            .filter(|&i| seq[i] == b'*')
            .collect();

        if stops.is_empty() {
            result.push((header.clone(), seq.clone()));
            continue;
        }

        // Helper: collect scorable (column, matches_ref) pairs.
        let collect_scorable = |range_iter: Box<dyn Iterator<Item = usize>>| -> Vec<(usize, bool)> {
            range_iter
                .filter_map(|i| {
                    let c = seq[i];
                    if c == b'-' || c == b'*' {
                        return None;
                    }
                    if ref_char_sets[i].is_empty() {
                        return None;
                    }
                    Some((i, ref_char_sets[i].contains(&c)))
                })
                .collect()
        };

        // Helper: scan outward from a stop, return the column where
        // good sequence resumes (or the data edge if it never does).
        let find_recovery = |positions: &[(usize, bool)], toward_start: bool| -> usize {
            if positions.len() < recovery_window {
                // Not enough scorable positions to fill a window.
                return if toward_start { data_start } else { data_end };
            }
            let required = (recovery_window as f64 * recovery_threshold).ceil() as usize;
            for w in 0..=(positions.len() - recovery_window) {
                let w_matches = positions[w..w + recovery_window]
                    .iter()
                    .filter(|(_, m)| *m)
                    .count();
                if w_matches >= required {
                    // Recovery window found.  Advance past leading
                    // mismatches so only matching residues are kept.
                    for j in w..w + recovery_window {
                        if positions[j].1 {
                            return positions[j].0;
                        }
                    }
                }
            }
            // No recovery found.
            if toward_start { data_start } else { data_end }
        };

        let mut cols_to_mask: HashSet<usize> = HashSet::new();

        for &star in &stops {
            // Always mask the stop itself.
            cols_to_mask.insert(star);

            // Scan RIGHT from stop.
            let right_positions = collect_scorable(
                Box::new((star + 1)..data_end),
            );
            let recovery_right = find_recovery(&right_positions, false);
            let right_had_garbage = recovery_right > star + 1;
            if right_had_garbage {
                cols_to_mask.extend((star + 1)..recovery_right);
                if debug >= 1 {
                    stop_log.push((
                        header.clone(),
                        star,
                        "right".to_string(),
                        recovery_right,
                    ));
                }
            }

            // Scan LEFT from stop.
            let left_positions = collect_scorable(
                Box::new((data_start..star).rev()),
            );
            let recovery_left = find_recovery(&left_positions, true);
            let left_had_garbage;
            if recovery_left == data_start {
                // No recovery: mask entire left side.
                left_had_garbage = true;
                cols_to_mask.extend(data_start..star);
                if debug >= 1 {
                    stop_log.push((
                        header.clone(),
                        star,
                        "left_full".to_string(),
                        data_start,
                    ));
                }
            } else if recovery_left < star && recovery_left + 1 < star {
                left_had_garbage = true;
                cols_to_mask.extend((recovery_left + 1)..star);
                if debug >= 1 {
                    stop_log.push((
                        header.clone(),
                        star,
                        "left".to_string(),
                        recovery_left,
                    ));
                }
            } else {
                left_had_garbage = false;
            }

            // If garbage was found on exactly one side, the stop is a
            // frameshift artifact.  Toss the smaller surviving fragment.
            if left_had_garbage != right_had_garbage {
                let left_residues = (data_start..star)
                    .filter(|&i| seq[i] != b'-' && seq[i] != b'*' && !cols_to_mask.contains(&i))
                    .count();
                let right_residues = (star + 1..data_end)
                    .filter(|&i| seq[i] != b'-' && seq[i] != b'*' && !cols_to_mask.contains(&i))
                    .count();

                if left_residues < right_residues && left_residues > 0 {
                    cols_to_mask.extend(data_start..star);
                    if debug >= 1 {
                        stop_log.push((
                            header.clone(),
                            star,
                            "left_orphan".to_string(),
                            left_residues,
                        ));
                    }
                } else if right_residues < left_residues && right_residues > 0 {
                    cols_to_mask.extend((star + 1)..data_end);
                    if debug >= 1 {
                        stop_log.push((
                            header.clone(),
                            star,
                            "right_orphan".to_string(),
                            right_residues,
                        ));
                    }
                }
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
            intron_columns
                .entry(header.clone())
                .or_default()
                .extend(&cols_to_mask);
            result.push((header.clone(), masked));
        }
    }

    // Debug log.
    if !stop_log.is_empty() && debug >= 1 {
        if let Some(ld) = log_dir {
            let log_path = Path::new(ld).join("stop_side_trim_debug.csv");
            let write_header = !log_path.exists();
            if let Ok(mut f) = std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(&log_path)
            {
                if write_header {
                    let _ = writeln!(
                        f,
                        "header,stop_col,side,recovery_col"
                    );
                }
                for (h, stop_col, side, recovery_col) in &stop_log {
                    let safe = if h.contains(',') {
                        format!("\"{}\"", h)
                    } else {
                        h.clone()
                    };
                    let _ = writeln!(
                        f,
                        "{},{},{},{}",
                        safe, stop_col, side, recovery_col
                    );
                }
            }
        }
    }

    result
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
    stop_codon_side_threshold = 0.50,
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
    stop_codon_side_threshold: f64,
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

    // Detect and mask retained introns.
    let (after_intron, mut intron_masked_cols) =
        mask_retained_introns(&normalized, ref_suffix, &nt_seqs, min_ref_supported_gap, debug, &log_dir);

    // Mask low-quality sides of internal stop codons.
    let after_stop_sides = if stop_codon_side_threshold > 0.0 {
        mask_stop_sides(
            &after_intron,
            &mut intron_masked_cols,
            ref_suffix,
            stop_codon_side_threshold,
            8, // recovery_window
            debug,
            &log_dir,
        )
    } else {
        after_intron
    };

    // Gap cull: remove all-gap columns. Use Vec<bool> instead of HashSet.
    let mut gap_cull_mask = vec![false; aln_len];
    let mut gap_cull_count = 0usize;
    if gap_cull_threshold > 0.0 {
        let mut col_non_gap = vec![0usize; aln_len];
        for (_, seq) in &after_stop_sides {
            for (i, &ch) in seq.iter().enumerate() {
                if ch != b'-' {
                    col_non_gap[i] += 1;
                }
            }
        }
        for i in 0..aln_len {
            if col_non_gap[i] == 0 {
                gap_cull_mask[i] = true;
                gap_cull_count += 1;
            }
        }
    }

    let after_gap_cull: Vec<(String, Vec<u8>)>;
    let orig_index_map: Option<Vec<usize>>;

    if gap_cull_count > 0 {
        let oimap: Vec<usize> = (0..aln_len)
            .filter(|&i| !gap_cull_mask[i])
            .collect();
        after_gap_cull = after_stop_sides
            .iter()
            .map(|(h, s)| {
                let new_seq: Vec<u8> = oimap.iter().map(|&i| s[i]).collect();
                (h.clone(), new_seq)
            })
            .collect();
        aln_len = oimap.len();
        orig_index_map = Some(oimap);
    } else {
        after_gap_cull = after_stop_sides;
        orig_index_map = None;
    }

    let current_to_orig: Vec<usize> = match &orig_index_map {
        Some(m) => m.clone(),
        None => (0..aln_len).collect(),
    };

    let ref_seqs: Vec<&[u8]> = after_gap_cull
        .iter()
        .filter(|(h, _)| h.ends_with(ref_suffix))
        .map(|(_, s)| s.as_slice())
        .collect();

    // Edge columns as Vec<bool>.
    let mut edge_mask = vec![false; aln_len];

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
                for i in 0..left {
                    edge_mask[i] = true;
                }
                for i in (right + 1)..aln_len {
                    edge_mask[i] = true;
                }
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

    let mut masked_records: Vec<(String, Vec<u8>)> = Vec::new();
    // Store trim masks as Vec<bool> instead of HashSet<usize>.
    let mut trim_masks: Vec<Vec<bool>> = Vec::new();

    // Process each record. Single combined mask per record instead of two mask_columns calls.
    for (idx, (header, _)) in records.iter().enumerate() {
        let seq = &after_gap_cull[idx].1;

        // Start with edge mask as base.
        let mut combined_mask = edge_mask.clone();

        if !header.ends_with(ref_suffix) && !ref_seqs.is_empty() {
            // Apply edge mask to find data boundaries.
            let first_data = seq.iter().enumerate().position(|(i, &c)| c != b'-' && !edge_mask[i]);
            let (start, end_data) = if let Some(fd) = first_data {
                let last_data = seq.iter().enumerate().rposition(|(i, &c)| c != b'-' && !edge_mask[i]).unwrap();
                (fd, last_data + 1)
            } else {
                (0usize, 0usize)
            };

            // Trim outside data region.
            for i in 0..start {
                combined_mask[i] = true;
            }
            for i in end_data..aln_len {
                combined_mask[i] = true;
            }

            // Build effective sequence for block picking (apply combined mask so far).
            let eb: Vec<u8> = seq.iter().enumerate().map(|(i, &c)| {
                if combined_mask[i] { b'-' } else { c }
            }).collect();

            let left_trim_point = pick_final_block_start(
                &eb,
                &ref_has_data,
                min_ref_supported_gap,
                anchor_run,
                tail_max_data,
            );

            let right_trim_point = trim_sparse_tail_end(
                &eb,
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

            for i in start..left_trim_point {
                combined_mask[i] = true;
            }
            for i in right_trim_point..aln_len {
                combined_mask[i] = true;
            }
        }

        // Apply combined mask in one pass.
        let mut masked_seq = seq.clone();
        mask_columns_bool(&mut masked_seq, &combined_mask);

        // Min length filter.
        if !header.ends_with(ref_suffix) && min_seq_len > 0 {
            let non_gap = masked_seq.iter().filter(|&&c| c != b'-').count();
            if non_gap < min_seq_len {
                continue;
            }
        }

        masked_records.push((header.clone(), masked_seq));
        trim_masks.push(combined_mask);
    }

    if masked_records.is_empty() {
        return Ok((Vec::new(), HashMap::new(), HashMap::new(), HashMap::new(), HashMap::new()));
    }

    // Global empty column removal.
    let masked_len = masked_records[0].1.len();
    let mut masked_non_gap = vec![0usize; masked_len];
    for (_, seq) in &masked_records {
        for (i, &ch) in seq.iter().enumerate() {
            if ch != b'-' {
                masked_non_gap[i] += 1;
            }
        }
    }
    let mut global_empty_mask = vec![false; masked_len];
    let mut has_global_empty = false;
    for i in 0..masked_len {
        if masked_non_gap[i] == 0 {
            global_empty_mask[i] = true;
            has_global_empty = true;
        }
    }

    let trimmed_records: Vec<(String, String)> = if has_global_empty {
        let keep_indices: Vec<usize> = (0..masked_len)
            .filter(|&i| !global_empty_mask[i])
            .collect();
        masked_records
            .iter()
            .map(|(h, s)| {
                let new_seq: String = keep_indices.iter().map(|&i| s[i] as char).collect();
                (h.clone(), new_seq)
            })
            .collect()
    } else {
        masked_records
            .iter()
            .map(|(h, s)| {
                let new_seq: String = s.iter().map(|&c| c as char).collect();
                (h.clone(), new_seq)
            })
            .collect()
    };

    let global_empty_original: Vec<bool> = {
        let mut v = vec![false; orig_aln_len];
        for (i, &is_empty) in global_empty_mask.iter().enumerate() {
            if is_empty && i < current_to_orig.len() {
                v[current_to_orig[i]] = true;
            }
        }
        v
    };

    // Build nt_trim_map, gff_trim_map.
    let mut nt_trim_map: HashMap<String, HashSet<usize>> = HashMap::new();
    let mut gff_trim_map: HashMap<String, (usize, usize)> = HashMap::new();

    for (rec_i, (header, _)) in masked_records.iter().enumerate() {
        if header.ends_with(ref_suffix) {
            continue;
        }

        let trim_mask = &trim_masks[rec_i];

        // Build removed_original as Vec<bool> over orig_aln_len.
        let mut removed_original = vec![false; orig_aln_len];

        // gap_cull_mask columns
        for i in 0..orig_aln_len {
            if gap_cull_mask[i] {
                removed_original[i] = true;
            }
        }

        // trim mask mapped back through current_to_orig
        for (i, &is_trimmed) in trim_mask.iter().enumerate() {
            if is_trimmed && i < current_to_orig.len() {
                removed_original[current_to_orig[i]] = true;
            }
        }

        // global empty
        for i in 0..orig_aln_len {
            if global_empty_original[i] {
                removed_original[i] = true;
            }
        }

        // intron masked cols
        let intron_cols = intron_masked_cols.get(header.as_str());
        if let Some(ic) = intron_cols {
            for &col in ic {
                if col < orig_aln_len {
                    removed_original[col] = true;
                }
            }
        }

        let orig_seq = orig_seq_map.get(header.as_str()).copied().unwrap_or(&[] as &[u8]);
        let drops = codon_drops_from_mask(orig_seq, &removed_original);
        if !drops.is_empty() {
            nt_trim_map.insert(header.clone(), drops);
        }

        // edge_removed = removed_original minus intron cols
        let mut edge_removed = removed_original.clone();
        if let Some(ic) = intron_cols {
            for &col in ic {
                if col < orig_aln_len {
                    edge_removed[col] = false;
                }
            }
        }
        if let Some(trims) = count_edge_residue_trims_from_removed(orig_seq, &edge_removed) {
            if trims.0 > 0 || trims.1 > 0 {
                gff_trim_map.insert(header.clone(), trims);
            }
        }
    }

    // Build gff_intron_map and exon_split_map in one pass over intron_masked_cols.
    let mut gff_intron_map: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    let mut exon_split_map: HashMap<String, Vec<(usize, usize)>> = HashMap::new();

    for (header, cols) in &intron_masked_cols {
        let orig_seq = match orig_seq_map.get(header.as_str()) {
            Some(s) => *s,
            None => continue,
        };

        // Compute blocks once, reuse for both maps.
        let blocks = group_contiguous_blocks(cols);

        // gff_intron_map
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

        // exon_split_map
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
