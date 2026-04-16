mod aligner;
mod column_cull;
mod consensus;
mod dedupe;
mod exon_dp;
mod flexcull;
mod identity;
mod interval_tree;
mod overlap;
mod translate;

use bio::alignment::distance::simd::hamming;
use dedupe::fast_dedupe as rust_fast_dedupe;
use flexcull::*;
use overlap::get_overlap;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::collections::HashSet;
use std::path::PathBuf;

#[pyfunction]
fn bio_revcomp(sequence: String) -> String {
    String::from_utf8(bio::alphabets::dna::revcomp(sequence.into_bytes())).unwrap()
}

fn find_indices(sequence: &[u8], gap: u8) -> (usize, usize) {
    (
        sequence.iter().position(|c: &u8| *c != gap).unwrap_or(0),
        sequence.iter().rposition(|c: &u8| *c != gap).unwrap_or(0) + 1,
    )
}

#[pyfunction]
fn find_index_pair(sequence: &str, gap: &str) -> (usize, usize) {
    if gap.len() != 1 {
        panic!("Gap character must be a single character");
    }
    let gap_char = gap.as_bytes()[0];
    find_indices(sequence.as_bytes(), gap_char)
}

#[pyfunction]
fn is_same_kmer(one: &str, two: &str) -> bool {
    hamming(one.as_ref(), two.as_ref()) == 0
}

#[pyfunction]
fn constrained_distance(consensus: &str, candidate: &str) -> u64 {
    let can = candidate.as_bytes();
    let con = consensus.as_bytes();

    let (start, end) = find_indices(candidate.as_bytes(), b'-');
    let con = &con[start..end];
    let can = &can[start..end];
    let hamming_distance = hamming(con, can);
    let blank_count = con.iter().filter(|c: &&u8| **c == b'X').count() as u64;
    match hamming_distance > blank_count {
        true => hamming_distance - blank_count,
        false => 0,
    }
}

#[pyfunction]
fn blosum62_distance(one: String, two: String) -> f64 {
    // only used for ref to ref
    let first: &[u8] = one.as_bytes();
    let second: &[u8] = two.as_bytes();
    let mut score = 0;
    let mut max_first = 0;
    let mut max_second = 0;
    let length = first.len();
    let allowed: HashSet<u8> = HashSet::from([
        65, 84, 67, 71, 73, 68, 82, 80, 87, 77, 69, 81, 83, 72, 86, 76, 75, 70, 89, 78, 88, 90, 74,
        66, 79, 85, 42,
    ]);
    for i in 0..length {
        let mut char1 = first[i];
        let mut char2 = second[i];
        if first[i] == 45 {
            char1 = b'*';
        }
        if second[i] == 45 {
            char2 = b'*';
        }
        if !(allowed.contains(&char1)) {
            panic!("first[i]  {} not in allowed\n{}", char1 as char, one);
        }
        if !(allowed.contains(&char2)) {
            panic!("second[i] {} not in allowed\n{}", char2 as char, two);
        }
        // score += bio::scores::blosum62(char1, char2);
        score += bio::scores::blosum62(char1, char2);
        max_first += bio::scores::blosum62(char1, char1);
        max_second += bio::scores::blosum62(char2, char2);
    }
    let maximum_score = std::cmp::max(max_first, max_second);
    1.0_f64 - (score as f64 / maximum_score as f64)
}

#[pyfunction]
fn blosum62_candidate_to_reference(candidate: &str, reference: &str) -> f64 {
    let cand_bytes: &[u8] = candidate.as_bytes();
    let ref_bytes: &[u8] = reference.as_bytes();
    let mut score = 0;
    let mut max_first = 0;
    let mut max_second = 0;
    let length = cand_bytes.len();
    assert!(length == ref_bytes.len());
    let allowed: HashSet<u8> = HashSet::from([
        65, 84, 67, 71, 73, 68, 82, 80, 87, 77, 69, 81, 83, 72, 86, 76, 75, 70, 89, 78, 88, 90, 74,
        66, 79, 85, 42,
    ]);
    for i in 0..length {
        let mut char1 = cand_bytes[i];
        let mut char2 = ref_bytes[i];
        if cand_bytes[i] == 45 {
            char1 = b'*';
        }
        if ref_bytes[i] == 45 || ref_bytes[i] == b'*' {
            continue;
        }
        if !(allowed.contains(&char1)) {
            panic!("first[i]  {} not in allowed\n{}", char1 as char, candidate);
        }
        if !(allowed.contains(&char2)) {
            panic!("second[i] {} not in allowed\n{}", char2 as char, reference);
        }

        if char1 == b'*' || char1 == b'-' {
            score += -11;
            max_first += -11;
            max_second += bio::scores::blosum62(char2, char2);
        } else {
            score += bio::scores::blosum62(char1, char2);
            max_first += bio::scores::blosum62(char1, char1);
            max_second += bio::scores::blosum62(char2, char2);
        }
    }
    let maximum_score = std::cmp::max(max_first, max_second);
    1.0_f64 - (score as f64 / maximum_score as f64)
}

/// Column-filter helper: drops alignment columns where every sequence holds
/// a gap.  Operates on `(header, sequence)` pairs the SAPPHYRE Python side
/// already has on hand.
#[pyfunction]
fn delete_empty_columns_pairs(records: Vec<(String, String)>) -> Vec<(String, String)> {
    if records.is_empty() {
        return Vec::new();
    }
    let len = records[0].1.len();
    let n_recs = records.len();

    // Mark each column as kept the moment we see a non-gap character.
    let mut keep = vec![false; len];
    let mut remaining = len;
    'cols: for col in 0..len {
        for (_, seq) in &records {
            // Defensive: ragged inputs shouldn't happen but skip past the end.
            if col >= seq.len() {
                continue;
            }
            if seq.as_bytes()[col] != b'-' {
                keep[col] = true;
                remaining -= 1;
                if remaining == 0 {
                    break 'cols;
                }
                break;
            }
        }
    }

    let kept_indices: Vec<usize> = keep
        .iter()
        .enumerate()
        .filter_map(|(i, &k)| if k { Some(i) } else { None })
        .collect();

    let mut out = Vec::with_capacity(n_recs);
    let mut buf = Vec::with_capacity(kept_indices.len());
    for (header, seq) in records {
        let bytes = seq.as_bytes();
        buf.clear();
        for &i in &kept_indices {
            if i < bytes.len() {
                buf.push(bytes[i]);
            }
        }
        out.push((header, String::from_utf8(buf.clone()).unwrap()));
    }
    out
}

/// Sliding-window low-complexity check on a nucleotide sequence.
///
/// Mirrors `sapphyre.utils.low_complexity.is_low_complexity_nt` but runs the
/// inner counting / entropy loops in Rust.  Returns `true` if any window in
/// any non-N segment falls below either entropy threshold.
#[pyfunction]
#[pyo3(signature = (
    seq,
    window_size = 60,
    step = 15,
    dinuc_entropy_threshold = 2.5,
    trinuc_entropy_threshold = 3.75,
    min_segment_length = 30,
))]
fn is_low_complexity_nt(
    seq: &str,
    window_size: usize,
    step: usize,
    dinuc_entropy_threshold: f64,
    trinuc_entropy_threshold: f64,
    min_segment_length: usize,
) -> bool {
    // Encode bases A/C/G/T -> 0..3, anything else -> 4 ("N").
    let mut codes = Vec::with_capacity(seq.len());
    for b in seq.as_bytes() {
        let c = match b.to_ascii_uppercase() {
            b'A' => 0u8,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        };
        codes.push(c);
    }

    // Collect inclusive (start, end) spans of contiguous non-N bases.
    let mut segments: Vec<(usize, usize)> = Vec::new();
    let n = codes.len();
    let mut i = 0;
    while i < n {
        if codes[i] == 4 {
            i += 1;
            continue;
        }
        let mut j = i;
        while j < n && codes[j] != 4 {
            j += 1;
        }
        segments.push((i, j - 1));
        i = j;
    }
    if segments.is_empty() {
        return false;
    }

    let entropy_from_counts = |counts: &[u32], total: u32| -> f64 {
        if total == 0 {
            return 0.0;
        }
        let inv_total = 1.0 / total as f64;
        let log_total = (total as f64).ln();
        let mut sum_clogc = 0.0;
        for &c in counts {
            if c > 0 {
                let cf = c as f64;
                sum_clogc += cf * cf.ln();
            }
        }
        // H = sum(-(c/T) log2(c/T)) = (log T - (1/T) sum(c log c)) / ln 2
        (log_total - sum_clogc * inv_total) / std::f64::consts::LN_2
    };

    for (seg_start, seg_end) in segments {
        let seg_len = seg_end - seg_start + 1;
        if seg_len < min_segment_length {
            continue;
        }

        // Helper that evaluates one window over the segment.
        let mut check_window = |wstart: usize, wend: usize| -> bool {
            let abs_start = seg_start + wstart;
            let abs_end = seg_start + wend; // inclusive

            // Dinucleotides: positions abs_start..abs_end (i to i+1)
            let mut di = [0u32; 16];
            let mut di_total = 0u32;
            for k in abs_start..abs_end {
                let a = codes[k];
                let b = codes[k + 1];
                if a < 4 && b < 4 {
                    di[(a as usize) * 4 + b as usize] += 1;
                    di_total += 1;
                }
            }
            let dinuc_h = if di_total == 0 {
                4.0
            } else {
                entropy_from_counts(&di, di_total)
            };
            if dinuc_h < dinuc_entropy_threshold {
                return true;
            }

            // Trinucleotides: positions abs_start..abs_end-1 (i to i+2)
            if abs_end >= abs_start + 1 {
                let mut tri = [0u32; 64];
                let mut tri_total = 0u32;
                let stop = abs_end.saturating_sub(1);
                for k in abs_start..stop {
                    let a = codes[k];
                    let b = codes[k + 1];
                    let c = codes[k + 2];
                    if a < 4 && b < 4 && c < 4 {
                        tri[(a as usize) * 16 + (b as usize) * 4 + c as usize] += 1;
                        tri_total += 1;
                    }
                }
                let trinuc_h = if tri_total == 0 {
                    6.0
                } else {
                    entropy_from_counts(&tri, tri_total)
                };
                if trinuc_h < trinuc_entropy_threshold {
                    return true;
                }
            }

            false
        };

        // Sliding windows over the segment.
        if window_size >= seg_len {
            if check_window(0, seg_len - 1) {
                return true;
            }
            continue;
        }

        let mut w = 0usize;
        while w + window_size <= seg_len {
            if check_window(w, w + window_size - 1) {
                return true;
            }
            w += step;
        }
        // Tail window aligned to the end if it isn't already covered.
        let last_start = seg_len - window_size;
        if last_start > w.saturating_sub(step) {
            if check_window(last_start, seg_len - 1) {
                return true;
            }
        }
    }

    false
}

fn find_question_marks(sequences: Vec<&[u8]>, start: usize, stop: usize) -> HashSet<usize> {
    fn check_index(index: usize, sequences: &Vec<&[u8]>) -> bool {
        for seq in sequences {
            if seq[index] != b'-' {
                return false;
            }
        }
        true
    }

    let mut output = HashSet::new();
    for i in start..stop {
        if check_index(i, &sequences) {
            output.insert(i);
        }
    }
    output
}

fn replace_small_gaps(mut consensus_array: Vec<u8>, limit: usize) -> Vec<u8> {
    let mut last_seen = Option::None;
    for i in 0..consensus_array.len() {
        let bp = consensus_array[i];
        if bp == 63 {
            match last_seen {
                None => last_seen = Some(i),
                _ => {}
            }
        } else {
            match last_seen {
                Some(index) => {
                    if i - index < limit {
                        for letter in index..i {
                            consensus_array[letter] = 45;
                        }
                    }
                }
                _ => {}
            }
            last_seen = None;
        }
    }
    return consensus_array;
}

#[pyfunction]
fn convert_consensus(sequences: Vec<String>, consensus: &str) -> String {
    let replacement = b'?';
    let seqs: Vec<&[u8]> = sequences.iter().map(|seq| seq.as_bytes()).collect();
    let con_list = consensus.as_bytes();
    let (start, stop) = find_indices(con_list, b'X');
    let mut con_list = con_list.to_vec();
    let to_replace = find_question_marks(seqs, start, stop);
    for index in to_replace {
        con_list[index] = replacement;
    }
    con_list = replace_small_gaps(con_list, 8_usize);
    String::from_utf8(con_list).unwrap()
}

#[pyfunction]
pub fn preprocess_n_chunks(
    // Accepting String here fixes the "Can't extract str to Vec" error
    data: Vec<(String, String)>,
    min_length: usize,
) -> PyResult<Vec<(String, String)>> {
    // Pre-allocate for performance
    let mut results = Vec::with_capacity(data.len());

    for (header, seq) in data {
        // Fast fail for length
        if seq.len() < min_length {
            continue;
        }

        // Optimization: No 'N' found
        // seq.as_bytes() allows us to use SIMD scanning on the string data
        if !seq.as_bytes().contains(&b'N') && !seq.as_bytes().contains(&b'n') {
            // Move the strings into the results (No new allocation for the data)
            results.push((header, seq));
        } else {
            // Dirty path: Split by 'N' or 'n'
            // We use .as_bytes().split() because it's faster than string splitting
            for chunk_bytes in seq.as_bytes().split(|&b| b == b'N' || b == b'n') {
                if chunk_bytes.len() >= min_length {
                    // Convert the valid slice back to an owned String
                    let chunk_str = String::from_utf8_lossy(chunk_bytes).into_owned();
                    // We must clone the header for each sub-chunk
                    results.push((header.clone(), chunk_str));
                }
            }
        }
    }

    Ok(results)
}

#[pyfunction]
fn fast_dedupe(
    inputs: Vec<String>, // Accept a list of strings from Python
    out: &str,
    sort_by_size: bool,
    min_size: u64,
) -> PyResult<()> {
    // Convert all input strings into PathBufs
    let input_paths: Vec<PathBuf> = inputs.into_iter().map(PathBuf::from).collect();

    let out_path = PathBuf::from(out);

    // Call the updated Rust function
    rust_fast_dedupe(input_paths, out_path, sort_by_size, min_size)
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    Ok(())
}
// A Python module implemented in Rust.
#[pymodule]

fn sapphyre_tools(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(preprocess_n_chunks, m)?)?;
    m.add_function(wrap_pyfunction!(fast_dedupe, m)?)?;
    m.add_function(wrap_pyfunction!(blosum62_distance, m)?)?;
    m.add_function(wrap_pyfunction!(bio_revcomp, m)?)?;
    m.add_function(wrap_pyfunction!(constrained_distance, m)?)?;
    m.add_function(wrap_pyfunction!(convert_consensus, m)?)?;
    m.add_function(wrap_pyfunction!(consensus::dumb_consensus, m)?)?;
    m.add_function(wrap_pyfunction!(consensus::dumb_consensus_dupe, m)?)?;
    m.add_function(wrap_pyfunction!(consensus::consensus_distance, m)?)?;
    m.add_function(wrap_pyfunction!(find_index_pair, m)?)?;
    m.add_function(wrap_pyfunction!(blosum62_candidate_to_reference, m)?)?;
    m.add_function(wrap_pyfunction!(identity::filter_nt, m)?)?;
    m.add_function(wrap_pyfunction!(delete_empty_columns_pairs, m)?)?;
    m.add_function(wrap_pyfunction!(is_low_complexity_nt, m)?)?;
    m.add_function(wrap_pyfunction!(exon_dp::exon_dp, m)?)?;
    m.add_function(wrap_pyfunction!(join_by_tripled_index, m)?)?;
    m.add_function(wrap_pyfunction!(join_with_exclusions, m)?)?;
    m.add_function(wrap_pyfunction!(join_triplets_with_exclusions, m)?)?;
    m.add_function(wrap_pyfunction!(get_overlap, m)?)?;
    m.add_function(wrap_pyfunction!(is_same_kmer, m)?)?;
    m.add_function(wrap_pyfunction!(translate::translate, m)?)?;
    m.add_function(wrap_pyfunction!(interval_tree::del_cols, m)?)?;
    m.add_function(wrap_pyfunction!(aligner::hmm_align, m)?)?;
    m.add_function(wrap_pyfunction!(column_cull::cull_columns, m)?)?;
    m.add_function(wrap_pyfunction!(column_cull::apply_gff_culls, m)?)?;

    m.add_class::<interval_tree::OverlapTree>()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::blosum62_distance;

    #[test]
    fn test_blosum62_gap_penalty() {
        let result = blosum62_distance(String::from("A"), String::from("-"));
        assert_eq!(result, 2.0)
    }
}
