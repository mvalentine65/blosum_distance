mod align;
mod consensus;
mod entropy;
mod flexcull;
mod hit;
mod identity;
mod interval_tree;
mod motif;
mod overlap;
mod sigclust;
mod translate;
mod utils;
mod dedupe;

use dedupe::fast_dedupe as rust_fast_dedupe;
use bio::alignment::distance::simd::hamming;
use flexcull::*;
use hit::Hit;
use hit::ReferenceHit;
use itertools::enumerate;
use overlap::{get_overlap, get_overlap_percent};
use pyo3::prelude::*;
use std::collections::HashSet;
use std::path::PathBuf;
use pyo3::exceptions::PyRuntimeError;

#[pyfunction]
fn asm_index_split(sequence: String) -> Vec<(usize, usize)> {
    let mut start = 0_usize;
    let mut in_data = false;
    let mut gap_count = 0_usize;
    let mut output = Vec::new();

    for (i, bp) in enumerate(sequence.as_bytes()) {
        if *bp == b'-' {
            gap_count += 1;
            if gap_count >= 20 {
                if in_data {
                    output.push((start, i - gap_count + 1))
                }
                in_data = false
            }
        } else {
            gap_count = 0;
            if !in_data {
                in_data = true;
                start = i
            }
        }
    }
    if in_data {
        output.push((start, sequence.len() - gap_count))
    }
    output
}

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

fn find_character_indices(sequence: &[u8], character: u8) -> Option<(usize, usize)> {
    match (
        sequence.iter().position(|c: &u8| *c == character),
        sequence.iter().rposition(|c: &u8| *c == character),
    ) {
        (Some(first), Some(last)) => Some((first, last + 1)),
        _ => None,
    }
}
// Searches the given sequence for the first and last occurrence of the given character.
// Returns an Option[tuple[int, int]].
// If found, returns tuple with the indices. Start is inclusive, stop is exclusive. [start, stop)
// If the character is not found, returns None.
// Requires ascii input for both arguments. Character must be a single ascii character. (length of 1)
// Search is case sensitive.
#[pyfunction]
fn find_first_last_character(sequence: &str, character: &str) -> Option<(usize, usize)> {
    if character.len() != 1 {
        panic!("Search character must be a single ascii character.");
    }
    find_character_indices(sequence.as_bytes(), character.as_bytes()[0])
}

#[pyfunction]
fn has_data(sequence: &str, gap: char) -> bool {
    sequence.as_bytes().iter().any(|x| *x != gap as u8)
}

#[pyfunction]
fn simd_hamming(one: &str, two: &str) -> u64 {
    return hamming(one.as_ref(), two.as_ref());
}

#[pyfunction]
fn is_same_kmer(one: &str, two: &str) -> bool {
    hamming(one.as_ref(), two.as_ref()) == 0
}

#[pyfunction]
fn score_splits(ref_slice: &str, splits: Vec<(String, u64)>) -> u64 {
    let mut min_distance: u64 = ref_slice.len() as u64;
    let mut winning_index: u64 = 0;
    let mut current: u64 = ref_slice.len() as u64;
    for (seq, index) in splits {
        current = simd_hamming(ref_slice, seq.as_ref());
        if current < min_distance {
            min_distance = current;
            winning_index = index;
        }
    }
    winning_index
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
fn len_without_gaps(x: String) -> usize {
    return x.as_bytes().iter().filter(|&&c| c != b'-').count();
}

fn _excise_consensus_tail(consensus: &String, limit: f64) -> (String, usize) {
    let step = 16_usize;
    let mut percent_invalid = 1.0_f64;
    let mut num_invalid = 0;
    let mut denominator = 0;
    let length = consensus.len();
    let mut end = length;
    let mut start = length - step;
    while start >= step && percent_invalid > limit {
        num_invalid = consensus[start..end].bytes().filter(|&x| x == b'X').count();
        denominator = step;
        percent_invalid = num_invalid as f64 / denominator as f64;
        if percent_invalid < limit {
            break;
        }
        start -= step;
        end -= step;
    }
    while end > 0 && end < usize::MAX && percent_invalid > limit {
        end -= 1;
        if &consensus[end..end + 1] == "X" {
            num_invalid += 1;
        }
        denominator += 1;
        percent_invalid = num_invalid as f64 / denominator as f64;
    }
    (consensus[0..end].to_string(), end)
}   

#[pyfunction]
fn excise_consensus_tail(consensus: String, limit: f64) -> (String, usize) {
    _excise_consensus_tail(&consensus, limit)
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

#[pyfunction]
pub fn debug_blosum62_candidate_to_reference(
    candidate: &str,
    reference: &str,
) -> (
    f64,
    i32,
    Vec<i32>,
    Vec<i32>,
    Vec<i32>,
    Vec<i32>,
    Vec<i32>,
    Vec<i32>,
) {
    let mut cand_locations = Vec::with_capacity(candidate.len());
    let mut cand_sums = Vec::with_capacity(candidate.len());

    let mut ref_locations = Vec::with_capacity(reference.len());
    let mut ref_sums = Vec::with_capacity(reference.len());

    let mut score_locations = Vec::with_capacity(candidate.len());
    let mut score_sums = Vec::with_capacity(candidate.len());

    let cand_bytes: &[u8] = candidate.as_bytes();
    let ref_bytes: &[u8] = reference.as_bytes();

    let mut this_score = 0;
    let mut this_max_first = 0;
    let mut this_max_second = 0;
    let mut score = 0;
    let mut max_first = 0;
    let mut max_second = 0;
    let length = cand_bytes.len();
    let allowed: HashSet<u8> = HashSet::from([
        65, 84, 67, 71, 73, 68, 82, 80, 87, 77, 69, 81, 83, 72, 86, 76, 75, 70, 89, 78, 88, 90, 74,
        66, 79, 85, 42,
    ]);
    for i in 0..length {
        let mut char1 = cand_bytes[i];
        let char2 = ref_bytes[i];
        if cand_bytes[i] == 45 {
            char1 = b'*';
        }
        if ref_bytes[i] == 45 {
            cand_locations.push(0);
            cand_sums.push(max_first);

            ref_locations.push(0);
            ref_sums.push(max_second);

            score_locations.push(0);
            score_sums.push(score);
            continue;
        }
        if !(allowed.contains(&char1)) {
            panic!("first[i]  {} not in allowed\n{}", char1 as char, candidate);
        }
        if !(allowed.contains(&char2)) {
            panic!("second[i] {} not in allowed\n{}", char2 as char, reference);
        }
        if char1 == b'*' || char1 == b'-' {
            this_score = -11;
            this_max_first = -11;
            this_max_second += bio::scores::blosum62(char2, char2);
        } else {
            this_score += bio::scores::blosum62(char1, char2);
            this_max_first += bio::scores::blosum62(char1, char1);
            this_max_second += bio::scores::blosum62(char2, char2);
        }

        score += this_score;
        max_first += this_max_first;
        max_second += this_max_second;

        cand_locations.push(this_max_first);
        cand_sums.push(max_first);

        ref_locations.push(this_max_second);
        ref_sums.push(max_second);

        score_locations.push(this_score);
        score_sums.push(score);
    }
    let maximum_score = std::cmp::max(max_first, max_second);
    let end_result = 1.0_f64 - (score as f64 / maximum_score as f64);

    (
        end_result,
        maximum_score,
        cand_locations,
        cand_sums,
        ref_locations,
        ref_sums,
        score_locations,
        score_sums,
    )
}

#[pyfunction]
fn delete_empty_columns(raw_fed_sequences: Vec<String>) -> (Vec<String>, Vec<usize>) {
    let mut result = Vec::<String>::with_capacity(raw_fed_sequences.len());
    let mut headers = Vec::<&String>::with_capacity(raw_fed_sequences.len() / 2 + 1);
    let mut sequences = Vec::<&[u8]>::with_capacity(raw_fed_sequences.len() / 2 + 1);
    for i in (1..raw_fed_sequences.len()).step_by(2) {
        headers.push(&raw_fed_sequences[i - 1]);
        sequences.push(raw_fed_sequences[i].as_bytes())
    }
    if sequences.is_empty() {
        return (result, Vec::<usize>::new());
    }
    let mut positions_to_keep: Vec<usize> = Vec::with_capacity(sequences[0].len());
    for i in 0..sequences[0].len() {
        for sequence in &sequences {
            if sequence[i] != b'-' {
                positions_to_keep.push(i);
                break;
            }
        }
    }
    let mut temp = Vec::<u8>::with_capacity(sequences[0].len());
    for i in 0..sequences.len() {
        temp = positions_to_keep
            .iter()
            .map(|index| sequences[i][*index])
            .collect();
        result.push(headers[i].to_string());
        result.push(String::from_utf8(temp).unwrap());
    }
    (result, positions_to_keep)
}

fn _check_index(sequences: &Vec<&[u8]>, i: usize, target: u8) -> bool {
    for seq in sequences.iter() {
        if seq[i] != target {
            return false;
        }
    }
    true
}

fn _find_replaceble(sequences: Vec<&[u8]>, target: u8) -> HashSet<usize> {
    let mut output = HashSet::new();
    for i in 0..sequences[0].len() {
        if _check_index(&sequences, i, target) {
            output.insert(i);
        }
    }
    output
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
    inputs: Vec<String>,      // Accept a list of strings from Python
    out: &str,
    sort_by_size: bool,
    min_size: u64,
) -> PyResult<()> {
    // Convert all input strings into PathBufs
    let input_paths: Vec<PathBuf> = inputs
        .into_iter()
        .map(PathBuf::from)
        .collect();
        
    let out_path = PathBuf::from(out);

    // Call the updated Rust function
    rust_fast_dedupe(
        input_paths,
        out_path,
        sort_by_size,
        min_size,
    ).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    Ok(())
}
// A Python module implemented in Rust.
#[pymodule]

fn sapphyre_tools(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(preprocess_n_chunks, m)?)?;
    m.add_function(wrap_pyfunction!(fast_dedupe, m)?)?;
    m.add_function(wrap_pyfunction!(asm_index_split, m)?)?;
    m.add_function(wrap_pyfunction!(blosum62_distance, m)?)?;
    m.add_function(wrap_pyfunction!(bio_revcomp, m)?)?;
    m.add_function(wrap_pyfunction!(constrained_distance, m)?)?;
    m.add_function(wrap_pyfunction!(convert_consensus, m)?)?;
    m.add_function(wrap_pyfunction!(len_without_gaps, m)?)?;
    m.add_function(wrap_pyfunction!(consensus::dumb_consensus, m)?)?;
    m.add_function(wrap_pyfunction!(consensus::dumb_consensus_dupe, m)?)?;
    m.add_function(wrap_pyfunction!(consensus::filter_regions, m)?)?;
    m.add_function(wrap_pyfunction!(consensus::consensus_distance, m)?)?;
    m.add_function(wrap_pyfunction!(excise_consensus_tail, m)?)?;
    m.add_function(wrap_pyfunction!(find_index_pair, m)?)?;
    m.add_function(wrap_pyfunction!(find_first_last_character, m)?)?;
    m.add_function(wrap_pyfunction!(has_data, m)?)?;
    m.add_function(wrap_pyfunction!(blosum62_candidate_to_reference, m)?)?;
    m.add_function(wrap_pyfunction!(debug_blosum62_candidate_to_reference, m)?)?;
    m.add_function(wrap_pyfunction!(identity::filter_nt, m)?)?;
    m.add_function(wrap_pyfunction!(score_splits, m)?)?;
    m.add_function(wrap_pyfunction!(simd_hamming, m)?)?;
    m.add_function(wrap_pyfunction!(delete_empty_columns, m)?)?;
    m.add_function(wrap_pyfunction!(join_by_tripled_index, m)?)?;
    m.add_function(wrap_pyfunction!(join_with_exclusions, m)?)?;
    m.add_function(wrap_pyfunction!(join_triplets_with_exclusions, m)?)?;
    m.add_function(wrap_pyfunction!(get_overlap, m)?)?;
    m.add_function(wrap_pyfunction!(is_same_kmer, m)?)?;
    m.add_function(wrap_pyfunction!(get_overlap_percent, m)?)?;
    m.add_function(wrap_pyfunction!(align::make_aligned_ingredients, m)?)?;
    m.add_function(wrap_pyfunction!(align::run_intermediate, m)?)?;
    m.add_function(wrap_pyfunction!(align::process_clusters, m)?)?;
    m.add_function(wrap_pyfunction!(align::align_remove_empty_columns, m)?)?;
    m.add_function(wrap_pyfunction!(align::align_remove_dashes, m)?)?;
    m.add_function(wrap_pyfunction!(entropy::entropy_filter, m)?)?;
    m.add_function(wrap_pyfunction!(entropy::entropy, m)?)?;
    m.add_function(wrap_pyfunction!(sigclust::sigclust, m)?)?;
    m.add_function(wrap_pyfunction!(sigclust::sigclust_with_sequence, m)?)?;
    m.add_function(wrap_pyfunction!(utils::write_fasta_compressed, m)?)?;
    m.add_function(wrap_pyfunction!(utils::write_fasta_uncompressed, m)?)?;
    m.add_function(wrap_pyfunction!(translate::translate, m)?)?;
    m.add_function(wrap_pyfunction!(interval_tree::del_cols, m)?)?;

    m.add_class::<Hit>()?;
    m.add_class::<ReferenceHit>()?;
    m.add_class::<interval_tree::OverlapTree>()?;
    m.add_class::<motif::ScoredPosition>()?;
    m.add_class::<motif::ProteinMotif>()?;

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
