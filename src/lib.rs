mod hit;
extern crate bio;

extern crate pyo3;

use bio::alignment::distance::simd::hamming;
use bio::scores::blosum62;
use itertools::{enumerate, Itertools};
use hit::Hit;
use hit::ReferenceHit;
// use hit::hit_from_series;
use pyo3::prelude::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::{prelude::*, BufReader};


#[pyfunction]
fn find_references_and_candidates(path: String) -> (Vec<String>, Vec<String>) {
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let mut candidates = Vec::new();
    let mut references = Vec::new();
    let mut reached_candidates = false;
    let mut line: String;
    for line_in in reader.lines() {
        line = line_in.unwrap();
        if reached_candidates {
            candidates.push(line);
        } else if line.starts_with('>') && !line.ends_with('.') {
            reached_candidates = true;
            candidates.push(line)
        } else {
            references.push(line)
        }
    }
    (references, candidates)
}

#[pyfunction]
fn batch_reverse_complement(list: Vec<String>) -> Vec<String> {
    list.into_iter().map(bio_revcomp).collect()
}

#[pyfunction]
fn bio_revcomp(sequence: String) -> String {
    String::from_utf8(bio::alphabets::dna::revcomp(sequence.into_bytes())).unwrap()
}

// fn not_skip_character(character: u8) -> bool {
//     const HYPHEN: u8 = 45;
//     const ASTERISK: u8 = 42;
//     character != HYPHEN && character != ASTERISK
// }

fn blosum62_check(a: u8, b: u8) -> i32 {
    if a == 45 || b == 45 {
        return 0_i32;
    }
    blosum62(a, b)
}

fn blosum62_self_check(a: u8, b: u8) -> i32 {
    if a == 45 || b == 45 {
        return 0_i32;
    }
    blosum62(a, a)
}

fn self_distance_vectors(seq1: &str, seq2: &str) -> (Vec<i32>, Vec<i32>) {
    (
        seq1.as_bytes()
            .iter()
            .zip(seq2.as_bytes().iter())
            .scan(0_i32, |state, (a, b)| {
                *state += blosum62_self_check(*a, *b);
                Some(*state)
            })
            .collect(),
        seq2.as_bytes()
            .iter()
            .zip(seq1.as_bytes().iter())
            .scan(0_i32, |state, (a, b)| {
                *state += blosum62_self_check(*a, *b);
                Some(*state)
            })
            .collect(),
    )
}

#[pyfunction]
fn make_ref_distance_vector_blosum62(seq1: &str, seq2: &str) -> (Vec<i32>, (Vec<i32>, Vec<i32>)) {
    let numerator = seq1
        .as_bytes()
        .iter()
        .zip(seq2.as_bytes().iter())
        .scan(0, |state, (a, b)| {
            *state += blosum62_check(*a, *b);
            Some(*state)
        })
        .collect();
    (numerator, self_distance_vectors(seq1, seq2))
}

#[pyfunction]
fn make_ref_distance_matrix_supervector(refs: Vec<&str>) -> Vec<(Vec<i32>, (Vec<i32>, Vec<i32>))> {
    refs.iter()
        .combinations(2)
        // .iter()
        .map(|seq_vec| make_ref_distance_vector_blosum62(&seq_vec[0], &seq_vec[1]))
        .collect()
}

fn extract_distance(
    cross_distance: &Vec<i32>,
    max_vec1: &Vec<i32>,
    max_vec2: &Vec<i32>,
    start_index: usize,
    stop_index: usize,
) -> f64 {
    let denominator = match start_index {
        0 => std::cmp::max(max_vec1[stop_index - 1], max_vec2[stop_index - 1]),
        _ => std::cmp::max(
            max_vec1[stop_index - 1] - max_vec1[start_index - 1],
            max_vec2[stop_index - 1] - max_vec2[start_index - 1],
        ),
    };
    let numerator = match start_index {
        0 => cross_distance[stop_index - 1],
        _ => cross_distance[stop_index - 1] - cross_distance[start_index - 1],
    };
    1.0_f64 - (numerator as f64 / denominator as f64)
}

#[pyfunction]
fn ref_index_vector(
    supervector: Vec<(Vec<i32>, (Vec<i32>, Vec<i32>))>,
    start: usize,
    end: usize,
) -> Vec<f64> {
    supervector
        .iter()
        .map(|(cross, (max1, max2))| extract_distance(cross, max1, max2, start, end))
        .collect()
}

fn find_indices(sequence: &[u8], gap: u8) -> (usize, usize) {
    (sequence.iter().position(|c: &u8| *c != gap).unwrap_or(0), sequence.iter().rposition(|c: &u8| *c != gap).unwrap_or(0)+1)
}

#[pyfunction]
fn find_index_pair(sequence: &str, gap: char) -> (usize, usize) {
    find_indices(sequence.as_bytes(), gap as u8)
}

#[pyfunction]
fn constrained_distance_bytes(consensus: &[u8], candidate: &[u8]) -> u64 {
    let (start, end) = find_indices(candidate,b'-');
    let con_slice = &consensus[start..end];
    let can_slice = &candidate[start..end];
    hamming(con_slice, can_slice) - con_slice.iter().filter(|c: &&u8| **c == b'X').count() as u64
}


#[pyfunction]
fn has_data(sequence: &str, gap: char) -> bool{
    sequence.as_bytes().iter().any(|x| *x != gap as u8)
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
        false => 0
    }
}


#[pyfunction]
fn dumb_consensus(sequences: Vec<&str>, threshold: f64) -> String {
    let first = &sequences[0];
    let mut total_at_position = vec![0_u32;first.len()];
    let mut counts_at_position = vec![[0_u32;27];first.len()];
    const ASCII_OFFSET: u8= 65;
    const HYPHEN:u8= 45;
    const ASTERISK:u8= 42;
    let mut min = usize::MAX;
    let mut max:usize = 0;
    for sequence in sequences.iter() {
        let seq = sequence.as_bytes();
        let (start, end) = find_indices(seq, b'-');
        if start < min {min = start;}
        if end > max {max = end;}
        // let seq = &seq[..];
        for index in start..end {
            if index == seq.len() {continue;}
            total_at_position[index] += 1;
            if !(seq[index] == HYPHEN || seq[index] == ASTERISK) {
                // if seq[index]-ASCII_OFFSET == 233 {println!("{}",seq[index]);}
                counts_at_position[index][(seq[index] - ASCII_OFFSET) as usize] += 1;
            } else {
                counts_at_position[index][26] += 1;
            }
        }
    }
    let mut output = Vec::<u8>::with_capacity(total_at_position.len());
    for ((_, total), counts) in enumerate(total_at_position).zip(counts_at_position.iter()) {
        if total == 0 {
                output.push(b'X');
                continue;
        } // if no characters at position, continue
        let mut max_count:u32 = 0;
        let mut winner = b'X';  // default to X if no winner found
        for (index, count) in enumerate(counts) {
            if *count as f64 / total as f64 > threshold {
                if *count > max_count {
                    max_count = *count;
                    if index != 26 { winner = index as u8 +ASCII_OFFSET;}
                    else {winner = HYPHEN;}
                }
            }
        }
        output.push(winner);
    }
    String::from_utf8(output).unwrap()
}


#[pyfunction]
fn blosum62_distance(one: String, two: String) -> f64 {
    let first: &[u8] = one.as_bytes();
    let second: &[u8] = two.as_bytes();
    let mut score = 0;
    let mut max_first = 0;
    let mut max_second = 0;
    let length = first.len();
    let allowed: HashSet<u8> = HashSet::from([
        65, 84, 67, 71, 73, 68, 82, 80, 87, 77, 69, 81, 83, 72, 86, 76, 75, 70, 89, 78, 88, 90, 74,
        66, 79, 85,42
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
    let allowed: HashSet<u8> = HashSet::from([
        65, 84, 67, 71, 73, 68, 82, 80, 87, 77, 69, 81, 83, 72, 86, 76, 75, 70, 89, 78, 88, 90, 74,
        66, 79, 85,42
    ]);
    for i in 0..length {
        let mut char1 = cand_bytes[i];
        let char2 = ref_bytes[i];
        if cand_bytes[i] == 45 {
            char1 = b'*';
        }
        if ref_bytes[i] == 45 {
            continue;
        }
        if !(allowed.contains(&char1)) {
            panic!("first[i]  {} not in allowed\n{}", char1 as char, candidate);
        }
        if !(allowed.contains(&char2)) {
            panic!("second[i] {} not in allowed\n{}", char2 as char, reference);
        }
        score += bio::scores::blosum62(char1, char2);
        max_first += bio::scores::blosum62(char1, char1);
        max_second += bio::scores::blosum62(char2, char2);
    }
    let maximum_score = std::cmp::max(max_first, max_second);
    1.0_f64 - (score as f64 / maximum_score as f64)
}

// A Python module implemented in Rust.
#[pymodule]
fn phymmr_tools(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(blosum62_distance, m)?)?;
    m.add_function(wrap_pyfunction!(batch_reverse_complement, m)?)?;
    m.add_function(wrap_pyfunction!(bio_revcomp, m)?)?;
    m.add_function(wrap_pyfunction!(find_references_and_candidates, m)?)?;
    m.add_function(wrap_pyfunction!(constrained_distance, m)?)?;
    m.add_function(wrap_pyfunction!(constrained_distance_bytes, m)?)?;
    m.add_function(wrap_pyfunction!(dumb_consensus, m)?)?;
    m.add_function(wrap_pyfunction!(find_index_pair, m)?)?;
    m.add_function(wrap_pyfunction!(has_data, m)?)?;
    m.add_function(wrap_pyfunction!(blosum62_candidate_to_reference,m)?)?;
    m.add_class::<Hit>()?;
    m.add_class::<ReferenceHit>()?;
    // m.add_function(wrap_pyfunction!(hit_from_series, m)?)?;

    Ok(())
}


#[cfg(test)]
mod tests {
    use crate::blosum62_distance;

    #[test]
    fn test_blosum62_gap_penalty() {
        let result = blosum62_distance(String::from("A"),String::from("-"));
        assert_eq!(result, 2.0)
    }
}