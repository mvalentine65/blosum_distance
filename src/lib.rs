mod hit;
extern crate bio;

extern crate pyo3;

use bio::alignment::distance::simd::hamming;
use itertools::enumerate;
use hit::Hit;
use hit::ReferenceHit;
use pyo3::prelude::*;
use std::collections::HashSet;

#[pyfunction]
fn bio_revcomp(sequence: String) -> String {
    String::from_utf8(bio::alphabets::dna::revcomp(sequence.into_bytes())).unwrap()
}

fn find_indices(sequence: &[u8], gap: u8) -> (usize, usize) {
    (sequence.iter().position(|c: &u8| *c != gap).unwrap_or(0), sequence.iter().rposition(|c: &u8| *c != gap).unwrap_or(0)+1)
}

#[pyfunction]
fn find_index_pair(sequence: &str, gap: char) -> (usize, usize) {
    find_indices(sequence.as_bytes(), gap as u8)
}

#[pyfunction]
fn has_data(sequence: &str, gap: char) -> bool{
    sequence.as_bytes().iter().any(|x| *x != gap as u8)
}


#[pyfunction]
fn simd_hamming(one: &str, two: &str) -> u64 {
    return hamming(one.as_ref(), two.as_ref())
}

#[pyfunction]
fn score_splits(ref_slice: &str, splits: Vec<(&str, u64)>) -> u64 {
    let mut min_distance: u64 = ref_slice.len() as u64;
    let mut winning_index: u64 = 0;
    let mut current: u64 = ref_slice.len() as u64;
    for (seq, index) in splits {
        current = simd_hamming(ref_slice, seq);
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
fn dumb_consensus_dupe(sequences: Vec<(&str, u32)>, threshold: f64) -> String {
    let (first, _) = &sequences[0];
    let mut total_at_position = vec![0_u32;first.len()];
    let mut counts_at_position = vec![[0_u32;27];first.len()];
    const ASCII_OFFSET: u8= 65;
    const HYPHEN:u8= 45;
    const ASTERISK:u8= 42;
    let mut min = usize::MAX;
    let mut max:usize = 0;
    for (sequence, count) in sequences.iter() {
        let seq = sequence.as_bytes();
        let (start, end) = find_indices(seq, b'-');
        if start < min {min = start;}
        if end > max {max = end;}
        // let seq = &seq[..];
        for index in start..end {
            if index == seq.len() {continue;}
            total_at_position[index] += count;
            if !(seq[index] == HYPHEN || seq[index] == ASTERISK) {
                // if seq[index]-ASCII_OFFSET == 233 {println!("{}",seq[index]);}
                counts_at_position[index][(seq[index] - ASCII_OFFSET) as usize] += count;
            } else {
                counts_at_position[index][26] += count;
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
    m.add_function(wrap_pyfunction!(bio_revcomp, m)?)?;
    m.add_function(wrap_pyfunction!(constrained_distance, m)?)?;
    m.add_function(wrap_pyfunction!(dumb_consensus, m)?)?;
    m.add_function(wrap_pyfunction!(dumb_consensus_dupe, m)?)?;
    m.add_function(wrap_pyfunction!(find_index_pair, m)?)?;
    m.add_function(wrap_pyfunction!(has_data, m)?)?;
    m.add_function(wrap_pyfunction!(blosum62_candidate_to_reference,m)?)?;
    m.add_function(wrap_pyfunction!(score_splits,m)?)?;
    m.add_function(wrap_pyfunction!(simd_hamming,m)?)?;
    m.add_class::<Hit>()?;
    m.add_class::<ReferenceHit>()?;
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