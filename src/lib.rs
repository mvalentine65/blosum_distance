mod hit;
mod flexcull;
extern crate bio;

extern crate pyo3;

use bio::alignment::distance::simd::hamming;
use flexcull::*;
use hit::Hit;
use hit::ReferenceHit;
use itertools::enumerate;
use pyo3::prelude::*;
use std::collections::HashSet;


#[pyclass]
struct LocationData {
    winning_count: u32,
    unique_chars: u8,
    mean: f64,
}


#[pyfunction]
fn bio_revcomp(sequence: String) -> String {
    String::from_utf8(bio::alphabets::dna::revcomp(sequence.into_bytes())).unwrap()
}

fn find_indices(sequence: &[u8], gap: u8) -> (usize, usize) {
    (sequence.iter().position(|c: &u8| *c != gap).unwrap_or(0), sequence.iter().rposition(|c: &u8| *c != gap).unwrap_or(0)+1)
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
fn make_location_data(letters: &[u32;27]) -> LocationData {
    let averages: Vec<&u32> = letters.iter()
        .filter(|x| **x > 0)
        .collect();
    // let mean = averages.iter().map(|x| **x).sum::<u32>() as f64 / averages.len() as f64;
    LocationData {
        winning_count: **averages.iter().max().unwrap(),
        unique_chars: averages.len() as u8,
        mean: averages.iter().map(|x| **x).sum::<u32>() as f64 / averages.len() as f64,
    }
}

fn _excise_consensus_tail(consensus: &String, limit: f64) -> (String, usize) {
    let STEP = 16_usize;
    let mut percent_invalid = 1.0_f64;
    let mut num_invalid = 0;
    let mut denominator = 0;
    let length = consensus.len();
    let mut end = length;
    let mut start = length - STEP;
    while percent_invalid > limit {
        num_invalid += consensus[start..end].bytes().filter(|x| *x == b'X').count();
        denominator += STEP;
        percent_invalid = num_invalid as f64 / denominator as f64;

        if start == 0 && percent_invalid > limit {
            break;
        } // prevent underflow
        if start < 16 {
            start = 0;
        } else {
            start -= 16;
        }
        end -= 16;
    }
    (consensus[0..start].to_string(), length-start)
}

#[pyfunction]
fn excise_consensus_tail(consensus: String, limit: f64) -> (String, usize) {
    let STEP = 16_usize;
    let mut percent_invalid = 1.0_f64;
    let mut num_invalid = 0;
    let mut denominator = 0;
    let length = consensus.len();
    let mut end = length;
    let mut start = length - STEP;
    while percent_invalid > limit {
        num_invalid += consensus[start..end].bytes().filter(|x| *x == b'X').count();
        denominator += STEP;
        percent_invalid = num_invalid as f64 / denominator as f64;

        if start == 0 && percent_invalid > limit {
            break;
        } // prevent underflow
        if start < 16 {
            start = 0;
        } else {
            start -= 16;
        }
        end -= 16;
    }
    (consensus[0..start].to_string(), length-start)
}



#[pyfunction]
fn dumb_consensus_with_excise(sequences: Vec<&str>, consensus_threshold: f64, excise_threshold: f64) -> (String, usize, String) {
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
    // let mut ratios = Vec::<u8>::with_capacity(total_at_position.len());
    // let mut ratio = Vec::<u8>::with_capacity(to)
    for ((_, total), counts) in enumerate(total_at_position).zip(counts_at_position.iter()) {
        if total == 0 {
            output.push(b'X');
            continue;
        } // if no characters at position, continue
        let mut max_count:u32 = 0;
        let mut winner = b'X';  // default to X if no winner found
        for (index, count) in enumerate(counts) {
            if *count as f64 / total as f64 > consensus_threshold {
                if *count > max_count {
                    max_count = *count;
                    if index != 26 { winner = index as u8 +ASCII_OFFSET;}
                    else {winner = HYPHEN;}
                }
            }
        }

        output.push(winner);
    }
    // let locations: Vec<LocationData> = counts_at_position.iter()
    //     .map(|letters| weigh_winner(letters))
    //     .collect();
    let consensus = String::from_utf8(output).unwrap();
    let (excised, cut_length) = _excise_consensus_tail(&consensus, excise_threshold);
    (excised, cut_length, consensus)
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


#[pyfunction]
fn delete_empty_columns(raw_fed_sequences: Vec<&str>) -> (Vec<String>, Vec<usize>) {
    let mut result = Vec::<String>::with_capacity(raw_fed_sequences.len());
    let mut headers = Vec::<&str>::with_capacity(raw_fed_sequences.len() / 2 + 1);
    let mut sequences= Vec::<&[u8]>::with_capacity(raw_fed_sequences.len() / 2 + 1);
    for i in (1..raw_fed_sequences.len()).step_by(2) {
        headers.push(raw_fed_sequences[i-1]);
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
        temp = positions_to_keep.iter()
            .map(|index|sequences[i][*index])
            .collect();
        result.push(headers[i].to_string());
        result.push(String::from_utf8(temp).unwrap());
    }
    (result, positions_to_keep)
}

// fn delete_empty_columns(raw_input_seqs: Vec<&str>, verbose: bool) -> (Vec<String>, Vec<u32>) {
//     let mut result = Vec::<String>::new();
//     let sequences = raw_input_seqs.iter()
//         .skip(1)
//         .step_by(2)
//         .map(|seq| seq.as_bytes())
//         .collect();
//     let length: usize = sequences[0].len();
//     let mut positions_to_keep = HashSet::<usize>::new();
//     for i in 0..length{
//         for sequence in sequences {
//             if sequence[i] != b'-' {
//                 positions_to_keep.insert(i);
//                 break;
//             }
//         }
//     }
//     let sequences = sequences.iter()
//         .
//     for (header, sequence) in raw_input_seqs.iter().step_by(2).zip(sequences.iter()) {
//         result.push(header.to_string());
//         result.push()
//     }
//     for i in 0..length {
//         for sequence in raw_input_seqs:
//
//     }
//     unimplemented!()
// }
// A Python module implemented in Rust.
#[pymodule]
fn phymmr_tools(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(blosum62_distance, m)?)?;
    m.add_function(wrap_pyfunction!(bio_revcomp, m)?)?;
    m.add_function(wrap_pyfunction!(constrained_distance, m)?)?;
    m.add_function(wrap_pyfunction!(dumb_consensus, m)?)?;
    m.add_function(wrap_pyfunction!(dumb_consensus_with_excise, m)?)?;
    m.add_function(wrap_pyfunction!(dumb_consensus_dupe, m)?)?;
    m.add_function(wrap_pyfunction!(excise_consensus_tail, m)?)?;
    m.add_function(wrap_pyfunction!(find_index_pair, m)?)?;
    m.add_function(wrap_pyfunction!(has_data, m)?)?;
    m.add_function(wrap_pyfunction!(blosum62_candidate_to_reference,m)?)?;
    m.add_function(wrap_pyfunction!(score_splits,m)?)?;
    m.add_function(wrap_pyfunction!(simd_hamming,m)?)?;
    m.add_function(wrap_pyfunction!(delete_empty_columns,m)?)?;
    m.add_function(wrap_pyfunction!(join_by_tripled_index,m)?)?;
    m.add_function(wrap_pyfunction!(join_with_exclusions,m)?)?;
    m.add_function(wrap_pyfunction!(join_triplets_with_exclusions,m)?)?;


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