use pyo3::pyfunction;
use crate::find_indices;

#[pyfunction]
pub fn dumb_consensus(sequences: Vec<String>, threshold: f64, min_depth: u32) -> String {
    consensus(sequences[0].len(), sequences.iter().map(|s| (s.as_str(), 1)), threshold, min_depth)
}

#[pyfunction]
pub fn dumb_consensus_dupe(sequences: Vec<(String, u32)>, threshold: f64, min_depth: u32) -> String {
    consensus(sequences[0].0.len(), sequences.iter().map(|(s, n)| (s.as_str(), *n)), threshold, min_depth)
}

/// The residue of each column that holds more than `threshold` of the column's
/// weight ('-' for gaps and stops), else X; X too where the weight is zero or
/// under `min_depth`. A sequence counts `weight` times from its first to its
/// last non-gap column.
fn consensus<'a>(
    len: usize,
    sequences: impl Iterator<Item = (&'a str, u32)>,
    threshold: f64,
    min_depth: u32,
) -> String {
    const ASCII_OFFSET: u8 = 65;
    const HYPHEN: u8 = 45;
    const ASTERISK: u8 = 42;
    let mut total_at_position = vec![0_u32; len];
    let mut counts_at_position = vec![[0_u32; 27]; len];
    for (sequence, weight) in sequences {
        let seq = sequence.as_bytes();
        let (start, end) = find_indices(seq, b'-');
        for index in start..end {
            if index == seq.len() {
                continue;
            }
            total_at_position[index] += weight;
            if !(seq[index] == HYPHEN || seq[index] == ASTERISK) {
                counts_at_position[index][(seq[index] - ASCII_OFFSET) as usize] += weight;
            } else {
                counts_at_position[index][26] += weight;
            }
        }
    }
    let min_depth = min_depth.max(1);
    let mut output = Vec::<u8>::with_capacity(len);
    for (total, counts) in total_at_position.into_iter().zip(counts_at_position.iter()) {
        if total < min_depth {
            output.push(b'X');
            continue;
        }
        let mut max_count: u32 = 0;
        let mut winner = b'X'; // no residue over the threshold
        for (index, count) in counts.iter().enumerate() {
            if *count as f64 / total as f64 > threshold {
                if *count > max_count {
                    max_count = *count;
                    if index != 26 {
                        winner = index as u8 + ASCII_OFFSET;
                    } else {
                        winner = HYPHEN;
                    }
                }
            }
        }
        output.push(winner);
    }
    String::from_utf8(output).unwrap()
}

fn _mask_small_regions(sequence: &str, min_length: usize) -> String {
    let sequence = sequence.as_bytes();
    let mut start = Option::None;
    let mut regions = Vec::<(usize, usize)>::new();
    let mut output = vec![b'X';sequence.len()];
    for (i, bp) in sequence.iter().enumerate() {
        if *bp == b'-' {
            match start {
                Some(index) => {
                    regions.push((index, i));
                    start = None;
                }
                None => continue
            }
        } else {
            match start {
                Some(_) => {},
                None => start = Some(i),
            }
        }
    }
    match start {
        Some(index) => {
            regions.push((index, sequence.len()));
        }
        None => {}
    }
    for (begin, end) in regions.iter() {
        if end - begin < min_length {
            continue;
        }
        for index in *begin..*end {
            output[index] = sequence[index]
        }
    }
    String::from_utf8(output).unwrap()
}

fn detect_regions(sequence: &str) -> Vec<(usize, usize)> {
    let sequence = sequence.as_bytes();
    let mut start = Option::None;
    let mut regions = Vec::<(usize, usize)>::new();
    // let mut output = vec![b'-';sequence.len()];
    for (i, bp) in sequence.iter().enumerate() {
        if *bp == b'-' {
            match start {
                Some(index) => {
                    regions.push((index, i));
                    start = None;
                }
                None => continue
            }
        } else {
            match start {
                Some(_) => {},
                None => start = Some(i),
            }
        }
    }
    match start {
        Some(index) => {
            regions.push((index, sequence.len()));
        }
        None => {}
    }
    regions
}

fn overlap_and_distance(con: &[u8], can: &[u8]) -> (usize, usize, usize) {
    let mut overlap = 0;
    let mut distance = 0;
    let mut total_length = 0;
    for (con_char,can_char) in con.iter().zip(can.iter()) {
        if !(*con_char == b'X' || *can_char == b'X') { overlap += 1;}
        if *con_char != b'X' && *con_char != *can_char {distance += 1;}
        total_length += 1;
    }
    (overlap, distance, total_length)
}
#[pyfunction]
pub fn consensus_distance(consensus: String, candidate: String, min_length: usize, min_overlap: usize) -> (usize,usize) {
    let can = &_mask_small_regions(&candidate, min_length);
    let regions = detect_regions(&candidate);
    let can = can.as_bytes();
    let con = consensus.as_bytes();
    let mut total_length = 0;
    let mut total_distance = 0;
    for (begin, end) in regions.iter() {
        let (overlap, distance, length) = overlap_and_distance(&con[*begin..*end], &can[*begin..*end]);
        if overlap >= min_overlap {
            total_distance += distance;
            total_length += length;
        }
    }
    (total_distance, total_length)
}
