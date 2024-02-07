use itertools::enumerate;
use pyo3::pyfunction;
use crate::{find_indices, simd_hamming};

#[pyfunction]
pub fn dumb_consensus(sequences: Vec<&str>, threshold: f64, min_depth: u32) -> String {
    match min_depth {
    0 => _dumb_consensus1(sequences,threshold),
    _ => _dumb_consensus2(sequences, threshold, min_depth),
    }
}
#[pyfunction]
pub fn dumb_consensus_dupe(sequences: Vec<(&str, u32)>, threshold: f64, min_depth: u32) -> String {
    match min_depth {
        0 => _dumb_consensus_dupe1(sequences, threshold),
        _ => _dumb_consensus_dupe2(sequences, threshold, min_depth),
    }
}
fn _dumb_consensus1(sequences: Vec<&str>, threshold: f64) -> String {
    let first = &sequences[0];
    let mut total_at_position = vec![0_u32; first.len()];
    let mut counts_at_position = vec![[0_u32; 27]; first.len()];
    const ASCII_OFFSET: u8 = 65;
    const HYPHEN: u8 = 45;
    const ASTERISK: u8 = 42;
    let mut min = usize::MAX;
    let mut max: usize = 0;
    for sequence in sequences.iter() {
        let seq = sequence.as_bytes();
        let (start, end) = find_indices(seq, b'-');
        if start < min {
            min = start;
        }
        if end > max {
            max = end;
        }
        // let seq = &seq[..];
        for index in start..end {
            if index == seq.len() {
                continue;
            }
            total_at_position[index] += 1;
            if !(seq[index] == HYPHEN || seq[index] == ASTERISK) {
                // if seq[index]-ASCII_OFFSET == 233 {println!("{}",seq[index]);}
                counts_at_position[index][(seq[index] - ASCII_OFFSET) as usize] += 1;
                // total_at_position[index] += 1;
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
        let mut max_count: u32 = 0;
        let mut winner = b'X'; // default to X if no winner found
        for (index, count) in enumerate(counts) {
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

fn _dumb_consensus2(sequences: Vec<&str>, threshold: f64, min_depth: u32) -> String {
    let first = &sequences[0];
    let mut total_at_position = vec![0_u32; first.len()];
    let mut counts_at_position = vec![[0_u32; 27]; first.len()];
    const ASCII_OFFSET: u8 = 65;
    const HYPHEN: u8 = 45;
    const ASTERISK: u8 = 42;
    let mut min = usize::MAX;
    let mut max: usize = 0;
    for sequence in sequences.iter() {
        let seq = sequence.as_bytes();
        let (start, end) = find_indices(seq, b'-');
        if start < min {
            min = start;
        }
        if end > max {
            max = end;
        }
        // let seq = &seq[..];
        for index in start..end {
            if index == seq.len() {
                continue;
            }
            // total_at_position[index] += 1;
            if !(seq[index] == HYPHEN || seq[index] == ASTERISK) {
                // if seq[index]-ASCII_OFFSET == 233 {println!("{}",seq[index]);}
                counts_at_position[index][(seq[index] - ASCII_OFFSET) as usize] += 1;
                total_at_position[index] += 1;
            } else {
                counts_at_position[index][26] += 1;
            }
        }
    }
    let mut output = Vec::<u8>::with_capacity(total_at_position.len());
    for ((_, total), counts) in enumerate(total_at_position).zip(counts_at_position.iter()) {
        if total < min_depth {
            output.push(b'X');
            continue;
        } // if no characters at position, continue
        let mut max_count: u32 = 0;
        let mut winner = b'X'; // default to X if no winner found
        for (index, count) in enumerate(counts) {
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

fn _dumb_consensus_dupe1(sequences: Vec<(&str, u32)>, threshold: f64) -> String {
    let (first, _) = &sequences[0];
    let mut total_at_position = vec![0_u32; first.len()];
    let mut counts_at_position = vec![[0_u32; 27]; first.len()];
    const ASCII_OFFSET: u8 = 65;
    const HYPHEN: u8 = 45;
    const ASTERISK: u8 = 42;
    let mut min = usize::MAX;
    let mut max: usize = 0;
    for (sequence, count) in sequences.iter() {
        let seq = sequence.as_bytes();
        let (start, end) = find_indices(seq, b'-');
        if start < min {
            min = start;
        }
        if end > max {
            max = end;
        }
        // let seq = &seq[..];
        for index in start..end {
            if index == seq.len() {
                continue;
            }
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
        let mut max_count: u32 = 0;
        let mut winner = b'X'; // default to X if no winner found
        for (index, count) in enumerate(counts) {
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
fn _dumb_consensus_dupe2(sequences: Vec<(&str, u32)>, threshold: f64, min_depth: u32) -> String {
    let (first, _) = &sequences[0];
    let mut total_at_position = vec![0_u32; first.len()];
    let mut counts_at_position = vec![[0_u32; 27]; first.len()];
    const ASCII_OFFSET: u8 = 65;
    const HYPHEN: u8 = 45;
    const ASTERISK: u8 = 42;
    let mut min = usize::MAX;
    let mut max: usize = 0;
    for (sequence, count) in sequences.iter() {
        let seq = sequence.as_bytes();
        let (start, end) = find_indices(seq, b'-');
        if start < min {
            min = start;
        }
        if end > max {
            max = end;
        }
        // let seq = &seq[..];
        for index in start..end {
            if index == seq.len() {
                continue;
            }
            if !(seq[index] == HYPHEN || seq[index] == ASTERISK) {
                // if seq[index]-ASCII_OFFSET == 233 {println!("{}",seq[index]);}
                counts_at_position[index][(seq[index] - ASCII_OFFSET) as usize] += count;
                total_at_position[index] += count;
            } else {
                counts_at_position[index][26] += count;
            }
        }
    }
    let mut output = Vec::<u8>::with_capacity(total_at_position.len());
    for ((_, total), counts) in enumerate(total_at_position).zip(counts_at_position.iter()) {
        if total < min_depth {
            output.push(b'X');
            continue;
        } // if no characters at position, continue
        let mut max_count: u32 = 0;
        let mut winner = b'X'; // default to X if no winner found
        for (index, count) in enumerate(counts) {
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


#[pyfunction]
pub fn filter_regions(sequence: String, min_length: usize) -> String {
    _mask_small_regions(&sequence, min_length)
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


fn overlap_and_distance(con: &[u8], can: &[u8]) -> (usize, usize) {
    let mut overlap = 0;
    let mut distance = 0;
    for (con_char,can_char) in con.iter().zip(can.iter()) {
        if *con_char == b'X' || *can_char == b'-' { continue;}
        overlap += 1;
        if *con_char != *can_char {distance += 1;}
    }
    (overlap, distance)
}
#[pyfunction]
pub fn consensus_distance(consensus: String, candidate: String, min_length: usize) -> u64 {
    let can = &filter_regions(candidate, min_length);
    simd_hamming(&can, &consensus)
}
// #[pyfunction]
// fn dumb_consensus_with_excise(
//     sequences: Vec<&str>,
//     consensus_threshold: f64,
//     min_depth: u32,
//     excise_threshold: f64,
// ) -> (String, usize, String) {
//     let first = &sequences[0];
//     let mut total_at_position = vec![0_u32; first.len()];
//     let mut counts_at_position = vec![[0_u32; 27]; first.len()];
//     const ASCII_OFFSET: u8 = 65;
//     const HYPHEN: u8 = 45;
//     const ASTERISK: u8 = 42;
//     let mut min = usize::MAX;
//     let mut max: usize = 0;
//     for sequence in sequences.iter() {
//         let seq = sequence.as_bytes();
//         let (start, end) = find_indices(seq, b'-');
//         if start < min {
//             min = start;
//         }
//         if end > max {
//             max = end;
//         }
//         // let seq = &seq[..];
//         for index in start..end {
//             if index == seq.len() {
//                 continue;
//             }
//             if !(seq[index] == HYPHEN || seq[index] == ASTERISK) {
//                 // if seq[index]-ASCII_OFFSET == 233 {println!("{}",seq[index]);}
//                 counts_at_position[index][(seq[index] - ASCII_OFFSET) as usize] += 1;
//                 total_at_position[index] += 1;
//             } else {
//                 counts_at_position[index][26] += 1;
//             }
//         }
//     }
//     let mut output = Vec::<u8>::with_capacity(total_at_position.len());
//     // let mut ratios = Vec::<u8>::with_capacity(total_at_position.len());
//     // let mut ratio = Vec::<u8>::with_capacity(to)
//     for ((_, total), counts) in enumerate(total_at_position).zip(counts_at_position.iter()) {
//         if total < min_depth {
//             output.push(b'X');
//             continue;
//         } // if no characters at position, continue
//         let mut max_count: u32 = 0;
//         let mut winner = b'X'; // default to X if no winner found
//         for (index, count) in enumerate(counts) {
//             if *count as f64 / total as f64 > consensus_threshold {
//                 if *count > max_count {
//                     max_count = *count;
//                     if index != 26 {
//                         winner = index as u8 + ASCII_OFFSET;
//                     } else {
//                         winner = HYPHEN;
//                     }
//                 }
//             }
//         }
//
//         output.push(winner);
//     }
//     // let locations: Vec<LocationData> = counts_at_position.iter()
//     //     .map(|letters| weigh_winner(letters))
//     //     .collect();
//     let consensus = String::from_utf8(output).unwrap();
//     let (excised, cut_length) = _excise_consensus_tail(&consensus, excise_threshold);
//     (excised, cut_length, consensus)
// }