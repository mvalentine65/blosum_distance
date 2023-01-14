#![allow(unused)]

extern crate bio;
extern crate hashbrown;
extern crate pyo3;

use core::num;
use std::borrow::BorrowMut;
// use statrs::statistics::OrderStatistics;
// use statrs::statistics::Data;
use std::fs::File;
use std::io::{prelude::*, BufReader};
// use std::sync::Arc;
//use bio::alphabets::dna::complement;
use bio::io::fasta;
use bio::scores::blosum62;
use itertools::Itertools;
use proc_utils::{add_x_members, repeat_across_x_times, repeat_across_x_times_and_create_struct};
use pyo3::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::IoSliceMut;
use std::ptr::addr_of_mut;

// use rayon::prelude::*;

// Lifetime annotations for object lifetime c
// and some lifetime s which is greater than c.
// If you prefer, c is the current scope,
// and s is the outer scope.
// #[pyclass]
// struct HammingCluster<'c>{ // do i really need this if its one vec?
//     leader: HammingRecord<'c>,
//     others: Vec<HammingRecord<'c>>
// #
// #[new]
//
// }

// #[derive(Hash, Debug, Clone)]
// struct HammingRecord<'c> {
//     header: &'c str,
//     sequence: &'c str,
//     dupe_count: u8,
// }

macro_rules! str_make {
    ($s:literal) => {
        String::from($s)
    };
}

macro_rules! hm_make {
    (
        $key_type:ty,
        $value_type:ty,
        $($key:expr => $value:expr), *) => {
            {
                let mut hm = HashMap::<$key_type, $value_type>::new();

                $
                (
                    hm.insert($key, $value);
                )
                *

                hm
            }

    };
    (
        $key_type:ty,
        $value_type:ty,
        $($key:literal => $value:literal), *) => {
            {
                let mut hm = HashMap::<$key_type, $value_type>::new();

                $
                (
                    hm.insert($key, $value);
                )
                *

                hm
            }

    };
}

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
fn make_indices(sequence: String) -> (usize, usize) {
    let start: usize = sequence
        .chars()
        .enumerate()
        .filter(|(_, bp)| *bp != '-')
        .map(|(index, _)| index)
        .take(1)
        .next()
        .unwrap();
    let delta_end: usize = sequence
        .chars()
        .rev()
        .enumerate()
        .filter(|(_, bp)| *bp != '-')
        .map(|(index, _)| index)
        .take(1)
        .next()
        .unwrap();
    (start, sequence.len() - delta_end)
}

#[pyfunction]
fn fasta_reader(path: String) -> Vec<String> {
    let mut result: Vec<String> = Vec::new();
    let fasta_reader = fasta::Reader::from_file(path).unwrap();
    for fasta in fasta_reader.records() {
        let record = &fasta.unwrap();
        result.push(format!(">{}", record.id()));
        result.push(String::from_utf8(record.seq().to_vec()).unwrap());
    }
    result
}

#[pyfunction]
fn cluster_distance_filter(lines: Vec<&str>) -> Vec<Vec<&str>> {
    let mut clusters = HashMap::new();
    //let lines = lines_vector.as_slice();
    let mut line_reader = lines.iter().cloned();
    // while we have more lines to pull records from...
    while let (Some(name), Some(seq)) = (line_reader.next(), line_reader.next()) {
        let key = (seq.len(), &seq[0..10]);
        let array = vec![name, seq];
        // let ham: HammingRecord = HammingRecord{
        //     header: name,
        //     sequence: seq,
        //     dupe_count: 0, // change later if sort is implemented
        // };
        // let mut ham_cluster = HammingCluster {
        //     leader:ham,
        //     others:Vec::new(),
        // };
        clusters.entry(key).or_insert(Vec::new()).push(array);
    }
    let mut output = Vec::new();
    for (_key, mut cluster) in clusters {
        // Pops Hamming records off the cluster vector and makes the new record
        // the leader in a HammingCluster. Compares leader to other values by distance
        // and
        loop {
            // if we pop a None, the vector is empty, so break
            let lead = cluster.pop();
            if lead.is_none() {
                break;
            }
            let mut lead = lead.unwrap();
            // if lead.is_none() { continue; }
            // let group = HammingCluster {
            //     leader: lead.unwrap(),
            //     others: Vec::new(),
            // };
            // Iterate over other seqs in reverse.
            // If the seq is within distance, swap
            // with the last element and pop it, then
            // append to the lead seq's vector.
            // This works because we start at the end of the vec anyway
            // by the time we hit an inner element within distance,
            // we have already verified the last element is not within distance.
            for index in (0..cluster.len()).rev() {
                let candidate = &cluster[index];
                // if candidate.is_none() {continue;}
                // let other = candidate.unwrap();
                if seqs_within_distance(candidate[1], lead[1], 1) {
                    // group.others.push(other)
                    lead.extend(cluster.swap_remove(index));
                }
            }
            output.push(lead)
        }
    }
    output
}

#[pyfunction]
fn batch_reverse_complement(list: Vec<String>) -> Vec<String> {
    list.into_iter().map(bio_revcomp).collect()
}

#[pyfunction]
fn bio_revcomp(sequence: String) -> String {
    String::from_utf8(bio::alphabets::dna::revcomp(sequence.into_bytes())).unwrap()
}

#[pyfunction]
fn seqs_within_distance(first: &str, second: &str, max_distance: u32) -> bool {
    let (array_one, array_two) = (first.as_bytes(), second.as_bytes());
    if array_one.len() != array_two.len() {
        return false;
    }
    let mut distance: u32 = 0;
    for i in 0..array_one.len() {
        if array_one[i] != array_two[i] {
            distance += 1;
            if distance > max_distance {
                return false;
            }
        }
    }
    true
}

#[add_x_members(buff_*: Vec<u8>, => 256)]
struct Buffers {}

impl Buffers {
    pub fn new() -> Self {
        repeat_across_x_times! {
            let buff_* = vec![0u8; 8192]; @ 256
        };

        repeat_across_x_times_and_create_struct!(
            Self { * } @ buff_*, @ 256
        )
    }
}

struct FastLineReader {
    bufreader: BufReader<File>,
    iovecs: Vec<IoSliceMut<'static>>,
    str_buffers: Vec<String>,
    remainder: String,
    remaining_vec: Vec<usize>,
}

impl FastLineReader {
    pub fn new(path: String, buffers_ptr: *mut Buffers) -> Self {
        println!("Starting new FastLineReader instance on path {path}");

        let buffers: &mut Buffers = unsafe { &mut *buffers_ptr };

        let bufreader = BufReader::new(File::open(path).expect("Error opening output file"));

        let mut iovecs: Vec<IoSliceMut<'static>> = Vec::new();

        repeat_across_x_times! {
            iovecs.push(IoSliceMut::new(&mut buffers.buff_*)); @ 256
        };

        let str_buffers = vec![str_make!(""); 256];

        Self {
            bufreader,
            iovecs,
            str_buffers,
            remaining_vec: vec![],
            remainder: str_make!(""),
        }
    }

    pub fn fill_and_get_lines(&mut self, index: usize) {
        let curr_buffer = self
            .str_buffers
            .get_mut(index)
            .expect("Error getting buffer");
        let curr_vio = self.iovecs.get(index).expect("Error getting vecio");

        curr_vio.iter().cloned().for_each(|c| {
            curr_buffer.push(c as char);
        });
    }

    pub fn read_batch(&mut self) -> Vec<String> {
        let rem = self
            .bufreader
            .read_vectored(&mut self.iovecs)
            .expect("Error with Vector Read");

        self.remaining_vec.push(rem);

        if rem == 0 {
            if self.remainder.len() != 0 {
                return vec![self.remainder.clone()];
            }

            return vec![];
        }

        let ret_list = (0..256)
            .into_iter()
            .for_each(|i| self.fill_and_get_lines(i));

        let strings_together = self.str_buffers.join("");
        let with_rem = format!("{}{}{}", self.remainder, strings_together, "\n");

        let mut split: Vec<String> = vec![];
        let mut buff = str_make!("");

        with_rem.chars().for_each(|c| {
            if c == '\n' {
                let trimmed = buff.trim();
                split.push(trimmed.to_string());
                buff.clear();
            } else {
                if c != '\0' {
                    buff.push(c);
                }
            }
        });

        if *self.remaining_vec.first().unwrap() == 8192 * 256
            && self.remaining_vec.len() > 1
            && *self.remaining_vec.last().unwrap() == 8192 * 256
        {
            self.remainder = split.pop().unwrap();
        }

        split
    }
}

#[pyclass]
struct VecIO {
    freader: Option<FastLineReader>,
    buffers: Buffers,
}

#[pymethods]
impl VecIO {
    #[new]
    pub fn new() -> Self {
        Self {
            buffers: Buffers::new(),
            freader: None,
        }
    }

    pub fn init(&mut self, path: String) {
        let buffers_ptr: *mut Buffers = addr_of_mut!(self.buffers);

        self.freader = Some(FastLineReader::new(path, buffers_ptr));
    }

    pub fn is_done(&mut self) -> bool {
        let derefed = self.freader.as_mut().unwrap();

        *derefed.remaining_vec.last().unwrap() == 0
    }

    pub fn get_next_batch(&mut self) -> Vec<String> {
        let derefed = self.freader.as_mut().unwrap();

        derefed.read_batch()
    }
}

#[pyfunction]
fn get_vecio(path: String) -> VecIO {
    let mut vio = VecIO::new();
    vio.init(path);

    vio
}

fn not_skip_character(character: u8) -> bool {
    const HYPHEN: u8 = 45;
    const ASTERISK: u8 = 42;
    character != HYPHEN && character != ASTERISK
}

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

#[pyfunction]
fn blosum62_distance(one: String, two: String) -> PyResult<f64> {
    let first: &[u8] = one.as_bytes();
    let second: &[u8] = two.as_bytes();
    let mut score = 0;
    let mut max_first = 0;
    let mut max_second = 0;
    let length = first.len();
    let allowed: HashSet<u8> = HashSet::from([
        65, 84, 67, 71, 73, 68, 82, 80, 87, 77, 69, 81, 83, 72, 86, 76, 75, 70, 89, 78, 88, 90, 74,
        66, 79, 85,
    ]);
    for i in 0..length {
        // if !(first[i] == HYPHEN || second[i] == HYPHEN) {
        if not_skip_character(first[i]) && not_skip_character(second[i]) {
            if !(allowed.contains(&first[i])) {
                panic!("first[i]  {} not in allowed\n{}", first[i] as char, one);
            }
            if !(allowed.contains(&second[i])) {
                panic!("second[i] {} not in allowed\n{}", second[i] as char, two);
            }
            // println!("score = {}\nmax_first = {}\nmax_second = {}\n first[i] = {}\n second[i] = {}\n",score,max_first,max_second,first[i],second[i]);
            score += bio::scores::blosum62(first[i], second[i]);
            max_first += bio::scores::blosum62(first[i], first[i]);
            max_second += bio::scores::blosum62(second[i], second[i]);
        }
    }
    // println!("score: {}", score);
    // println!("max_first: {}", max_first);
    // println!("max_second: {}", max_second);
    let maximum_score = std::cmp::max(max_first, max_second);
    Ok(1.0 - (score as f64 / maximum_score as f64))
}

// A Python module implemented in Rust.
#[pymodule]
fn phymmr_tools(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(blosum62_distance, m)?)?;
    m.add_function(wrap_pyfunction!(batch_reverse_complement, m)?)?;
    m.add_function(wrap_pyfunction!(bio_revcomp, m)?)?;
    m.add_function(wrap_pyfunction!(fasta_reader, m)?)?;
    m.add_function(wrap_pyfunction!(cluster_distance_filter, m)?)?;
    m.add_function(wrap_pyfunction!(seqs_within_distance, m)?)?;
    m.add_function(wrap_pyfunction!(find_references_and_candidates, m)?)?;
    m.add_function(wrap_pyfunction!(make_indices, m)?)?;
    m.add_function(wrap_pyfunction!(make_ref_distance_matrix_supervector, m)?)?;
    m.add_function(wrap_pyfunction!(make_ref_distance_vector_blosum62, m)?)?;
    m.add_function(wrap_pyfunction!(ref_index_vector, m)?)?;
    m.add_function(wrap_pyfunction!(get_vecio, m)?)?;

    Ok(())
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//
//     fn test_ref_distance_array
// }
