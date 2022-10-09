extern crate bio;
extern crate pyo3;

use bio::alphabets::dna::complement;
use pyo3::prelude::*;
use std::collections::HashSet;

#[pyfunction]
fn batch_reverse_complement(list: Vec<String>) -> Vec<String> {
    list.into_iter()
        .map(|sequence| reverse_complement(sequence))
        .collect()
}

#[pyfunction]
fn reverse_complement(sequence: String) -> String {
    let array: Vec<u8> = sequence.into_bytes();
    //let max: usize = array.len() -1;
    let mut complement: Vec<u8> = Vec::new();
    for ascii in array.iter().rev() {
        complement.push( match ascii {
            65 => 84, // A => T
            84 => 65, // T => A
            67 => 71, // C => G
            71 => 67, // G => C
            78 => 78, // C => G
            _ => panic!("Invalid character found")
            // panic! is fine for now, but an error is the proper response for production
        }); // switch case should be at least as fast as a HashMap
    }
    String::from_utf8(complement).unwrap()
}
#[pyfunction]
fn bio_revcomp(sequence: String) -> String {
    String::from_utf8(bio::alphabets::dna::revcomp(sequence.into_bytes())).unwrap()
}

fn not_skip_character(character: u8) -> bool {
    const HYPHEN: u8 = 45;
    const ASTERISK: u8 = 42;
    character != HYPHEN && character != ASTERISK
}

#[pyfunction]
fn blosum62_distance(one: String, two: String) -> PyResult<f64>{
    let first: &[u8] = one.as_bytes();
    let second: &[u8] = two.as_bytes();
    let mut score: i128 = 0;
    let mut max_first: i128 = 0;
    let mut max_second: i128 = 0;
    let length = first.iter().count();
    let allowed: HashSet<u8> = HashSet::from([65,84,67,71,73,68,82,
        80,87,77,69,81,83,72,86,76,75,70,89,78,88]);
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
            score += bio::scores::blosum62(first[i], second[i]) as i128;
            max_first += bio::scores::blosum62(first[i], first[i]) as i128;
            max_second += bio::scores::blosum62(second[i], second[i]) as i128;
        }
    };
    let maximum_score = std::cmp::max(max_first, max_second);
    Ok(1.0 - (score as f64 / maximum_score as f64))
}


// A Python module implemented in Rust.
#[pymodule]
fn blosum_distance(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(blosum62_distance, m)?)?;
    m.add_function(wrap_pyfunction!(reverse_complement, m)?)?;
    m.add_function(wrap_pyfunction!(batch_reverse_complement, m)?)?;
    m.add_function(wrap_pyfunction!(bio_revcomp, m)?)?;
    Ok(())
}
