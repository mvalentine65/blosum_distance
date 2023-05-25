use std::collections::HashSet;
use pyo3::prelude::*;

#[pyfunction]
pub fn join_by_tripled_index(string: &str, positions_to_keep: Vec<usize>) -> String {
    let mut result = Vec::with_capacity(string.len());
    for i in positions_to_keep.iter() {
        result.push(&string[i*3..i*3+3]);
    }
    result.join("")
}

#[pyfunction]
pub fn join_with_exclusions(string: &str, column_cull: HashSet<usize>) -> String {
    let chars = string.as_bytes();
    let mut result = Vec::with_capacity(string.len());
    for i in 0..string.len() {
        match column_cull.contains(&(i * 3)) {
            true => result.push(b'-'),
            false => result.push(chars[i]),
        }
    }
    String::from_utf8(result).unwrap()
}

#[pyfunction]
pub fn join_triplets_with_exclusions(string: &str, exclusion1: HashSet<usize>, exclusion2: HashSet<usize>) -> String {
    let mut result = Vec::with_capacity(string.len());
    for i in (0..string.len()).step_by(3) {
        match exclusion1.contains(&i) || exclusion2.contains(&i){
            true => result.push("---"),
            false => result.push(&string[i..i+3]),
        }
    }
    result.join("")
}