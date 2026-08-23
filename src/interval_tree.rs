use std::collections::HashSet;
use pyo3::prelude::*;

#[pyfunction]
pub fn del_cols(sequence: String, x_positions: HashSet<usize>, is_nt: bool) -> String {
    match is_nt {
        true => nt_delete_columns(sequence, x_positions),
        false => aa_delete_columns(sequence, x_positions)
    }
}

fn nt_delete_columns(sequence: String, x_positions: HashSet<usize>) -> String {
    let mut bytes = sequence.as_bytes().to_vec();
    for index in x_positions {
        bytes[index] = b'-';
        bytes[index+1] = b'-';
        bytes[index+2] = b'-';
    }
    String::from_utf8(bytes).unwrap()
}

fn aa_delete_columns(sequence: String, x_positions: HashSet<usize>) -> String {
    let mut bytes = sequence.as_bytes().to_vec();
    for index in x_positions {
        bytes[index] = b'-';
    }
    String::from_utf8(bytes).unwrap()
}