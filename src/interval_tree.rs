use std::collections::HashSet;
use meminterval::{Interval, IntervalTree};
use pyo3::prelude::*;

#[pyclass]
pub struct OverlapTree {
    pub tree: IntervalTree<i32, ()>,
}

#[pymethods]
impl OverlapTree {
    #[new]
    fn new() -> Self {
        OverlapTree {
            tree: IntervalTree::new()
        }
    }
    
    fn insert(&mut self, tuple: (i32, i32)) {
        self.tree.insert(Interval::new(tuple.0, tuple.1), ());
    }
    
    fn insert_vector(&mut self, iterable: Vec<(i32, i32)>) {
        iterable.iter()
            .for_each(|tup| self.insert(*tup));
    }
    
    fn query_overlap(&mut self, tuple: (i32, i32)) -> Vec<(i32, i32)>{
        self.tree
            .query(Interval::new(tuple.0, tuple.1))
            .map(|entry| (entry.interval.start, entry.interval.end))
            .collect()
    }
}

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