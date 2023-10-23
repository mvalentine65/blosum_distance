use pyo3::prelude::*;
use std::borrow::Cow;
use std::cmp::{max, min};

struct Seq1d<'a> {
    header: Cow<'a, str>,
    seq: Cow<'a, str>,
    start: usize,
    stop: usize,
}

#[pyfunction]
pub fn get_overlap(
    start1: i32,
    end1: i32,
    start2: i32,
    end2: i32,
    min_overlap: i32,
) -> PyResult<Option<(i32, i32)>> {
    let begin = max(start1, start2);
    let end = min(end1, end2);
    match (end - begin) >= min_overlap {
        true => Ok(Some((begin, end))),
        false => Ok(None),
    }
}

#[pyfunction]
pub fn get_overlap_percent(start1: i32, end1: i32, start2: i32, end2: i32) -> f32 {
    let distance = min(&end1, &end2) - max(&start1, &start2);
    match distance < 0 {
        false => distance as f32 / min(&end1 - &start1, &end2 - &start2) as f32,
        true => 0.0,
    }
}
