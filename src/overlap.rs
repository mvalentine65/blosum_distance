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
