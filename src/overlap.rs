use pyo3::prelude::*;
use std::cmp::{max, min};

#[pyfunction]
pub fn get_overlap(
    start1: i64,
    end1: i64,
    start2: i64,
    end2: i64,
    min_overlap: i64,
) -> PyResult<Option<(i64, i64)>> {
    let begin = max(start1, start2);
    let end = min(end1, end2);
    match (end - begin) >= min_overlap {
        true => Ok(Some((begin, end))),
        false => Ok(None),
    }
}
