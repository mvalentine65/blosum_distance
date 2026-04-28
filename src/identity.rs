use pyo3::pyfunction;
use std::collections::HashSet;

#[pyfunction]
pub fn filter_nt(
    candidates: Vec<(String, String)>,
    failed: HashSet<String>,
) -> Vec<(String, String)> {
    let mut output = Vec::new();
    for (header, seq) in candidates {
        if !failed.contains(&header) {
            output.push((header, seq))
        }
    }
    output
}
