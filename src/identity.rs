use std::collections::HashSet;
use bio::io::fasta;
use pyo3::pyfunction;

#[pyfunction]
pub fn read_file(input_path: &str) -> (Vec<(Vec<u8>, Vec<u8>)>, Vec<(Vec<u8>, Vec<u8>)>) {
    let mut candidates = Vec::new();
    let mut references = Vec::new();
    let reader = fasta::Reader::from_file(input_path).unwrap();

    for record in reader.records() {
        if record.is_err() {continue;}
        let record = record.unwrap();
        let header = record.id().as_bytes().to_owned();
        let sequence = record.seq().to_ascii_uppercase();
        if header[header.len()-1] != b'.' {
            candidates.push((header, sequence));
        } else {
            references.push((header, sequence))
        }
    }
    (candidates, references)
}

#[pyfunction]
pub fn filter_nt(candidates: Vec<(String, String)>, failed:HashSet<String>) -> Vec<(String, String)>{
    let mut output = Vec::new();
    for (header, seq) in candidates{
        if !failed.contains(&header) {
            output.push((header, seq))
        }
    }
    output
}
