use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};
use tempfile::Builder;

use fastx::FastX;

fn parse_fasta_file(path: &str) -> Vec<(String, String)> {
    let mut records = Vec::new();
    let mut reader = FastX::reader_from_path(Path::new(path)).unwrap();
    let mut rec = FastX::from_reader(&mut reader).unwrap();
    while let Ok(1..=usize::MAX) = rec.read(&mut reader) {
        records.push((
            rec.id().to_string(),
            String::from_utf8(rec.seq()).unwrap(),
        ));
    }
    records
}

fn write_fasta_to(path: &str, records: &[(String, String)]) {
    let mut file = std::fs::File::create(path).unwrap();
    for (header, seq) in records {
        write!(file, ">{}\n{}\n", header, seq).unwrap();
    }
    file.flush().unwrap();
}

/// Run hmmbuild + hmmalign on candidate sequences against a reference alignment.
///
/// Returns aligned (header, sequence) tuples with insertion dots normalised to
/// dashes and all residues uppercased.
#[pyfunction]
pub fn hmm_align(
    candidates: Vec<(String, String)>,
    references: Vec<(String, String)>,
) -> PyResult<Vec<(String, String)>> {
    let temp_aln = Builder::new()
        .prefix("aln_")
        .suffix(".fa")
        .tempfile()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let temp_hmm = Builder::new()
        .prefix("hmm_")
        .suffix(".hmm")
        .tempfile()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let temp_cand = Builder::new()
        .prefix("cand_")
        .suffix(".fa")
        .tempfile()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let temp_result = Builder::new()
        .prefix("res_")
        .suffix(".afa")
        .tempfile()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    let aln_path = temp_aln.path().to_str().unwrap().to_string();
    let hmm_path = temp_hmm.path().to_str().unwrap().to_string();
    let cand_path = temp_cand.path().to_str().unwrap().to_string();
    let result_path = temp_result.path().to_str().unwrap().to_string();

    // Write reference alignment
    write_fasta_to(&aln_path, &references);

    // hmmbuild
    let output = Command::new("hmmbuild")
        .args(["-n", "hmm"])
        .arg(&hmm_path)
        .arg(&aln_path)
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .output()
        .map_err(|e| PyRuntimeError::new_err(format!("hmmbuild failed to start: {}", e)))?;
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(PyRuntimeError::new_err(format!(
            "hmmbuild failed: {}",
            stderr
        )));
    }

    // Write candidate sequences
    write_fasta_to(&cand_path, &candidates);

    // hmmalign --mapali
    let output = Command::new("hmmalign")
        .args([
            "--mapali",
            &aln_path,
            "--outformat",
            "afa",
            "-o",
            &result_path,
            &hmm_path,
            &cand_path,
        ])
        .stderr(Stdio::piped())
        .output()
        .map_err(|e| PyRuntimeError::new_err(format!("hmmalign failed to start: {}", e)))?;
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(PyRuntimeError::new_err(format!(
            "hmmalign failed: {}",
            stderr
        )));
    }

    // Parse and normalise output
    let recs = parse_fasta_file(&result_path);
    Ok(recs
        .into_iter()
        .map(|(header, seq)| {
            let normalised = seq.replace('.', "-").to_ascii_uppercase();
            (header, normalised)
        })
        .collect())
}
