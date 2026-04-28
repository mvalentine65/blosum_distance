use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::io::{BufWriter, Write};
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

/// Build a FASTA blob from records and write it to `out` in a single buffered
/// stream.  Mirrors the Python `fp.write("".join(...).encode())` pattern but
/// avoids the intermediate Python string allocation.
fn write_fasta_to_writer<W: Write>(out: &mut W, records: &[(String, String)]) -> std::io::Result<()> {
    let total: usize = records
        .iter()
        .map(|(h, s)| h.len() + s.len() + 3)
        .sum();
    let mut buf = Vec::with_capacity(total);
    for (header, seq) in records {
        buf.push(b'>');
        buf.extend_from_slice(header.as_bytes());
        buf.push(b'\n');
        buf.extend_from_slice(seq.as_bytes());
        buf.push(b'\n');
    }
    out.write_all(&buf)
}

fn run_command(cmd: &mut Command, name: &str) -> Result<(), String> {
    // Default to discarding stderr; only re-run with a captured pipe if the
    // command failed, so the happy path avoids pipe creation entirely.
    let status = cmd
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map_err(|e| format!("{} failed to start: {}", name, e))?;
    if status.success() {
        return Ok(());
    }

    let output = cmd
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .output()
        .map_err(|e| format!("{} failed to start (rerun): {}", name, e))?;
    let stderr = String::from_utf8_lossy(&output.stderr);
    Err(format!("{} failed: {}", name, stderr))
}

/// Run hmmbuild + hmmalign on candidate sequences against a reference alignment.
///
/// Returns aligned (header, sequence) tuples with insertion dots normalised to
/// dashes and all residues uppercased.
///
/// `tmpdir` selects where the four scratch files live.  Pass the same path
/// SAPPHYRE's Python side resolves with `get_temp_dir()` (``/dev/shm`` when
/// available) so we keep the I/O off spinning disks.
///
/// `gene_name`, when provided, is embedded in each scratch file's prefix and
/// used as the HMM's internal name.  This avoids collisions when many workers
/// run hmmbuild/hmmalign concurrently against the same tmpdir.
#[pyfunction]
#[pyo3(signature = (candidates, references, tmpdir = None, gene_name = None))]
pub fn hmm_align(
    py: Python<'_>,
    candidates: Vec<(String, String)>,
    references: Vec<(String, String)>,
    tmpdir: Option<String>,
    gene_name: Option<String>,
) -> PyResult<Vec<(String, String)>> {
    py.detach(move || {
        hmm_align_inner(candidates, references, tmpdir, gene_name)
    })
    .map_err(PyRuntimeError::new_err)
}

fn hmm_align_inner(
    candidates: Vec<(String, String)>,
    references: Vec<(String, String)>,
    tmpdir: Option<String>,
    gene_name: Option<String>,
) -> Result<Vec<(String, String)>, String> {
    // Sanitise the gene name so it's safe inside a filename — strip anything
    // that isn't alphanumeric/_/-/. so we don't accidentally inject path
    // separators or shell metacharacters into the tempfile prefix.
    let slug: String = gene_name
        .as_deref()
        .unwrap_or("")
        .chars()
        .filter(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.'))
        .collect();
    let tag = if slug.is_empty() {
        String::new()
    } else {
        format!("{}_", slug)
    };

    let make_temp = |kind: &str, suffix: &str| {
        let prefix = format!("{}{}_", tag, kind);
        let mut b = Builder::new();
        b.prefix(&prefix).suffix(suffix);
        match &tmpdir {
            Some(d) => b.tempfile_in(d),
            None => b.tempfile(),
        }
        .map_err(|e| e.to_string())
    };

    let mut temp_aln = make_temp("aln", ".fa")?;
    let temp_hmm = make_temp("hmm", ".hmm")?;
    let mut temp_cand = make_temp("cand", ".fa")?;
    let temp_result = make_temp("res", ".afa")?;

    // Write through the existing tempfile handles - no second open().
    {
        let mut w = BufWriter::with_capacity(1 << 20, temp_aln.as_file_mut());
        write_fasta_to_writer(&mut w, &references).map_err(|e| e.to_string())?;
        w.flush().map_err(|e| e.to_string())?;
    }
    {
        let mut w = BufWriter::with_capacity(1 << 20, temp_cand.as_file_mut());
        write_fasta_to_writer(&mut w, &candidates).map_err(|e| e.to_string())?;
        w.flush().map_err(|e| e.to_string())?;
    }

    let aln_path = temp_aln.path().to_str().unwrap().to_string();
    let hmm_path = temp_hmm.path().to_str().unwrap().to_string();
    let cand_path = temp_cand.path().to_str().unwrap().to_string();
    let result_path = temp_result.path().to_str().unwrap().to_string();

    // hmmbuild
    //
    // --cpu 1 keeps each invocation single-threaded; SAPPHYRE drives this
    // function from a multiprocessing pool, so HMMER's default of 2 pthreads
    // per call would oversubscribe (N_workers x 2 threads).
    let hmm_name = if slug.is_empty() { "hmm" } else { slug.as_str() };
    let mut hmmbuild = Command::new("hmmbuild");
    hmmbuild
        .args(["-n", hmm_name])
        .args(["--cpu", "1"])
        .arg(&hmm_path)
        .arg(&aln_path);
    run_command(&mut hmmbuild, "hmmbuild")?;

    // hmmalign --mapali (hmmalign is single-threaded in HMMER3 — no --cpu).
    let mut hmmalign = Command::new("hmmalign");
    hmmalign.args([
        "--mapali",
        &aln_path,
        "--outformat",
        "afa",
        "-o",
        &result_path,
        &hmm_path,
        &cand_path,
    ]);
    run_command(&mut hmmalign, "hmmalign")?;

    // Parse and normalise output in place: '.' -> '-', uppercase residues.
    let recs = parse_fasta_file(&result_path);
    Ok(recs
        .into_iter()
        .map(|(header, seq)| {
            let mut bytes = seq.into_bytes();
            for b in bytes.iter_mut() {
                if *b == b'.' {
                    *b = b'-';
                } else {
                    b.make_ascii_uppercase();
                }
            }
            (header, unsafe { String::from_utf8_unchecked(bytes) })
        })
        .collect())
}
