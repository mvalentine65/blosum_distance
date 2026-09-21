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

/// Build a FASTA blob from records and write it to `out` in a single buffered
/// stream.  Mirrors the Python `fp.write("".join(...).encode())` pattern but
/// avoids the intermediate Python string allocation.
///
/// When `mask_stops` is true, stop codons (`*`) are rewritten to `X` on the
/// fly so bathbuild doesn't choke on them.  We only do this for the
/// reference alignment; candidate `*`s are passed through verbatim.
fn write_fasta_to_writer<W: Write>(
    out: &mut W,
    records: &[(String, String)],
    mask_stops: bool,
) -> std::io::Result<()> {
    let total: usize = records
        .iter()
        .map(|(h, s)| h.len() + s.len() + 3)
        .sum();
    let mut buf = Vec::with_capacity(total);
    for (header, seq) in records {
        buf.push(b'>');
        buf.extend_from_slice(header.as_bytes());
        buf.push(b'\n');
        let start = buf.len();
        buf.extend_from_slice(seq.as_bytes());
        if mask_stops {
            for b in &mut buf[start..] {
                if *b == b'*' {
                    *b = b'X';
                }
            }
        }
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

/// Run bathbuild + bathalign on candidate sequences against a reference alignment.
///
/// Returns aligned (header, sequence) tuples with insertion dots normalised to
/// dashes and all residues uppercased.
///
/// `tmpdir` selects where the four scratch files live.  Pass the same path
/// SAPPHYRE's Python side resolves with `get_temp_dir()` (``/dev/shm`` when
/// available) so we keep the I/O off spinning disks.
///
/// `gene_name` and `taxa`, when provided, are embedded in each scratch file's
/// prefix; `gene_name` is also used as the HMM's internal name.  Tagging by
/// taxa+gene avoids collisions when many workers run bathbuild/bathalign
/// concurrently against the same tmpdir for different runs of the same gene.
///
/// `cached_hmm` + `cached_template` opt into the bhmm cache hmmsearch fills:
/// the model and the exact MSA it was built from.  Both must be given and both
/// must exist, or we fall back to building.  When they are used `references` is
/// ignored -- `--mapali` re-derives the template's checksum and bathalign dies
/// on a mismatch, so model and template can only travel as a pair.
#[pyfunction]
#[pyo3(signature = (candidates, references, tmpdir = None, gene_name = None, taxa = None, cached_hmm = None, cached_template = None))]
pub fn hmm_align(
    py: Python<'_>,
    candidates: Vec<(String, String)>,
    references: Vec<(String, String)>,
    tmpdir: Option<String>,
    gene_name: Option<String>,
    taxa: Option<String>,
    cached_hmm: Option<String>,
    cached_template: Option<String>,
) -> PyResult<Vec<(String, String)>> {
    py.detach(move || {
        hmm_align_inner(
            candidates,
            references,
            tmpdir,
            gene_name,
            taxa,
            cached_hmm,
            cached_template,
        )
    })
    .map_err(PyRuntimeError::new_err)
}

fn hmm_align_inner(
    candidates: Vec<(String, String)>,
    references: Vec<(String, String)>,
    tmpdir: Option<String>,
    gene_name: Option<String>,
    taxa: Option<String>,
    cached_hmm: Option<String>,
    cached_template: Option<String>,
) -> Result<Vec<(String, String)>, String> {
    // Sanitise tag components so they're safe inside a filename — strip anything
    // that isn't alphanumeric/_/-/. so we don't accidentally inject path
    // separators or shell metacharacters into the tempfile prefix.
    fn slugify(s: Option<&str>) -> String {
        s.unwrap_or("")
            .chars()
            .filter(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.'))
            .collect()
    }
    let gene_slug = slugify(gene_name.as_deref());
    let taxa_slug = slugify(taxa.as_deref());
    let tag = match (taxa_slug.is_empty(), gene_slug.is_empty()) {
        (true, true) => String::new(),
        (true, false) => format!("{}_", gene_slug),
        (false, true) => format!("{}_", taxa_slug),
        (false, false) => format!("{}_{}_", taxa_slug, gene_slug),
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

    // A cache entry is only usable as a pair; a half-populated cache dir falls
    // back to building rather than handing bathalign a checksum it will reject.
    let cache_hit = match (&cached_hmm, &cached_template) {
        (Some(h), Some(t)) if Path::new(h).is_file() && Path::new(t).is_file() => {
            Some((h.clone(), t.clone()))
        }
        _ => None,
    };

    let mut temp_cand = make_temp("cand", ".fa")?;
    let temp_result = make_temp("res", ".afa")?;

    // Write through the existing tempfile handles - no second open().
    // Candidate `*`s are preserved so downstream stages still see the stop
    // signal; only the reference template masks them (see below).
    //
    // write_fasta_to_writer builds the complete Vec<u8> internally and emits
    // it in a single write_all, so wrapping the file in a BufWriter would
    // just add another copy-and-flush layer for no benefit.
    write_fasta_to_writer(temp_cand.as_file_mut(), &candidates, false)
        .map_err(|e| e.to_string())?;

    let cand_path = temp_cand.path().to_str().unwrap().to_string();
    let result_path = temp_result.path().to_str().unwrap().to_string();

    let hmm_name = if gene_slug.is_empty() { "hmm" } else { gene_slug.as_str() };

    // These two keep the scratch files alive for the length of the call; on a
    // cache hit neither is created and the cached paths are used instead.
    let mut _aln_guard = None;
    let mut _hmm_guard = None;

    let (aln_path, hmm_path) = match cache_hit {
        Some((cached_hmm_path, cached_template_path)) => (cached_template_path, cached_hmm_path),
        None => {
            let mut temp_aln = make_temp("aln", ".fa")?;
            let temp_hmm = make_temp("hmm", ".bhmm")?;

            // Mask stop codons (`*` -> `X`) in references only; the model is
            // built from these and bathbuild won't accept `*`.
            write_fasta_to_writer(temp_aln.as_file_mut(), &references, true)
                .map_err(|e| e.to_string())?;

            let aln_path = temp_aln.path().to_str().unwrap().to_string();
            let hmm_path = temp_hmm.path().to_str().unwrap().to_string();

            // bathbuild
            //
            // --cpu 1 keeps each invocation single-threaded; SAPPHYRE drives
            // this function from a multiprocessing pool, so a default of 2
            // pthreads per call would oversubscribe (N_workers x 2 threads).
            //
            // --informat afa is required: bathbuild refuses to guess between
            // aligned and unaligned FASTA and errors out rather than picking.
            //
            // The E-value calibration fits (200 sampled sequences each by
            // default) are ~18% of build runtime and only populate the STATS
            // lines. Alignment never reads those, and this model is built here,
            // consumed by the bathalign below and then discarded, so trim the
            // fits. 25 still fits the Gumbel/exponential tails; 1 fails
            // outright with "failed to determine msv mu".  Do NOT copy this
            // trim to a model the bhmm cache will hand to bathsearch.
            let mut bathbuild = Command::new("bathbuild");
            bathbuild
                .args(["--informat", "afa"])
                .args(["-n", hmm_name])
                .args(["--cpu", "1"])
                .args(["--EmN", "25", "--EvN", "25", "--EfN", "25"])
                .arg(&hmm_path)
                .arg(&aln_path);
            run_command(&mut bathbuild, "bathbuild")?;

            _aln_guard = Some(temp_aln);
            _hmm_guard = Some(temp_hmm);
            (aln_path, hmm_path)
        }
    };

    // bathalign --mapali (single-threaded, like HMMER3 hmmalign — no --cpu).
    let mut bathalign = Command::new("bathalign");
    bathalign.args([
        "--mapali",
        &aln_path,
        "--outformat",
        "afa",
        "-o",
        &result_path,
        &hmm_path,
        &cand_path,
    ]);
    run_command(&mut bathalign, "bathalign")?;

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
