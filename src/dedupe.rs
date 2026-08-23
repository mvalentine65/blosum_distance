use anyhow::{bail, Context, Result};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use seq_io::fastq::{Reader as FastqReader, Record as FastqRecord};
use seq_io::fasta::{Reader as FastaReader, Record as FastaRecord};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::error::Error;
use flate2::read::MultiGzDecoder;

// --- Optimized DedupTable from main_bestrs.rs ---

#[derive(Clone, Copy)]
struct Bucket {
    hash: u64,
    id: u32,
    len: u32,
}

/// `id` of an unoccupied slot. Doubles as the occupancy flag.
const EMPTY: u32 = u32::MAX;

impl Default for Bucket {
    fn default() -> Self {
        Bucket { hash: 0, id: EMPTY, len: 0 }
    }
}

/// Per-unique data, indexed by id.
#[derive(Clone, Copy, Default)]
struct Record {
    off: u64,
    count: u32,
}

struct DedupTable {
    buckets: Vec<Bucket>,
    mask: usize,
    arena: Vec<u8>,
    used: usize,
    records: Vec<Record>,
    next_id: usize,
    /// Bucket slot holding each id, so the writer can emit in insertion order
    /// without building an inverse index over every bucket afterwards.
    id_to_slot: Vec<u32>,
}

impl DedupTable {
    fn with_capacity(expected_uniques: usize, expected_bytes: usize) -> Self {
        let mut cap = (expected_uniques * 10) / 7;
        cap = cap.next_power_of_two().max(1024);
        Self {
            buckets: vec![Bucket::default(); cap],
            mask: cap - 1,
            arena: Vec::with_capacity(expected_bytes.max(64 * 1024 * 1024)),
            used: 0,
            records: Vec::with_capacity(expected_uniques),
            next_id: 0,
            id_to_slot: Vec::with_capacity(expected_uniques),
        }
    }

    #[inline]
    fn arena_get(&self, off: usize, len: usize) -> &[u8] {
        &self.arena[off..off + len]
    }

    fn maybe_grow(&mut self) {
        if (self.used * 10) < (self.buckets.len() * 7) { return; }
        let old = std::mem::take(&mut self.buckets);
        let new_cap = old.len() * 2;
        self.buckets = vec![Bucket::default(); new_cap];
        self.mask = new_cap - 1;
        self.used = 0;
        for b in old.into_iter().filter(|b| b.id != EMPTY) {
            let mut i = (b.hash as usize) & self.mask;
            while self.buckets[i].id != EMPTY {
                i = (i + 1) & self.mask;
            }
            self.id_to_slot[b.id as usize] = i as u32;
            self.buckets[i] = b;
            self.used += 1;
        }
    }

    #[inline]
    fn hash64(seq: &[u8]) -> u64 {
        use std::hash::Hasher;
        let mut h = ahash::AHasher::default();
        h.write(seq);
        h.finish()
    }

    fn add(&mut self, seq: &[u8]) {
        self.maybe_grow();
        let h = Self::hash64(seq);
        let len = seq.len() as u32;
        let mut i = (h as usize) & self.mask;
        loop {
            let bucket = self.buckets[i];
            if bucket.id == EMPTY {
                let off = self.arena.len() as u64;
                self.arena.extend_from_slice(seq);
                let id = self.next_id as u32;
                self.next_id += 1;
                self.id_to_slot.push(i as u32);
                self.records.push(Record { off, count: 1 });
                self.buckets[i] = Bucket { hash: h, id, len };
                self.used += 1;
                return;
            }
            // hash and length reject almost every collision without touching
            // records or the arena.
            if bucket.hash == h && bucket.len == len {
                let rec = self.records[bucket.id as usize];
                if self.arena_get(rec.off as usize, len as usize) == seq {
                    self.records[bucket.id as usize].count += 1;
                    return;
                }
            }
            i = (i + 1) & self.mask;
        }
    }
}

// --- Canonicalization Helpers ---

#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        b'N' | b'n' => b'N',
        _ => b'N',
    }
}

#[inline]
fn forward_is_canonical(seq: &[u8]) -> bool {
    let n = seq.len();
    for i in 0..n {
        let f = seq[i];
        let r = complement(seq[n - 1 - i]);
        if f < r { return true; }
        else if f > r { return false; }
    }
    true
}

fn open_reader(path: &PathBuf) -> Result<Box<dyn Read + Send>> {
    let f = File::open(path).with_context(|| format!("Failed to open {}", path.display()))?;
    let br = BufReader::with_capacity(8 << 20, f);
    if path.extension().and_then(|e| e.to_str()).map(|e| e.eq_ignore_ascii_case("gz")).unwrap_or(false) {
        Ok(Box::new(MultiGzDecoder::new(br)))
    } else {
        Ok(Box::new(br))
    }
}

#[derive(Clone, Copy)]
enum InputFormat {
    Fasta,
    Fastq,
}

fn data_extension(path: &PathBuf) -> Option<&str> {
    let ext = path.extension().and_then(|e| e.to_str());
    if ext.map(|e| e.eq_ignore_ascii_case("gz")).unwrap_or(false) {
        path.file_stem()
            .and_then(|s| Path::new(s).extension())
            .and_then(|e| e.to_str())
    } else {
        ext
    }
}

fn detect_format(path: &PathBuf) -> Result<InputFormat> {
    let ext = data_extension(path)
        .map(|e| e.to_ascii_lowercase())
        .unwrap_or_default();
    match ext.as_str() {
        "fq" | "fastq" => Ok(InputFormat::Fastq),
        "fa" | "fasta" | "fna" | "fas" | "fsa_nt" => Ok(InputFormat::Fasta),
        _ => bail!("Unsupported input extension for {} (expected FASTQ: .fq/.fastq or FASTA: .fa/.fasta/.fna/.fas/.fsa_nt)", path.display()),
    }
}


/// Rough record count and sequence-byte estimate for the inputs, used only to
/// size the table and arena. Growing into 112M uniques from the old 2M default
/// cost six full rehashes and repeated multi-GB arena copies; being wrong here
/// just restores that behaviour, so the estimate is deliberately cheap.
fn size_hint(paths: &[PathBuf]) -> (usize, usize) {
    const GZ_RATIO: f64 = 3.5;
    let mut records = 0usize;
    for path in paths {
        let Ok(meta) = std::fs::metadata(path) else { continue };
        let compressed = path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.eq_ignore_ascii_case("gz"))
            .unwrap_or(false);
        let bytes = meta.len() as f64 * if compressed { GZ_RATIO } else { 1.0 };
        // A fastq record carries the sequence twice plus two header lines; a
        // fasta record carries it once.
        let per_record = match detect_format(path) {
            Ok(InputFormat::Fastq) => 335.0,
            _ => 190.0,
        };
        records += (bytes / per_record) as usize;
    }
    (records.max(1_000_000), records.saturating_mul(150))
}

pub fn fast_dedupe(
    mut input_paths: Vec<PathBuf>, // Changed from r1, r2 to a Vec
    out: PathBuf,
    sort_by_size: bool,
    min_size: u64,
) -> Result<(), Box<dyn Error>> {
    // Stabilize input order by file name so IDs are deterministic across runs.
    input_paths.sort_by(|a, b| {
        let a_name = a.file_name().and_then(|s| s.to_str()).unwrap_or("");
        let b_name = b.file_name().and_then(|s| s.to_str()).unwrap_or("");
        a_name
            .cmp(b_name)
            .then_with(|| a.as_os_str().cmp(b.as_os_str()))
    });

    let table = build_table(input_paths)?;

    // PHASE 2: WRITING
    let mut writer: Box<dyn Write> = if out.as_os_str() == "-" {
        Box::new(BufWriter::with_capacity(8 << 20, io::stdout().lock()))
    } else {
        let f = File::create(&out).with_context(|| format!("Failed to create {}", out.display()))?;
        Box::new(BufWriter::with_capacity(8 << 20, f))
    };

    // Input order is preserved regardless of sort_by_size, to avoid
    // nondeterministic ties. id_to_slot already gives insertion order, so no
    // inverse index over the buckets is built.
    let _ = sort_by_size;
    for id in 0..table.next_id {
        let b = table.buckets[table.id_to_slot[id] as usize];
        let rec = table.records[id];
        if u64::from(rec.count) < min_size {
            continue;
        }
        writeln!(writer, ">{}|{}", id + 1, rec.count)?;
        writer.write_all(table.arena_get(rec.off as usize, b.len as usize))?;
        writeln!(writer)?;
    }
    
    writer.flush()?;
    Ok(())
}

/// Phase 1: read every input and fill the dedup table. Shared so `fast_dedupe`
/// and `dedupe_reads` can never disagree about ids.
fn build_table(input_paths: Vec<PathBuf>) -> Result<DedupTable> {
    let (expected_uniques, expected_bytes) = size_hint(&input_paths);
    let mut table = DedupTable::with_capacity(expected_uniques, expected_bytes);
    let mut scratch = Vec::with_capacity(300);

    for path in input_paths {
        match detect_format(&path)? {
            InputFormat::Fastq => {
                let mut reader = FastqReader::new(open_reader(&path)?);
                let mut records = reader.records();
                while let Some(record) = records.next() {
                    let rec = record.context("Error reading FASTQ record")?;
                    let seq = rec.seq();
                    if forward_is_canonical(seq) {
                        table.add(seq);
                    } else {
                        scratch.clear();
                        scratch.extend(seq.iter().rev().map(|&b| complement(b)));
                        table.add(&scratch);
                    }
                }
            }
            InputFormat::Fasta => {
                let mut reader = FastaReader::new(open_reader(&path)?);
                let mut records = reader.records();
                while let Some(record) = records.next() {
                    let rec = record.context("Error reading FASTA record")?;
                    let seq = rec.seq();
                    if forward_is_canonical(seq) {
                        table.add(seq);
                    } else {
                        scratch.clear();
                        scratch.extend(seq.iter().rev().map(|&b| complement(b)));
                        table.add(&scratch);
                    }
                }
            }
        }
    }
    Ok(table)
}

const fn valid_nt_table() -> [bool; 256] {
    let mut t = [false; 256];
    let v = b"ACGTUNRYSWKMBDHVacgtunryswkmbdhv";
    let mut i = 0;
    while i < v.len() {
        t[v[i] as usize] = true;
        i += 1;
    }
    t
}
static VALID_NT: [bool; 256] = valid_nt_table();

const DUPES_MAGIC: &[u8; 8] = b"SPKD1\0\0\0";

/// Deduped reads, filtered and chunked, ready for the sequences db.
///
/// Replaces the round trip that wrote a multi-GB temp fasta, re-parsed it in
/// Python, marshalled every record back through `preprocess_n_chunks`, and
/// reassembled the batches a record at a time.
#[pyclass]
pub struct PreparedReads {
    arena: Vec<u8>,
    /// (node id, arena offset, length) per emitted record, in id order.
    records: Vec<(u64, usize, u32)>,
    batch_size: usize,
    dupes_blob: Vec<u8>,
    total_dupes: u64,
}

#[pymethods]
impl PreparedReads {
    /// Number of records that survived the length filter and N-splitting.
    fn record_count(&self) -> usize {
        self.records.len()
    }

    /// Total duplicate observations, i.e. sum of (count - 1).
    fn total_dupes(&self) -> u64 {
        self.total_dupes
    }

    fn batch_count(&self) -> usize {
        self.records.len().div_ceil(self.batch_size)
    }

    /// Batch `index` as concatenated `>id\nSEQ\n` records.
    fn batch<'py>(&self, py: Python<'py>, index: usize) -> PyResult<Bound<'py, PyBytes>> {
        let start = index * self.batch_size;
        if start >= self.records.len() && !self.records.is_empty() {
            return Err(PyValueError::new_err(format!("batch {index} out of range")));
        }
        let end = ((index + 1) * self.batch_size).min(self.records.len());
        let mut out = Vec::with_capacity((end - start) * 160);
        for &(node_id, off, len) in &self.records[start..end] {
            out.push(b'>');
            out.extend_from_slice(node_id.to_string().as_bytes());
            out.push(b'\n');
            out.extend_from_slice(&self.arena[off..off + len as usize]);
            out.push(b'\n');
        }
        Ok(PyBytes::new(py, &out))
    }

    /// The packed duplicate counts, in the same layout `packed_dupes` writes.
    fn packed_dupes<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.dupes_blob)
    }
}

/// Dedupe the inputs and return the records the sequences db should hold.
///
/// Applies the same length filter and N-splitting as `preprocess_n_chunks`, and
/// the same residue validation prepare did per record, so an illegal character
/// still names the offending read rather than surfacing later as a bathsearch
/// parse error.
#[pyfunction]
pub fn dedupe_reads(
    inputs: Vec<String>,
    min_length: usize,
    batch_size: usize,
) -> PyResult<PreparedReads> {
    if batch_size == 0 {
        return Err(PyValueError::new_err("batch_size must be non-zero"));
    }
    let mut input_paths: Vec<PathBuf> = inputs.into_iter().map(PathBuf::from).collect();
    input_paths.sort_by(|a, b| {
        let a_name = a.file_name().and_then(|s| s.to_str()).unwrap_or("");
        let b_name = b.file_name().and_then(|s| s.to_str()).unwrap_or("");
        a_name.cmp(b_name).then_with(|| a.as_os_str().cmp(b.as_os_str()))
    });

    let table = build_table(input_paths).map_err(|e| PyValueError::new_err(e.to_string()))?;

    let mut records: Vec<(u64, usize, u32)> = Vec::with_capacity(table.next_id);
    let mut dupe_keys: Vec<i64> = Vec::with_capacity(table.next_id);
    let mut dupe_vals: Vec<i64> = Vec::with_capacity(table.next_id);
    let mut total_dupes: u64 = 0;

    for id in 0..table.next_id {
        let b = table.buckets[table.id_to_slot[id] as usize];
        let rec = table.records[id];
        let node_id = (id + 1) as u64;
        let off = rec.off as usize;
        let seq = &table.arena[off..off + b.len as usize];

        // Mirrors preprocess_n_chunks: whole read when it holds no N, else the
        // N-free runs, each kept only at or above the minimum length.
        let mut emitted = false;
        if !seq.contains(&b'N') && !seq.contains(&b'n') {
            if seq.len() >= min_length {
                validate_residues(node_id, seq)?;
                records.push((node_id, off, seq.len() as u32));
                emitted = true;
            }
        } else {
            let mut chunk_start = 0usize;
            for i in 0..=seq.len() {
                let is_end = i == seq.len();
                if is_end || seq[i] == b'N' || seq[i] == b'n' {
                    let len = i - chunk_start;
                    if len >= min_length {
                        let chunk_off = off + chunk_start;
                        validate_residues(node_id, &table.arena[chunk_off..chunk_off + len])?;
                        records.push((node_id, chunk_off, len as u32));
                        emitted = true;
                    }
                    chunk_start = i + 1;
                }
            }
        }

        // prepare recorded the dupe count once per surviving header.
        if emitted {
            let dupes = u64::from(rec.count) - 1;
            total_dupes += dupes;
            if dupes != 0 {
                dupe_keys.push(node_id as i64);
                dupe_vals.push(dupes as i64);
            }
        }
    }

    let mut dupes_blob =
        Vec::with_capacity(DUPES_MAGIC.len() + 8 + dupe_keys.len() * 16);
    dupes_blob.extend_from_slice(DUPES_MAGIC);
    dupes_blob.extend_from_slice(&(dupe_keys.len() as u64).to_le_bytes());
    for k in &dupe_keys {
        dupes_blob.extend_from_slice(&k.to_le_bytes());
    }
    for v in &dupe_vals {
        dupes_blob.extend_from_slice(&v.to_le_bytes());
    }

    Ok(PreparedReads {
        arena: table.arena,
        records,
        batch_size,
        dupes_blob,
        total_dupes,
    })
}

fn validate_residues(node_id: u64, seq: &[u8]) -> PyResult<()> {
    for (i, &byte) in seq.iter().enumerate() {
        if !VALID_NT[byte as usize] {
            let start = i.saturating_sub(20);
            let end = (i + 21).min(seq.len());
            let preview = String::from_utf8_lossy(&seq[start..end]).into_owned();
            return Err(PyValueError::new_err(format!(
                "{node_id}: nucleotide sequence contains character {:?} that is \
                 not a valid IUPAC nucleotide code and would make bathsearch \
                 fail. First at offset {i} of {}: ...{preview}... Check the \
                 input file for non-sequence data (e.g. a stray '|') on this \
                 record.",
                byte as char,
                seq.len()
            )));
        }
    }
    Ok(())
}
