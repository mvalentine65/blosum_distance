use anyhow::{bail, Context, Result};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use seq_io::fastq::{Reader as FastqReader, Record as FastqRecord};
use seq_io::fasta::{Reader as FastaReader, Record as FastaRecord};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use flate2::read::MultiGzDecoder;
// Same mapping this file used to define inline, but table-driven: it runs once
// per base of every read through `forward_is_canonical`.
use crate::reads::ingest::{group_inputs, Ingest, InputGroup};
use crate::reads::pipeline::{TrimOptions, TrimStats};
use crate::reads::seqops::complement;

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
    grows: usize,
    grow_time: std::time::Duration,
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
            grows: 0,
            grow_time: std::time::Duration::ZERO,
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
        let _t = std::time::Instant::now();
        self.grows += 1;
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
        self.grow_time += _t.elapsed();
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

/// Phase 1: read every input and fill the dedup table. Shared so `fast_dedupe`
/// and `dedupe_reads` can never disagree about ids.
///
/// With `trim` on, every record first goes through the fastp-derived chain in
/// [`crate::reads`]: adapters are detected from the head of the input and then
/// trimmed, quality filters drop what fails, and only the survivors are
/// canonicalised and hashed. That ordering is the point — trimming is what lets
/// two copies of the same fragment carrying different adapter remnants collapse
/// onto one entry, which they cannot do once hashed.
fn build_table(input_paths: Vec<PathBuf>, trim: bool) -> Result<(DedupTable, TrimSummary)> {
    let prof = std::env::var("SAPPHYRE_PROFILE").is_ok();
    let mut t_read = std::time::Duration::ZERO;
    let mut t_trim = std::time::Duration::ZERO;
    let mut t_hash = std::time::Duration::ZERO;
    let t_all = std::time::Instant::now();
    let (expected_uniques, expected_bytes) = size_hint(&input_paths);
    let mut table = DedupTable::with_capacity(expected_uniques, expected_bytes);

    if !trim {
        let mut scratch = Vec::with_capacity(300);
        for path in input_paths {
            let mut src = RecordSource::open(&path)?;
            loop {
                let t0 = std::time::Instant::now();
                let rec = src.next_record()?;
                if prof { t_read += t0.elapsed(); }
                let Some((seq, _qual)) = rec else { break };
                let t1 = std::time::Instant::now();
                add_canonical(&mut table, &seq, &mut scratch);
                if prof { t_hash += t1.elapsed(); }
            }
        }
        if prof {
            eprintln!(
                "[profile] no-trim total {:.1}s | read+decompress {:.1}s | canonicalise+hash {:.1}s",
                t_all.elapsed().as_secs_f64(), t_read.as_secs_f64(), t_hash.as_secs_f64());
        }
        return Ok((table, TrimSummary::default()));
    }

    let mut ingest = Ingest::new(TrimOptions::default());
    let mut scratch = Vec::with_capacity(300);

    for group in group_inputs(&input_paths) {
        match group {
            InputGroup::Paired(p1, p2) => {
                // Read the mates in lockstep so the overlap analyser can see
                // both halves of a fragment.
                let mut a = RecordSource::open(&p1)?;
                let mut b = RecordSource::open(&p2)?;
                loop {
                    let t0 = std::time::Instant::now();
                    let pair = (a.next_record()?, b.next_record()?);
                    if prof { t_read += t0.elapsed(); }
                    let (Some((s1, q1)), Some((s2, q2))) = pair else { break };
                    let t1 = std::time::Instant::now();
                    ingest.push_pair(&s1, &q1, &s2, &q2, &mut |kept: &[u8]| {
                        let t2 = std::time::Instant::now();
                        add_canonical(&mut table, kept, &mut scratch);
                        if prof { t_hash += t2.elapsed(); }
                    });
                    if prof { t_trim += t1.elapsed(); }
                }
            }
            InputGroup::Single(path) => {
                let mut src = RecordSource::open(&path)?;
                while let Some((seq, qual)) = src.next_record()? {
                    ingest.push_single(&seq, &qual, &mut |kept: &[u8]| {
                        add_canonical(&mut table, kept, &mut scratch)
                    });
                }
            }
        }
    }
    // Anything still staged awaiting detection has to be drained, or every
    // input shorter than the detection sample would vanish.
    ingest.finish(&mut |kept: &[u8]| add_canonical(&mut table, kept, &mut scratch));

    if prof {
        eprintln!(
            "[table] buckets {} ({:.2} GB) | uniques {} | load {:.1}% | grows {} ({:.1}s) | arena {:.2} GB of {:.2} GB reserved",
            table.buckets.len(),
            (table.buckets.len() * std::mem::size_of::<Bucket>()) as f64 / 1e9,
            table.next_id,
            100.0 * table.next_id as f64 / table.buckets.len() as f64,
            table.grows,
            table.grow_time.as_secs_f64(),
            table.arena.len() as f64 / 1e9,
            table.arena.capacity() as f64 / 1e9);
        eprintln!(
            "[profile] trim total {:.1}s | read+decompress {:.1}s | trim chain {:.1}s (of which hash {:.1}s)",
            t_all.elapsed().as_secs_f64(), t_read.as_secs_f64(),
            (t_trim - t_hash).as_secs_f64(), t_hash.as_secs_f64());
    }
    let summary = TrimSummary::from(ingest.stats(), ingest.detected());
    Ok((table, summary))
}

/// Canonicalise and hash one surviving sequence.
#[inline]
fn add_canonical(table: &mut DedupTable, seq: &[u8], scratch: &mut Vec<u8>) {
    if seq.is_empty() {
        return;
    }
    if forward_is_canonical(seq) {
        table.add(seq);
    } else {
        scratch.clear();
        scratch.extend(seq.iter().rev().map(|&b| complement(b)));
        table.add(scratch);
    }
}

/// A FASTQ or FASTA file yielding `(sequence, quality)`; quality is empty for
/// FASTA, which makes every quality-dependent stage a no-op downstream.
enum RecordSource {
    Fastq(FastqReader<Box<dyn Read + Send>>),
    Fasta(FastaReader<Box<dyn Read + Send>>),
}

impl RecordSource {
    fn open(path: &PathBuf) -> Result<Self> {
        Ok(match detect_format(path)? {
            InputFormat::Fastq => RecordSource::Fastq(FastqReader::new(open_reader(path)?)),
            InputFormat::Fasta => RecordSource::Fasta(FastaReader::new(open_reader(path)?)),
        })
    }

    #[allow(clippy::type_complexity)]
    fn next_record(&mut self) -> Result<Option<(Vec<u8>, Vec<u8>)>> {
        match self {
            RecordSource::Fastq(r) => match r.next() {
                None => Ok(None),
                Some(rec) => {
                    let rec = rec.context("Error reading FASTQ record")?;
                    Ok(Some((rec.seq().to_vec(), rec.qual().to_vec())))
                }
            },
            RecordSource::Fasta(r) => match r.next() {
                None => Ok(None),
                Some(rec) => {
                    let rec = rec.context("Error reading FASTA record")?;
                    Ok(Some((rec.seq().to_vec(), Vec::new())))
                }
            },
        }
    }
}

/// Stream one file, handing each record to `f`. Used by the untrimmed path,
/// which has no need to pair anything up.
fn for_each_record(
    path: &PathBuf,
    mut f: impl FnMut(&[u8], &[u8]) -> Result<()>,
) -> Result<()> {
    let mut src = RecordSource::open(path)?;
    while let Some((seq, qual)) = src.next_record()? {
        f(&seq, &qual)?;
    }
    Ok(())
}

/// What trimming did, for prepare's summary line.
#[derive(Default, Clone)]
pub struct TrimSummary {
    pub reads_in: u64,
    pub reads_passed: u64,
    pub adapter_trimmed_reads: u64,
    pub adapter_trimmed_bases: u64,
    pub failed_quality: u64,
    pub failed_n_base: u64,
    pub failed_length: u64,
    pub failed_adapter_dimer: u64,
    pub adapter_r1: Option<String>,
    pub adapter_r2: Option<String>,
}

impl TrimSummary {
    fn from(stats: &TrimStats, detected: (Option<&[u8]>, Option<&[u8]>)) -> Self {
        TrimSummary {
            reads_in: stats.reads_in,
            reads_passed: stats.reads_passed,
            adapter_trimmed_reads: stats.adapter_trimmed_reads,
            adapter_trimmed_bases: stats.adapter_trimmed_bases,
            failed_quality: stats.failed_quality,
            failed_n_base: stats.failed_n_base,
            failed_length: stats.failed_length,
            failed_adapter_dimer: stats.failed_adapter_dimer,
            adapter_r1: detected.0.map(|s| String::from_utf8_lossy(s).into_owned()),
            adapter_r2: detected.1.map(|s| String::from_utf8_lossy(s).into_owned()),
        }
    }
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
/// Python, marshalled every record back through Rust to split on N, and
/// reassembled the batches a record at a time.
#[pyclass]
pub struct PreparedReads {
    arena: Vec<u8>,
    trim: TrimSummary,
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

    // --- Trimming summary. All zero when trimming was switched off. ---

    /// Reads read in, before any filtering.
    fn reads_in(&self) -> u64 {
        self.trim.reads_in
    }

    /// Reads that survived trimming and filtering.
    fn reads_passed(&self) -> u64 {
        self.trim.reads_passed
    }

    fn adapter_trimmed_reads(&self) -> u64 {
        self.trim.adapter_trimmed_reads
    }

    fn adapter_trimmed_bases(&self) -> u64 {
        self.trim.adapter_trimmed_bases
    }

    /// Dropped reads by reason: quality, N content, length, adapter dimer.
    fn filtered_counts(&self) -> (u64, u64, u64, u64) {
        (
            self.trim.failed_quality,
            self.trim.failed_n_base,
            self.trim.failed_length,
            self.trim.failed_adapter_dimer,
        )
    }

    /// Adapters detection settled on, or None when nothing convincing was
    /// found -- which is the expected answer for an already-clean library.
    fn detected_adapters(&self) -> (Option<String>, Option<String>) {
        (self.trim.adapter_r1.clone(), self.trim.adapter_r2.clone())
    }
}

/// Dedupe the inputs and return the records the sequences db should hold.
///
/// Applies the length filter and N-splitting prepare used to run as its own
/// pass, plus the residue validation it did per record, so an illegal character
/// still names the offending read rather than surfacing later as a bathsearch
/// parse error.
#[pyfunction]
pub fn dedupe_reads(
    inputs: Vec<String>,
    min_length: usize,
    batch_size: usize,
    trim: bool,
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

    let (table, trim_summary) =
        build_table(input_paths, trim).map_err(|e| PyValueError::new_err(e.to_string()))?;

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

        // Whole read when it holds no N, else the N-free runs, each kept only
        // at or above the minimum length.
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
        trim: trim_summary,
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
