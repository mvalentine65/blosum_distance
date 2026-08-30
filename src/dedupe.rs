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
use crate::reads::detect::{detect_adapter, Detected};
use crate::reads::ingest::{group_inputs, Ingest, InputGroup, DETECT_SAMPLE};
use crate::reads::periodicity::{canonical_unit, repeat_period};
use crate::reads::pipeline::{
    process_pair, process_single, TrimOptions, TrimScratch, TrimStats,
};
use crate::reads::seed::PreparedAdapter;
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

/// Per-unique data, indexed by id. `len` rides in the tail padding `off` and
/// `count` already forced, so the writer never has to go back to the buckets.
#[derive(Clone, Copy, Default)]
struct Record {
    off: u64,
    count: u32,
    len: u32,
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
        let h = Self::hash64(seq);
        self.add_hashed(seq, h);
    }

    /// `add` with the hash already computed. Lets a worker thread do the
    /// hashing off the critical path; the insert itself stays single-threaded
    /// because `next_id` is the encounter order the output depends on.
    fn add_hashed(&mut self, seq: &[u8], h: u64) {
        self.maybe_grow();
        let len = seq.len() as u32;
        let mut i = (h as usize) & self.mask;
        loop {
            let bucket = self.buckets[i];
            if bucket.id == EMPTY {
                let off = self.arena.len() as u64;
                self.arena.extend_from_slice(seq);
                let id = self.next_id as u32;
                self.next_id += 1;
                self.records.push(Record { off, count: 1, len });
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
/// One input record: a mate pair, or a single read with an empty mate.
struct Item {
    s1: Vec<u8>,
    q1: Vec<u8>,
    s2: Vec<u8>,
    q2: Vec<u8>,
}

/// Every input group as one ordered stream, so the encounter order the output
/// depends on is a property of the stream rather than of the reader loop.
struct ItemStream {
    groups: std::vec::IntoIter<InputGroup>,
    cur: Option<Cursor>,
}

enum Cursor {
    Paired(Box<RecordSource>, Box<RecordSource>),
    Single(Box<RecordSource>),
}

/// One mate's records, read on its own thread. The final batch is short (and
/// may be empty), which is how the zipper learns the file ended.
type RawBatch = Vec<(Vec<u8>, Vec<u8>)>;

fn read_batches(
    mut src: Box<RecordSource>,
    batch: usize,
    tx: std::sync::mpsc::SyncSender<Result<RawBatch>>,
) {
    let mut cur: RawBatch = Vec::with_capacity(batch);
    loop {
        match src.next_record() {
            Err(e) => {
                let _ = tx.send(Err(e));
                return;
            }
            Ok(None) => {
                let _ = tx.send(Ok(cur));
                return;
            }
            Ok(Some(rec)) => {
                cur.push(rec);
                if cur.len() == batch {
                    let full = std::mem::replace(&mut cur, Vec::with_capacity(batch));
                    if tx.send(Ok(full)).is_err() {
                        return;
                    }
                }
            }
        }
    }
}

/// Read a mate pair with one thread per file. Inflate is 84% of the reader's
/// cost and the two gz streams are independent, so this is where the reader's
/// time actually goes. The zipper pairs batch k of one mate with batch k of
/// the other, which keeps the pairing and the order identical to reading them
/// in lockstep on one thread.
fn read_pair_parallel<F>(
    a: Box<RecordSource>,
    b: Box<RecordSource>,
    batch: usize,
    emit: &mut F,
) -> Result<()>
where
    F: FnMut(Vec<Item>) -> bool,
{
    std::thread::scope(|s| -> Result<()> {
        let (txa, rxa) = std::sync::mpsc::sync_channel::<Result<RawBatch>>(2);
        let (txb, rxb) = std::sync::mpsc::sync_channel::<Result<RawBatch>>(2);
        s.spawn(move || read_batches(a, batch, txa));
        s.spawn(move || read_batches(b, batch, txb));
        while let (Ok(ra), Ok(rb)) = (rxa.recv(), rxb.recv()) {
            let (va, vb) = (ra?, rb?);
            // zip stops at the shorter, which is the lockstep loop's rule for
            // mates of unequal length.
            let short = va.len() < batch || vb.len() < batch;
            let items: Vec<Item> = va
                .into_iter()
                .zip(vb)
                .map(|((s1, q1), (s2, q2))| Item { s1, q1, s2, q2 })
                .collect();
            if !items.is_empty() && !emit(items) {
                return Ok(());
            }
            if short {
                break;
            }
        }
        Ok(())
    })
}

impl ItemStream {
    fn new(input_paths: &[PathBuf]) -> Self {
        ItemStream { groups: group_inputs(input_paths).into_iter(), cur: None }
    }

    /// Hand back the open cursor and the groups not started yet, so the reader
    /// can resume where detection stopped instead of reopening the files.
    fn into_parts(self) -> (Option<Cursor>, std::vec::IntoIter<InputGroup>) {
        (self.cur, self.groups)
    }

    fn next_item(&mut self) -> Result<Option<Item>> {
        loop {
            match &mut self.cur {
                // Read the mates in lockstep so the overlap analyser can see
                // both halves of a fragment.
                Some(Cursor::Paired(a, b)) => match (a.next_record()?, b.next_record()?) {
                    (Some((s1, q1)), Some((s2, q2))) => {
                        return Ok(Some(Item { s1, q1, s2, q2 }))
                    }
                    _ => self.cur = None,
                },
                Some(Cursor::Single(src)) => match src.next_record()? {
                    Some((s1, q1)) => {
                        return Ok(Some(Item { s1, q1, s2: Vec::new(), q2: Vec::new() }))
                    }
                    None => self.cur = None,
                },
                None => match self.groups.next() {
                    None => return Ok(None),
                    Some(InputGroup::Paired(p1, p2)) => {
                        self.cur = Some(Cursor::Paired(
                            Box::new(RecordSource::open(&p1)?),
                            Box::new(RecordSource::open(&p2)?),
                        ))
                    }
                    Some(InputGroup::Single(p)) => {
                        self.cur = Some(Cursor::Single(Box::new(RecordSource::open(&p)?)))
                    }
                },
            }
        }
    }
}

/// Survivors of one batch, packed back to back. `spans` carries each kept
/// sequence's length and hash in emission order, so the consumer can insert
/// them without rehashing.
#[derive(Default)]
struct Kept {
    buf: Vec<u8>,
    spans: Vec<(u32, u64)>,
}

/// Canonicalise and hash one survivor into `out`. Mirrors `add_canonical`,
/// including its empty-sequence early return, but defers the table insert.
#[inline]
fn push_kept(out: &mut Kept, seq: &[u8], canon: &mut Vec<u8>) {
    if seq.is_empty() {
        return;
    }
    let bytes: &[u8] = if forward_is_canonical(seq) {
        seq
    } else {
        canon.clear();
        canon.extend(seq.iter().rev().map(|&b| complement(b)));
        canon
    };
    out.spans.push((bytes.len() as u32, DedupTable::hash64(bytes)));
    out.buf.extend_from_slice(bytes);
}

/// The trim chain over one batch. Pure given `opts`: the only state is the
/// worker's own scratch and stats, which is what makes the batch shardable.
fn run_batch(
    items: &[Item],
    opts: &TrimOptions,
    stats: &mut TrimStats,
    scratch: &mut TrimScratch,
    canon: &mut Vec<u8>,
) -> Kept {
    let mut out = Kept::default();
    for it in items {
        if it.s2.is_empty() {
            if let Some(r) = process_single(&it.s1, &it.q1, opts, stats, scratch) {
                push_kept(&mut out, r.bases(), canon);
            }
        } else {
            let o = process_pair(&it.s1, &it.q1, &it.s2, &it.q2, opts, stats, scratch);
            if let Some(r) = o.r1 {
                push_kept(&mut out, r.bases(), canon);
            }
            if let Some(r) = o.r2 {
                push_kept(&mut out, r.bases(), canon);
            }
        }
    }
    out
}

/// Fold one worker's stats into the total. Every field is a running count or a
/// histogram, so the sum does not depend on which worker saw which read.
fn merge_stats(dst: &mut TrimStats, src: &TrimStats) {
    dst.reads_in += src.reads_in;
    dst.reads_passed += src.reads_passed;
    dst.adapter_trimmed_reads += src.adapter_trimmed_reads;
    dst.adapter_trimmed_bases += src.adapter_trimmed_bases;
    dst.polyx_trimmed_reads += src.polyx_trimmed_reads;
    dst.polyx_trimmed_bases += src.polyx_trimmed_bases;
    dst.failed_quality += src.failed_quality;
    dst.failed_n_base += src.failed_n_base;
    dst.failed_length += src.failed_length;
    dst.failed_too_long += src.failed_too_long;
    dst.failed_complexity += src.failed_complexity;
    dst.failed_adapter_dimer += src.failed_adapter_dimer;
    dst.corrected_bases += src.corrected_bases;
    dst.pairs_merged += src.pairs_merged;
    dst.bases_in += src.bases_in;
    dst.bases_out += src.bases_out;
    for (d, s) in [
        (&mut dst.length_hist, &src.length_hist),
        (&mut dst.insert_hist, &src.insert_hist),
    ] {
        if d.len() < s.len() {
            d.resize(s.len(), 0);
        }
        for (i, v) in s.iter().enumerate() {
            d[i] += v;
        }
    }
}

/// Trim path, sharded across threads.
///
/// The reader emits items in input order and deals batches round robin; the
/// consumer takes them back in the same round robin, so `next_id` is assigned
/// in exactly the order the serial path assigns it. Only the trim chain and
/// the hashing move off the main thread — the table insert stays serial
/// because the encounter order is the output.
fn build_table_parallel(
    input_paths: Vec<PathBuf>,
    table: &mut DedupTable,
    workers: usize,
    prof: bool,
) -> Result<TrimSummary> {
    const BATCH: usize = 4096;

    let mut stream = ItemStream::new(&input_paths);

    // Detection sees the same head of the stream as the serial path.
    let mut staged: Vec<Item> = Vec::new();
    while staged.len() < DETECT_SAMPLE {
        match stream.next_item()? {
            Some(it) => staged.push(it),
            None => break,
        }
    }
    let t_detect = std::time::Instant::now();
    let mut opts = TrimOptions::default();
    let (mut det1, mut det2) = (None, None);
    if opts.adapter_trimming {
        let sample1: Vec<&[u8]> = staged.iter().map(|p| p.s1.as_slice()).collect();
        if let Detected::Known { seq, .. } | Detected::Novel { seq } = detect_adapter(&sample1, 1)
        {
            opts.adapter_r1 = Some(PreparedAdapter::new(&seq));
            det1 = Some(seq);
        }
        if staged.first().is_some_and(|p| !p.s2.is_empty()) {
            let sample2: Vec<&[u8]> = staged.iter().map(|p| p.s2.as_slice()).collect();
            if let Detected::Known { seq, .. } | Detected::Novel { seq } =
                detect_adapter(&sample2, 1)
            {
                opts.adapter_r2 = Some(PreparedAdapter::new(&seq));
                det2 = Some(seq);
            }
        }
    }

    if prof {
        eprintln!("[detect] staged {} items, detection took {:.1}s",
                  staged.len(), t_detect.elapsed().as_secs_f64());
    }
    let mut totals = TrimStats::default();
    let opts = &opts;

    std::thread::scope(|scope| -> Result<()> {
        let mut to_worker = Vec::with_capacity(workers);
        let mut from_worker = Vec::with_capacity(workers);
        let mut handles = Vec::with_capacity(workers);

        for _ in 0..workers {
            let (tx_in, rx_in) = std::sync::mpsc::sync_channel::<Vec<Item>>(2);
            let (tx_out, rx_out) = std::sync::mpsc::sync_channel::<Kept>(2);
            to_worker.push(tx_in);
            from_worker.push(rx_out);
            handles.push(scope.spawn(move || {
                let mut stats = TrimStats::default();
                let mut scratch = TrimScratch::default();
                let mut canon = Vec::with_capacity(300);
                while let Ok(items) = rx_in.recv() {
                    let kept = run_batch(&items, opts, &mut stats, &mut scratch, &mut canon);
                    if tx_out.send(kept).is_err() {
                        break;
                    }
                }
                stats
            }));
        }

        // Reader: staged head first, then the rest of the input.
        let reader = scope.spawn(move || -> Result<()> {
            let mut seq = 0usize;
            let mut emit = |items: Vec<Item>| -> bool {
                if items.is_empty() {
                    return true;
                }
                let ok = to_worker[seq % workers].send(items).is_ok();
                seq += 1;
                ok
            };

            let mut staged = staged.into_iter();
            loop {
                let chunk: Vec<Item> = staged.by_ref().take(BATCH).collect();
                if chunk.is_empty() {
                    break;
                }
                if !emit(chunk) {
                    return Ok(());
                }
            }

            // Everything detection did not consume, resuming the cursor it
            // left open before moving on to any later groups.
            let (cur, groups) = stream.into_parts();
            let mut pending = cur;
            let mut groups = groups;
            loop {
                let cursor = match pending.take() {
                    Some(c) => c,
                    None => match groups.next() {
                        None => break,
                        Some(InputGroup::Paired(p1, p2)) => Cursor::Paired(
                            Box::new(RecordSource::open(&p1)?),
                            Box::new(RecordSource::open(&p2)?),
                        ),
                        Some(InputGroup::Single(p)) => {
                            Cursor::Single(Box::new(RecordSource::open(&p)?))
                        }
                    },
                };
                match cursor {
                    Cursor::Paired(a, b) => read_pair_parallel(a, b, BATCH, &mut emit)?,
                    Cursor::Single(mut src) => {
                        let mut batch = Vec::with_capacity(BATCH);
                        while let Some((s1, q1)) = src.next_record()? {
                            batch.push(Item { s1, q1, s2: Vec::new(), q2: Vec::new() });
                            if batch.len() == BATCH {
                                let full =
                                    std::mem::replace(&mut batch, Vec::with_capacity(BATCH));
                                if !emit(full) {
                                    return Ok(());
                                }
                            }
                        }
                        if !emit(batch) {
                            return Ok(());
                        }
                    }
                }
            }
            Ok(())
        });

        // Consumer: same round robin, so survivors reach the table in input
        // order and `next_id` matches the serial path exactly.
        let mut seq = 0usize;
        let (mut t_wait, mut t_insert) = (
            std::time::Duration::ZERO,
            std::time::Duration::ZERO,
        );
        loop {
            let t0 = prof.then(std::time::Instant::now);
            let Ok(kept) = from_worker[seq % workers].recv() else { break };
            if let Some(t) = t0 {
                t_wait += t.elapsed();
            }
            let t1 = prof.then(std::time::Instant::now);
            let mut pos = 0usize;
            for &(len, h) in &kept.spans {
                let end = pos + len as usize;
                table.add_hashed(&kept.buf[pos..end], h);
                pos = end;
            }
            if let Some(t) = t1 {
                t_insert += t.elapsed();
            }
            seq += 1;
        }
        if prof {
            eprintln!(
                "[parallel] consumer blocked on recv {:.1}s | inserting {:.1}s",
                t_wait.as_secs_f64(),
                t_insert.as_secs_f64()
            );
        }

        for h in handles {
            if let Ok(stats) = h.join() {
                merge_stats(&mut totals, &stats);
            }
        }
        reader.join().unwrap_or(Ok(()))
    })?;

    Ok(TrimSummary::from(&totals, (det1.as_deref(), det2.as_deref())))
}

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
                let t0 = prof.then(std::time::Instant::now);
                let rec = src.next_record()?;
                if let Some(t) = t0 { t_read += t.elapsed(); }
                let Some((seq, _qual)) = rec else { break };
                let t1 = prof.then(std::time::Instant::now);
                add_canonical(&mut table, &seq, &mut scratch);
                if let Some(t) = t1 { t_hash += t.elapsed(); }
            }
        }
        if prof {
            eprintln!(
                "[profile] no-trim total {:.1}s | read+decompress {:.1}s | canonicalise+hash {:.1}s",
                t_all.elapsed().as_secs_f64(), t_read.as_secs_f64(), t_hash.as_secs_f64());
        }
        return Ok((table, TrimSummary::default()));
    }

    // Throughput plateaus once the workers can keep the serial insert fed,
    // which measured at three on DDM8637; past six only peak RSS moves, so the
    // cap is deliberate rather than a core count.
    //
    // SAPPHYRE_THREADS overrides it. A caller running several taxa at once
    // should divide the machine between them -- one worker still pipelines the
    // reader and the insert, so dividing down stays on this path. Zero selects
    // the serial path below, which is the reference the parallel one is
    // validated against.
    const MAX_WORKERS: usize = 6;
    let workers = std::env::var("SAPPHYRE_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or_else(|| {
            std::thread::available_parallelism()
                .map(|n| n.get().saturating_sub(2).min(MAX_WORKERS))
                .unwrap_or(1)
        });

    if workers >= 1 {
        let summary = build_table_parallel(input_paths, &mut table, workers, prof)?;
        if prof {
            eprintln!(
                "[profile] parallel trim total {:.1}s | {} workers",
                t_all.elapsed().as_secs_f64(),
                workers
            );
        }
        return Ok((table, summary));
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
                    let t0 = prof.then(std::time::Instant::now);
                    let pair = (a.next_record()?, b.next_record()?);
                    if let Some(t) = t0 { t_read += t.elapsed(); }
                    let (Some((s1, q1)), Some((s2, q2))) = pair else { break };
                    let t1 = prof.then(std::time::Instant::now);
                    ingest.push_pair(&s1, &q1, &s2, &q2, &mut |kept: &[u8]| {
                        let t2 = prof.then(std::time::Instant::now);
                        add_canonical(&mut table, kept, &mut scratch);
                        if let Some(t) = t2 { t_hash += t.elapsed(); }
                    });
                    if let Some(t) = t1 { t_trim += t.elapsed(); }
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

/// How many repeat units the summary carries. Far more than a report shows;
/// the tail stays visible as the gap between the listed counts and the total.
const REPEAT_UNITS_REPORTED: usize = 50;

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
    pub bases_in: u64,
    pub bases_out: u64,
    pub polyx_trimmed_reads: u64,
    pub polyx_trimmed_bases: u64,
    pub length_hist: Vec<u64>,
    pub insert_hist: Vec<u64>,
    pub repeat_sequences: u64,
    pub repeat_reads: u64,
    pub repeat_bases: u64,
    /// Repeat unit and count, commonest first, capped at the head of the tail.
    pub repeat_units: Vec<(String, u64)>,
    pub repeat_units_distinct: u64,
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
            bases_in: stats.bases_in,
            bases_out: stats.bases_out,
            polyx_trimmed_reads: stats.polyx_trimmed_reads,
            polyx_trimmed_bases: stats.polyx_trimmed_bases,
            length_hist: stats.length_hist.clone(),
            insert_hist: stats.insert_hist.clone(),
            repeat_sequences: 0,
            repeat_reads: 0,
            repeat_bases: 0,
            repeat_units: Vec::new(),
            repeat_units_distinct: 0,
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

    /// Bases before trimming and after it.
    fn base_counts(&self) -> (u64, u64) {
        (self.trim.bases_in, self.trim.bases_out)
    }

    /// polyG/polyX tails removed: (reads, bases).
    fn polyx_trimmed(&self) -> (u64, u64) {
        (self.trim.polyx_trimmed_reads, self.trim.polyx_trimmed_bases)
    }

    /// Length of every surviving read, indexed by length.
    fn length_histogram(&self) -> Vec<u64> {
        self.trim.length_hist.clone()
    }

    /// Tandem repeats dropped: (sequences, reads they represented, bases).
    fn repeats_removed(&self) -> (u64, u64, u64) {
        (
            self.trim.repeat_sequences,
            self.trim.repeat_reads,
            self.trim.repeat_bases,
        )
    }

    /// Commonest repeat units, and how many distinct units were seen. Folded
    /// across rotation and strand, so `(AG)n` and `(CT)n` are one satellite.
    fn repeat_units(&self) -> (Vec<(String, u64)>, u64) {
        (self.trim.repeat_units.clone(), self.trim.repeat_units_distinct)
    }

    /// Fragment length per pair, indexed by size. The final bucket holds pairs
    /// whose overlap could not be placed. Empty for unpaired input.
    fn insert_size_histogram(&self) -> Vec<u64> {
        self.trim.insert_hist.clone()
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
    repeat_filter: bool,
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
    // Serialised as they are found; the blob needs every key before every value,
    // which is the only reason both runs are held at all.
    let mut dupe_key_bytes: Vec<u8> = Vec::new();
    let mut dupe_val_bytes: Vec<u8> = Vec::new();
    let mut dupe_count: u64 = 0;
    let mut total_dupes: u64 = 0;

    // Ids are per emitted record, not per unique sequence: a read split at an
    // N can yield several runs, and downstream keys them in a dict by id.
    let mut next_node_id: u64 = 0;
    let mut runs: Vec<(usize, usize)> = Vec::new();
    let mut repeat_sequences: u64 = 0;
    let mut repeat_reads: u64 = 0;
    let mut repeat_bases: u64 = 0;
    let mut repeat_units: ahash::AHashMap<Vec<u8>, u64> = ahash::AHashMap::new();

    for id in 0..table.next_id {
        let rec = table.records[id];
        let off = rec.off as usize;
        let seq = &table.arena[off..off + rec.len as usize];

        // Whole read when it holds no N, else the N-free runs, each kept only
        // at or above the minimum length.
        runs.clear();
        if !seq.contains(&b'N') && !seq.contains(&b'n') {
            if seq.len() >= min_length {
                runs.push((off, seq.len()));
            }
        } else {
            let mut chunk_start = 0usize;
            for i in 0..=seq.len() {
                let is_end = i == seq.len();
                if is_end || seq[i] == b'N' || seq[i] == b'n' {
                    let len = i - chunk_start;
                    if len >= min_length {
                        runs.push((off + chunk_start, len));
                    }
                    chunk_start = i + 1;
                }
            }
        }

        // Dropped here, not in the trim chain: once per unique sequence rather
        // than per read, and only on runs that already cleared the floor.
        if repeat_filter {
            runs.retain(|&(off, len)| {
                let seq = &table.arena[off..off + len];
                let hit = repeat_period(
                    seq,
                    crate::reads::periodicity::DEFAULT_THRESHOLD,
                );
                let Some(period) = hit else {
                    return true;
                };
                repeat_sequences += 1;
                repeat_reads += u64::from(rec.count);
                repeat_bases += len as u64;
                *repeat_units.entry(canonical_unit(seq, period)).or_insert(0) += 1;
                false
            });
        }

        if runs.is_empty() {
            continue;
        }

        // Duplicates are identical reads, so they split identically and each
        // run inherits the parent's count. The total still counts them once.
        let dupes = u64::from(rec.count) - 1;
        total_dupes += dupes;

        for &(chunk_off, len) in &runs {
            next_node_id += 1;
            validate_residues(next_node_id, &table.arena[chunk_off..chunk_off + len])?;
            records.push((next_node_id, chunk_off, len as u32));
            if dupes != 0 {
                dupe_key_bytes.extend_from_slice(&(next_node_id as i64).to_le_bytes());
                dupe_val_bytes.extend_from_slice(&(dupes as i64).to_le_bytes());
                dupe_count += 1;
            }
        }
    }

    let mut dupes_blob = Vec::with_capacity(
        DUPES_MAGIC.len() + 8 + dupe_key_bytes.len() + dupe_val_bytes.len(),
    );
    dupes_blob.extend_from_slice(DUPES_MAGIC);
    dupes_blob.extend_from_slice(&dupe_count.to_le_bytes());
    dupes_blob.extend_from_slice(&dupe_key_bytes);
    dupes_blob.extend_from_slice(&dupe_val_bytes);

    // Commonest first, then alphabetically so ties are stable across runs.
    let mut units: Vec<(Vec<u8>, u64)> = repeat_units.into_iter().collect();
    units.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let units_distinct = units.len() as u64;
    units.truncate(REPEAT_UNITS_REPORTED);

    let mut trim_summary = trim_summary;
    trim_summary.repeat_sequences = repeat_sequences;
    trim_summary.repeat_reads = repeat_reads;
    trim_summary.repeat_bases = repeat_bases;
    trim_summary.repeat_units_distinct = units_distinct;
    trim_summary.repeat_units = units
        .into_iter()
        .map(|(u, n)| (String::from_utf8_lossy(&u).into_owned(), n))
        .collect();

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
