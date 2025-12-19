use anyhow::{bail, Context, Result};
use crossbeam_channel as chan;
use rayon::prelude::*;
use seq_io::fastq::{Reader, Record};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::error::Error;
use flate2::read::MultiGzDecoder;

// --- DedupTable Structures ---

#[derive(Clone, Copy)]
struct Bucket {
    hash: u64,
    off: u32,
    len: u32,
    count: u64,
    used: bool,
}

impl Default for Bucket {
    fn default() -> Self {
        Self {
            hash: 0,
            off: 0,
            len: 0,
            count: 0,
            used: false,
        }
    }
}

struct DedupTable {
    buckets: Vec<Bucket>,
    mask: u64,
    arena: Vec<u8>,
    used: usize,
}

impl DedupTable {
    fn with_capacity(expected_uniques: usize, arena_bytes: usize) -> Self {
        let mut cap = (expected_uniques * 10) / 7;
        cap = cap.next_power_of_two().max(1024);
        Self {
            buckets: vec![Bucket::default(); cap],
            mask: (cap as u64) - 1,
            arena: Vec::with_capacity(arena_bytes),
            used: 0,
        }
    }

    #[inline]
    fn arena_get(&self, off: u32, len: u32) -> &[u8] {
        let o = off as usize;
        let l = len as usize;
        &self.arena[o..o + l]
    }

    #[inline]
    fn maybe_grow(&mut self) {
        if (self.used * 10) < (self.buckets.len() * 7) {
            return;
        }
        let old = std::mem::take(&mut self.buckets);
        let new_cap = old.len() * 2;
        self.buckets = vec![Bucket::default(); new_cap];
        self.mask = (new_cap as u64) - 1;
        self.used = 0;

        for b in old.into_iter().filter(|b| b.used) {
            self.insert_existing(b);
        }
    }

    #[inline]
    fn insert_existing(&mut self, b: Bucket) {
        let mut i = (b.hash & self.mask) as usize;
        loop {
            let slot = &mut self.buckets[i];
            if !slot.used {
                *slot = b;
                self.used += 1;
                return;
            }
            i = (i + 1) & (self.mask as usize);
        }
    }

    #[inline]
    fn hash64(seq: &[u8]) -> u64 {
        use std::hash::Hasher;
        let mut h = ahash::AHasher::default();
        h.write(seq);
        h.finish()
    }

    fn add_with_count(&mut self, seq: &[u8], add: u64) {
        self.maybe_grow();
        let h = Self::hash64(seq);
        let len = seq.len() as u32;
        let mut i = (h & self.mask) as usize;
        loop {
            if !self.buckets[i].used {
                let off = self.arena.len() as u32;
                self.arena.extend_from_slice(seq);
                self.buckets[i] = Bucket {
                    hash: h,
                    off,
                    len,
                    count: add,
                    used: true,
                };
                self.used += 1;
                return;
            }

            let (bh, bl, bo) = {
                let b = self.buckets[i];
                (b.hash, b.len, b.off)
            };

            if bh == h && bl == len {
                let existing = self.arena_get(bo, bl);
                if existing == seq {
                    self.buckets[i].count += add;
                    return;
                }
            }
            i = (i + 1) & (self.mask as usize);
        }
    }

    #[inline]
    fn add(&mut self, seq: &[u8]) {
        self.add_with_count(seq, 1);
    }

    fn iter_buckets(&self) -> impl Iterator<Item = Bucket> + '_ {
        self.buckets.iter().copied().filter(|b| b.used)
    }
}

fn merge_tables(mut a: DedupTable, b: DedupTable) -> DedupTable {
    for buck in b.iter_buckets() {
        let seq = b.arena_get(buck.off, buck.len);
        a.add_with_count(seq, buck.count);
    }
    a
}

// --- Helper Logic ---

#[derive(Debug)]
struct Batch {
    data: Vec<u8>,
    spans: Vec<(usize, usize)>,
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

#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A',
        b'a' => b'T', b'c' => b'G', b'g' => b'C', b't' => b'A',
        b'N' | b'n' => b'N', _ => b'N',
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

#[inline]
fn append_canonical_into(batch: &mut Batch, seq: &[u8]) -> (usize, usize) {
    let start = batch.data.len();
    if forward_is_canonical(seq) {
        batch.data.extend_from_slice(seq);
    } else {
        batch.data.reserve(seq.len());
        for &b in seq.iter().rev() {
            batch.data.push(complement(b));
        }
    }
    let end = batch.data.len();
    (start, end)
}

// --- Main Entry Point ---

pub fn fast_dedupe(
    r1: PathBuf,
    r2: PathBuf,
    out: PathBuf,
    threads: usize,
    batch_pairs: usize,
    sort_by_size: bool,
    min_size: usize,
) -> Result<(), Box<dyn Error>> {
    let r1_path = r1;
    let r2_path = r2;
    let out_path = out;

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .ok();

    let reader1 = open_reader(&r1_path)?;
    let reader2 = open_reader(&r2_path)?;

    let mut r1r = Reader::new(reader1);
    let mut r2r = Reader::new(reader2);

    let (tx_batch, rx_batch) = chan::bounded::<Batch>(threads * 2);
    let (tx_tab, rx_tab) = chan::bounded::<DedupTable>(threads * 2);

    // Workers
    for _ in 0..threads {
        let rx_batch = rx_batch.clone();
        let tx_tab = tx_tab.clone();
        std::thread::spawn(move || {
            let mut local = DedupTable::with_capacity(2_000_000, 512 * 1024 * 1024);
            while let Ok(batch) = rx_batch.recv() {
                for (s, e) in batch.spans {
                    let key = &batch.data[s..e];
                    local.add(key);
                }
            }
            let _ = tx_tab.send(local);
        });
    }
    drop(tx_tab);

    // Producer
    let producer = std::thread::spawn({
        let tx_batch = tx_batch.clone();
        move || -> Result<()> {
            let mut it1 = r1r.records();
            let mut it2 = r2r.records();
            let mut batch = Batch {
                data: Vec::with_capacity(64 * 1024 * 1024),
                spans: Vec::with_capacity(batch_pairs * 2),
            };
            let mut pairs_in_batch = 0usize;

            loop {
                match (it1.next(), it2.next()) {
                    (None, None) => break,
                    (Some(_), None) | (None, Some(_)) => {
                        bail!("R1 and R2 have different numbers of records (desynced input).");
                    }
                    (Some(ra), Some(rb)) => {
                        let ra = ra.context("Error reading R1 record")?;
                        let rb = rb.context("Error reading R2 record")?;
                        let (s1, e1) = append_canonical_into(&mut batch, ra.seq());
                        batch.spans.push((s1, e1));
                        let (s2, e2) = append_canonical_into(&mut batch, rb.seq());
                        batch.spans.push((s2, e2));
                        pairs_in_batch += 1;

                        if pairs_in_batch >= batch_pairs {
                            tx_batch.send(batch).context("Failed to send batch")?;
                            batch = Batch {
                                data: Vec::with_capacity(64 * 1024 * 1024),
                                spans: Vec::with_capacity(batch_pairs * 2),
                            };
                            pairs_in_batch = 0;
                        }
                    }
                }
            }
            if !batch.spans.is_empty() {
                tx_batch.send(batch).context("Failed to send final batch")?;
            }
            Ok(())
        }
    });

    drop(tx_batch);
    producer.join().unwrap()?;

    let mut partials: Vec<DedupTable> = Vec::new();
    while let Ok(t) = rx_tab.recv() {
        partials.push(t);
    }

    let merged: DedupTable = partials.into_par_iter().reduce(
        || DedupTable::with_capacity(2_000_000, 512 * 1024 * 1024),
        merge_tables,
    );

    let mut out: Box<dyn Write> = if out_path.as_os_str() == "-" {
        Box::new(BufWriter::with_capacity(8 << 20, io::stdout().lock()))
    } else {
        let f = File::create(&out_path).with_context(|| format!("Failed to create {}", out_path.display()))?;
        Box::new(BufWriter::with_capacity(8 << 20, f))
    };

    if sort_by_size {
        let mut items: Vec<(u32, u32, u64)> = merged
            .iter_buckets()
            .map(|b| (b.off, b.len, b.count))
            .filter(|&(_, _, cnt)| cnt >= min_size as u64)
            .collect();

        items.sort_unstable_by(|a, b| b.2.cmp(&a.2));

        for (i, (off, len, cnt)) in items.into_iter().enumerate() {
            let seq = merged.arena_get(off, len);
            writeln!(out, ">NODE_{}|{}", i + 1, cnt)?;
            out.write_all(seq)?;
            writeln!(out)?;
        }
    } else {
        let mut idx: usize = 0;
        for b in merged.iter_buckets() {
            if b.count < min_size as u64 { continue; }
            idx += 1;
            let seq = merged.arena_get(b.off, b.len);
            writeln!(out, ">NODE_{}|{}", idx, b.count)?;
            out.write_all(seq)?;
            writeln!(out)?;
        }
    }

    out.flush()?;
    Ok(())
}