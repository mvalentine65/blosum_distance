use anyhow::{bail, Context, Result};
use seq_io::fastq::{Reader, Record};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::error::Error;
use flate2::read::MultiGzDecoder;

// --- Optimized DedupTable from main_bestrs.rs ---

#[derive(Clone, Copy, Default)]
struct Bucket {
    hash: u64,
    off: usize, 
    len: usize,
    count: u64,
    used: bool,
}

struct DedupTable {
    buckets: Vec<Bucket>,
    mask: usize,
    arena: Vec<u8>,
    used: usize,
}

impl DedupTable {
    fn with_capacity(expected_uniques: usize) -> Self {
        let mut cap = (expected_uniques * 10) / 7;
        cap = cap.next_power_of_two().max(1024);
        Self {
            buckets: vec![Bucket::default(); cap],
            mask: cap - 1,
            arena: Vec::with_capacity(512 * 1024 * 1024),
            used: 0,
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
        for b in old.into_iter().filter(|b| b.used) {
            let mut i = (b.hash as usize) & self.mask;
            while self.buckets[i].used {
                i = (i + 1) & self.mask;
            }
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
        let len = seq.len();
        let mut i = (h as usize) & self.mask;
        loop {
            if !self.buckets[i].used {
                let off = self.arena.len();
                self.arena.extend_from_slice(seq);
                self.buckets[i] = Bucket { hash: h, off, len, count: 1, used: true };
                self.used += 1;
                return;
            }
            if self.buckets[i].hash == h && self.buckets[i].len == len {
                if self.arena_get(self.buckets[i].off, len) == seq {
                    self.buckets[i].count += 1;
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

#[inline]
fn get_canonical(seq: &[u8]) -> Vec<u8> {
    if forward_is_canonical(seq) {
        seq.to_vec()
    } else {
        seq.iter().rev().map(|&b| complement(b)).collect()
    }
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


pub fn fast_dedupe(
    input_paths: Vec<PathBuf>, // Changed from r1, r2 to a Vec
    out: PathBuf,
    sort_by_size: bool,
    min_size: u64,
) -> Result<(), Box<dyn Error>> {
    // Single-threaded table initialization
    let mut table = DedupTable::with_capacity(2_000_000);
    
    // Reuse a buffer for canonicalization to avoid massive allocation overhead
    let mut scratch = Vec::with_capacity(300);

    // PHASE 1: DEDUPLICATION
    for path in input_paths {
        let mut reader = Reader::new(open_reader(&path)?);
        let mut records = reader.records();

        while let Some(record) = records.next() {
            let rec = record.context("Error reading FASTQ record")?;
            let seq = rec.seq();
            
            // In-place canonicalization logic to save memory/time
            if forward_is_canonical(seq) {
                table.add(seq);
            } else {
                scratch.clear();
                scratch.extend(seq.iter().rev().map(|&b| complement(b)));
                table.add(&scratch);
            }
        }
    }

    // PHASE 2: WRITING
    let mut writer: Box<dyn Write> = if out.as_os_str() == "-" {
        Box::new(BufWriter::with_capacity(8 << 20, io::stdout().lock()))
    } else {
        let f = File::create(&out).with_context(|| format!("Failed to create {}", out.display()))?;
        Box::new(BufWriter::with_capacity(8 << 20, f))
    };

    let mut items: Vec<Bucket> = table.buckets.iter().copied()
        .filter(|b| b.used && b.count >= min_size)
        .collect();

    if sort_by_size {
        items.sort_unstable_by(|a, b| b.count.cmp(&a.count));
    }

    for (i, b) in items.into_iter().enumerate() {
        writeln!(writer, ">{}|{}", i + 1, b.count)?;
        writer.write_all(table.arena_get(b.off, b.len))?;
        writeln!(writer)?;
    }
    
    writer.flush()?;
    Ok(())
}