use anyhow::{bail, Context, Result};
use crossbeam_channel as chan;
use hashbrown::hash_map::RawEntryMut;
use hashbrown::HashMap;
use rayon::prelude::*;
use seq_io::fastq::{Reader, Record};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::error::Error;
use flate2::read::MultiGzDecoder;

type FastMap = HashMap<Vec<u8>, u64, ahash::RandomState>;

#[derive(Debug)]
struct Batch {
    // All sequences concatenated back-to-back
    data: Vec<u8>,
    // Offsets into data for each sequence
    spans: Vec<(usize, usize)>,
}

fn open_reader(path: &PathBuf) -> Result<Box<dyn Read + Send>> {
    let f = File::open(path).with_context(|| format!("Failed to open {}", path.display()))?;
    let br = BufReader::with_capacity(8 << 20, f);

    // Handle gzipped FASTQ if the filename ends with ".gz"
    if path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        Ok(Box::new(MultiGzDecoder::new(br)))
    } else {
        Ok(Box::new(br))
    }
}


// Tiny helper so we don't need an extra crate.
mod num_cpus {
    pub fn get() -> usize {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    }
}
#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b'T',
        b'c' => b'G',
        b'g' => b'C',
        b't' => b'A',
        b'N' | b'n' => b'N',
        _ => b'N',
    }
}

/// Decide if forward <= revcomp(seq) lexicographically, without building revcomp.
/// Returns true if forward is canonical, false if revcomp is canonical.
#[inline]
fn forward_is_canonical(seq: &[u8]) -> bool {
    let n = seq.len();
    for i in 0..n {
        let f = seq[i];
        let r = complement(seq[n - 1 - i]);
        if f < r {
            return true;
        } else if f > r {
            return false;
        }
    }
    true // equal, pick forward
}

/// Append the RC-canonical form of seq into batch.data and return (start, end).
#[inline]
fn append_canonical_into(batch: &mut Batch, seq: &[u8]) -> (usize, usize) {
    let start = batch.data.len();
    if forward_is_canonical(seq) {
        batch.data.extend_from_slice(seq);
    } else {
        // Write revcomp directly, no intermediate Vec
        batch.data.reserve(seq.len());
        for &b in seq.iter().rev() {
            batch.data.push(complement(b));
        }
    }
    let end = batch.data.len();
    (start, end)
}

fn merge_maps(mut a: FastMap, b: FastMap) -> FastMap {
    for (k, v) in b {
        *a.entry(k).or_insert(0) += v;
    }
    a
}

/// Update map using borrowed lookup; allocate Vec only on first insertion.
#[inline]
fn add_count(map: &mut FastMap, key: &[u8]) {
    match map.raw_entry_mut().from_key(key) {
        RawEntryMut::Occupied(mut occ) => {
            *occ.get_mut() += 1;
        }
        RawEntryMut::Vacant(vac) => {
            vac.insert(key.to_vec(), 1);
        }
    }
}


pub fn fast_dedupe(
    r1: PathBuf,
    r2: PathBuf,
    out: PathBuf,
    threads: usize,
    batch_pairs: usize,
    sort_by_size: bool,
    min_size: usize,
) -> Result<(), Box<dyn Error>> {
    // Rename path args so we can reuse r1/r2 for readers later without shadowing.
    let r1_path = r1;
    let r2_path = r2;
    let out_path = out;

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .ok();

    let r1 = open_reader(&r1_path)?;
    let r2 = open_reader(&r2_path)?;

    let mut r1r = Reader::new(r1);
    let mut r2r = Reader::new(r2);

    let (tx_batch, rx_batch) = chan::bounded::<Batch>(threads * 2);
    let (tx_map, rx_map) = chan::bounded::<FastMap>(threads * 2);

    // Workers: build local maps with borrowed lookup, allocate only on new unique sequences.
    for _ in 0..threads {
        let rx_batch = rx_batch.clone();
        let tx_map = tx_map.clone();
        std::thread::spawn(move || {
            let mut local: FastMap = FastMap::with_hasher(ahash::RandomState::new());
            let mut seen_since_spill: usize = 0;

            while let Ok(batch) = rx_batch.recv() {
                for (s, e) in batch.spans {
                    let key = &batch.data[s..e];
                    add_count(&mut local, key);
                    seen_since_spill += 1;
                }

                // Spill occasionally to cap peak memory and improve parallel merge.
                if local.len() > 2_000_000 && seen_since_spill > 2_000_000 {
                    let spill = std::mem::take(&mut local);
                    let _ = tx_map.send(spill);
                    seen_since_spill = 0;
                }
            }

            let _ = tx_map.send(local);
        });
    }
    drop(tx_map);

    // Producer: read paired FASTQ in lockstep, but count mates independently (pooled).
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
                let a = it1.next();
                let b = it2.next();

                match (a, b) {
                    (None, None) => break,
                    (Some(_), None) | (None, Some(_)) => {
                        bail!("R1 and R2 have different numbers of records (desynced input).");
                    }
                    (Some(ra), Some(rb)) => {
                        let ra = ra.context("Error reading R1 record")?;
                        let rb = rb.context("Error reading R2 record")?;

                        // R1
                        let (s1, e1) = append_canonical_into(&mut batch, ra.seq());
                        batch.spans.push((s1, e1));

                        // R2
                        let (s2, e2) = append_canonical_into(&mut batch, rb.seq());
                        batch.spans.push((s2, e2));

                        pairs_in_batch += 1;

                        if pairs_in_batch >= batch_pairs {
                            tx_batch
                                .send(batch)
                                .context("Failed to send batch to workers")?;

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

    // Collect partial maps
    let mut partials: Vec<FastMap> = Vec::new();
    while let Ok(m) = rx_map.recv() {
        partials.push(m);
    }

    // Parallel reduce merge
    let merged: FastMap = partials.into_par_iter().reduce(
        || FastMap::with_hasher(ahash::RandomState::new()),
        merge_maps,
    );

    // Output
    let mut out: Box<dyn Write> = if out_path.as_os_str() == "-" {
        Box::new(BufWriter::with_capacity(8 << 20, io::stdout().lock()))
    } else {
        let f = File::create(&out_path)
            .with_context(|| format!("Failed to create {}", out_path.display()))?;
        Box::new(BufWriter::with_capacity(8 << 20, f))
    };

    if sort_by_size {
        let mut items: Vec<(&Vec<u8>, &u64)> = merged.iter().collect();
        items.retain(|(_, &cnt)| cnt >= min_size as u64);
        items.sort_unstable_by(|a, b| b.1.cmp(a.1));

        for (i, (seq, &cnt)) in items.into_iter().enumerate() {
            writeln!(out, ">NODE_{}|{}", i + 1, cnt)?;
            out.write_all(seq)?;
            writeln!(out)?;
        }
    } else {
        let mut idx: usize = 0;
        for (seq, cnt) in merged.iter() {
            if *cnt < min_size as u64 {
                continue;
            }
            idx += 1;
            writeln!(out, ">NODE_{}|{}", idx, cnt)?;
            out.write_all(seq)?;
            writeln!(out)?;
        }
    }

    out.flush()?;
    Ok(())
}
