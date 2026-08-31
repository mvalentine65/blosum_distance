//! Adapter detection from the reads themselves.
//!
//! Ported from fastp 1.3.6 `src/evaluator.cpp` and `src/nucleotidetree.cpp`
//! (MIT, (c) 2016 OpenGene).
//!
//! Two strategies, in fastp's order. First look for a known adapter riding on
//! enough reads to be real. Failing that, find the most over-represented
//! 10-mer that isn't low-complexity, then grow it outward by following the
//! dominant path through a prefix tree of what surrounds it — which recovers
//! an adapter nobody told us about, the case that matters for SRA downloads.

use crate::reads::known_adapters::match_known_adapter;

const KEYLEN: usize = 10;

/// Shortest overlap `check_known_adapters` will consider.
const MATCH_REQ: usize = 8;
/// One mismatch is forgiven per this many compared bases.
const ALLOW_ONE_MISMATCH_FOR_EACH: usize = 16;

// ---------------------------------------------------------------------------
// Nucleotide tree
// ---------------------------------------------------------------------------

/// A prefix tree over the bases following (or preceding) a seed.
///
/// fastp indexes children by `base & 0x07`, which maps A=1, C=3, G=7, T=4 and
/// gives it 8 slots; we keep the same 8-slot indexing so the traversal order
/// and therefore the dominant path are identical.
struct NucleotideNode {
    count: u32,
    base: u8,
    children: [Option<Box<NucleotideNode>>; 8],
}

impl NucleotideNode {
    fn new(base: u8) -> Self {
        NucleotideNode { count: 0, base, children: Default::default() }
    }
}

pub struct NucleotideTree {
    root: Box<NucleotideNode>,
}

impl Default for NucleotideTree {
    fn default() -> Self {
        Self::new()
    }
}

impl NucleotideTree {
    pub fn new() -> Self {
        NucleotideTree { root: Box::new(NucleotideNode::new(b'N')) }
    }

    pub fn add_seq(&mut self, seq: &[u8]) {
        let mut cur = &mut self.root;
        for &b in seq {
            if b == b'N' {
                break;
            }
            let idx = (b & 0x07) as usize;
            if cur.children[idx].is_none() {
                cur.children[idx] = Some(Box::new(NucleotideNode::new(b)));
            }
            let child = cur.children[idx].as_mut().unwrap();
            child.count += 1;
            cur = child;
        }
    }

    /// Walk while one child holds at least 95% of the traffic and the node has
    /// been visited at least 50 times. `reached_leaf` is false when the walk
    /// stopped because the path forked rather than because it ran dry.
    pub fn dominant_path(&self, reached_leaf: &mut bool) -> Vec<u8> {
        const RATIO_THRESHOLD: f64 = 0.95;
        const NUM_THRESHOLD: u32 = 50;

        let mut out = Vec::new();
        let mut cur = &self.root;
        loop {
            let total: u32 = cur.children.iter().flatten().map(|c| c.count).sum();
            if total < NUM_THRESHOLD {
                break;
            }
            let mut next = None;
            for child in cur.children.iter().flatten() {
                if f64::from(child.count) / f64::from(total) >= RATIO_THRESHOLD {
                    out.push(child.base);
                    next = Some(child);
                    break;
                }
            }
            match next {
                Some(child) => cur = child,
                None => {
                    *reached_leaf = false;
                    break;
                }
            }
        }
        out
    }
}

// ---------------------------------------------------------------------------
// k-mer packing
// ---------------------------------------------------------------------------

/// fastp's `seq2int`: 2 bits per base, A=0 T=1 C=2 G=3, -1 on anything else.
/// Pass the previous value to roll the window forward by one base.
pub fn seq2int(seq: &[u8], pos: usize, keylen: usize, last_val: i32) -> i32 {
    let code = |b: u8| -> i32 {
        match b {
            b'A' => 0,
            b'T' => 1,
            b'C' => 2,
            b'G' => 3,
            _ => -1,
        }
    };
    if last_val >= 0 {
        let mask = (1i32 << (keylen * 2)) - 1;
        let key = (last_val << 2) & mask;
        let c = code(seq[pos + keylen - 1]);
        if c < 0 {
            return -1;
        }
        key + c
    } else {
        let mut key = 0i32;
        for i in pos..pos + keylen {
            let c = code(seq[i]);
            if c < 0 {
                return -1;
            }
            key = (key << 2) + c;
        }
        key
    }
}

/// fastp's `int2seq`.
pub fn int2seq(mut val: u32, seqlen: usize) -> Vec<u8> {
    const BASES: [u8; 4] = [b'A', b'T', b'C', b'G'];
    let mut out = vec![b'N'; seqlen];
    for done in 0..seqlen {
        out[seqlen - done - 1] = BASES[(val & 0x03) as usize];
        val >>= 2;
    }
    out
}

// ---------------------------------------------------------------------------
// Detection
// ---------------------------------------------------------------------------

/// Look for a known adapter carried by enough of the sample, fastp's
/// `checkKnownAdapters`.

/// Eight bases packed two bits each, or None if any base is not ACGT. A block
/// that matches the adapter exactly cannot hold anything else, since the
/// adapter table is pure ACGT, so refusing to index those positions loses no
/// true match -- soft-masked and N bases included.
#[inline]
fn pack8(b: &[u8]) -> Option<u16> {
    let mut key = 0u16;
    for &c in &b[..8] {
        let two = match c {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        key = (key << 2) | two;
    }
    Some(key)
}

#[inline]
fn base2bits(c: u8) -> Option<u16> {
    match c {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Which adapter blocks carry each eight-base key.
///
/// `check_known_adapters` forgives `cmplen / ALLOW_ONE_MISMATCH_FOR_EACH`
/// mismatches over a window of `cmplen` bases, so cutting the adapter into
/// `allowed + 1` disjoint eight-base blocks leaves at least one of them
/// mismatch-free whenever the window matches. Blocks sit at offsets 0, 8,
/// 16 ... and `8 * (allowed + 1) <= cmplen` holds for every window the scan can
/// reach (its loop bound keeps `cmplen >= 9`), so the block a given window is
/// guaranteed to match on is always one this table carries. Sizing from the
/// adapter's own length makes it a superset for the shorter windows at a
/// read's tail, which is equally safe: extra blocks only add candidates, and
/// every candidate is still verified by the original window test.
///
/// Indexed by key so a read costs one lookup per position, rather than every
/// adapter walking every position of every read.
struct SeedTable {
    /// 1 << 16 buckets, each the head of a chain into `entries`.
    head: Vec<u32>,
    next: Vec<u32>,
    /// (adapter index, block offset)
    entries: Vec<(u16, u8)>,
    /// Adapters too short to carry a block; these still need a full scan.
    seedless: Vec<u16>,
}

const NO_ENTRY: u32 = u32::MAX;
/// Marks an adapter that has no seed block and must be scanned in full.
const SCAN_ALL: u32 = u32::MAX;

fn seed_table() -> &'static SeedTable {
    use crate::reads::known_adapters::KNOWN_ADAPTERS;
    static TABLE: std::sync::OnceLock<SeedTable> = std::sync::OnceLock::new();
    TABLE.get_or_init(|| {
        let mut t = SeedTable {
            head: vec![NO_ENTRY; 1 << 16],
            next: Vec::new(),
            entries: Vec::new(),
            seedless: Vec::new(),
        };
        for (ai, &(adapter, _)) in KNOWN_ADAPTERS.iter().enumerate() {
            let alen = adapter.len();
            let mut planted = false;
            if alen >= 8 {
                let allowed_max = alen / ALLOW_ONE_MISMATCH_FOR_EACH;
                for j in 0..=allowed_max {
                    let off = j * 8;
                    if off + 8 > alen {
                        break;
                    }
                    let Some(key) = pack8(&adapter[off..off + 8]) else {
                        continue;
                    };
                    let k = key as usize;
                    t.entries.push((ai as u16, off as u8));
                    t.next.push(t.head[k]);
                    t.head[k] = (t.entries.len() - 1) as u32;
                    planted = true;
                }
            }
            if !planted {
                t.seedless.push(ai as u16);
            }
        }
        t
    })
}

/// fastp's window test: the adapter's prefix against the read from `pos`,
/// forgiving one mismatch per ALLOW_ONE_MISMATCH_FOR_EACH bases compared.
/// Returns the mismatch count when it passes.
#[inline]
fn window_match(adapter: &[u8], r: &[u8], pos: usize) -> Option<usize> {
    let cmplen = (r.len() - pos).min(adapter.len());
    let allowed = cmplen / ALLOW_ONE_MISMATCH_FOR_EACH;
    let mut mismatch = 0;
    for i in 0..cmplen {
        if adapter[i] != r[i + pos] {
            mismatch += 1;
            if mismatch > allowed {
                return None;
            }
        }
    }
    Some(mismatch)
}

pub fn check_known_adapters(reads: &[&[u8]]) -> Option<(Vec<u8>, String)> {
    use crate::reads::known_adapters::KNOWN_ADAPTERS;

    const MAX_CHECK_READS: usize = 100_000;
    const MAX_CHECK_BASES: usize = MAX_CHECK_READS * 1000;
    const MAX_HIT: usize = 1000;

    let mut counts = vec![0usize; KNOWN_ADAPTERS.len()];
    let mut mismatches = vec![0usize; KNOWN_ADAPTERS.len()];

    let mut checked_reads = 0usize;
    let mut checked_bases = 0usize;
    let mut cur_max = 0usize;

    let table = seed_table();
    // (adapter index, candidate offset), collected per read.
    let mut candidates: Vec<(u16, u32)> = Vec::new();

    for r in reads {
        let rlen = r.len();
        checked_reads += 1;
        checked_bases += rlen;
        if checked_reads > MAX_CHECK_READS || checked_bases > MAX_CHECK_BASES {
            break;
        }
        if cur_max > MAX_HIT {
            break;
        }
        let limit = rlen.saturating_sub(MATCH_REQ);

        // One rolling pass over the read proposes every adapter's candidate
        // offsets. A true match always shows up here, by the pigeonhole
        // argument on SeedTable; each candidate is still verified below.
        candidates.clear();
        if rlen >= 8 {
            let mut key = 0u16;
            let mut run = 0usize;
            for (pos, &c) in r.iter().enumerate() {
                match base2bits(c) {
                    Some(two) => {
                        key = (key << 2) | two;
                        run += 1;
                    }
                    None => {
                        run = 0;
                        continue;
                    }
                }
                if run < 8 {
                    continue;
                }
                let start = pos + 1 - 8;
                let mut e = table.head[key as usize];
                while e != NO_ENTRY {
                    let (ai, off) = table.entries[e as usize];
                    let off = off as usize;
                    if start >= off && start - off < limit {
                        candidates.push((ai, (start - off) as u32));
                    }
                    e = table.next[e as usize];
                }
            }
        }
        // Adapters with no block to seed on still need a full scan. Marking
        // them here rather than looping separately keeps adapters visited in
        // index order, which matters because cur_max gates later ones.
        for &ai in &table.seedless {
            candidates.push((ai, SCAN_ALL));
        }
        // The old scan took the lowest matching offset for each adapter, so
        // order by adapter then offset before testing.
        candidates.sort_unstable();
        candidates.dedup();

        let mut i = 0usize;
        while i < candidates.len() {
            let ai = candidates[i].0;
            let mut j = i;
            while j < candidates.len() && candidates[j].0 == ai {
                j += 1;
            }
            let (adapter, _) = KNOWN_ADAPTERS[ai as usize];
            if adapter.len() < rlen && !(cur_max > 20 && counts[ai as usize] < cur_max / 10) {
                let scan_all = candidates[j - 1].1 == SCAN_ALL;
                let mut hit = None;
                if scan_all {
                    for pos in 0..limit {
                        if let Some(mismatch) = window_match(adapter, r, pos) {
                            hit = Some(mismatch);
                            break;
                        }
                    }
                } else {
                    for &(_, pos) in &candidates[i..j] {
                        if let Some(mismatch) = window_match(adapter, r, pos as usize) {
                            hit = Some(mismatch);
                            break;
                        }
                    }
                }
                if let Some(mismatch) = hit {
                    counts[ai as usize] += 1;
                    if cur_max < counts[ai as usize] {
                        cur_max = counts[ai as usize];
                    }
                    mismatches[ai as usize] += mismatch;
                }
            }
            i = j;
        }

    }

    let mut best = None;
    let mut max_count = 0usize;
    for (ai, &c) in counts.iter().enumerate() {
        if c > max_count {
            max_count = c;
            best = Some(ai);
        }
    }
    let ai = best?;
    // Either carried by 2% of reads, or by 0.5% and cleanly matched.
    if max_count > checked_reads / 50
        || (max_count > checked_reads / 200 && mismatches[ai] < checked_reads)
    {
        let (seq, name) = KNOWN_ADAPTERS[ai];
        return Some((seq.to_vec(), name.to_string()));
    }
    None
}

/// Grow a seed k-mer into a full adapter, fastp's `getAdapterWithSeed`.
fn adapter_with_seed(seed: i32, reads: &[&[u8]], shift_tail: usize) -> Option<Detected> {
    const MAX_SEARCH_LENGTH: usize = 500;

    let mut forward = NucleotideTree::new();
    let mut backward = NucleotideTree::new();

    for r in reads {
        let rlen = r.len();
        if rlen < KEYLEN + shift_tail + 20 {
            continue;
        }
        let last = rlen - KEYLEN - shift_tail;
        let mut key = -1i32;
        for pos in 20..=last.min(MAX_SEARCH_LENGTH) {
            key = seq2int(r, pos, KEYLEN, key);
            if key == seed {
                forward.add_seq(&r[pos + KEYLEN..rlen - shift_tail]);
                let mut rev: Vec<u8> = r[..pos].to_vec();
                rev.reverse();
                backward.add_seq(&rev);
            }
        }
    }

    let mut reached_leaf = true;
    let forward_path = forward.dominant_path(&mut reached_leaf);
    let backward_path = backward.dominant_path(&mut reached_leaf);

    let mut adapter: Vec<u8> = backward_path.into_iter().rev().collect();
    adapter.extend_from_slice(&int2seq(seed as u32, KEYLEN));
    adapter.extend_from_slice(&forward_path);
    if adapter.len() > 60 {
        adapter.truncate(60);
    }

    if let Some((known, name)) = match_known_adapter(&adapter) {
        return Some(Detected::Known { seq: known.to_vec(), name: name.to_string() });
    }
    if reached_leaf {
        Some(Detected::Novel { seq: adapter })
    } else {
        None
    }
}

/// What detection concluded.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Detected {
    /// Matched an entry in the known table.
    Known { seq: Vec<u8>, name: String },
    /// Assembled from the reads, matching nothing known.
    Novel { seq: Vec<u8> },
    /// Nothing convincing — the library looks clean, or the sample is too small.
    None,
}

/// Detect the adapter carried by a sample of reads, fastp's
/// `evalAdapterAndReadNum` minus its read-count estimation.
///
/// `shift_tail` drops the last cycle from consideration, which fastp does
/// because it is reliably the noisiest.
pub fn detect_adapter(reads: &[&[u8]], shift_tail: usize) -> Detected {
    // fastp refuses to guess from a small sample, and so do we: below this the
    // over-representation test has no power and would invent an adapter.
    if reads.len() < 10_000 {
        return Detected::None;
    }

    if let Some((seq, name)) = check_known_adapters(reads) {
        if seq.len() > 8 {
            return Detected::Known { seq, name };
        }
    }

    let shift_tail = shift_tail.max(1);
    let size = 1usize << (KEYLEN * 2);
    let mut counts = vec![0u32; size];
    for r in reads {
        let rlen = r.len();
        if rlen < KEYLEN + shift_tail + 20 {
            continue;
        }
        let last = rlen - KEYLEN - shift_tail;
        let mut key = -1i32;
        for pos in 20..=last {
            key = seq2int(r, pos, KEYLEN, key);
            if key >= 0 {
                counts[key as usize] += 1;
            }
        }
    }
    counts[0] = 0; // AAAAAAAAAA

    // Rank candidate keys, skipping the ones that are structurally junk.
    const TOPNUM: usize = 10;
    let mut topkeys = [0usize; TOPNUM];
    let mut total: u64 = 0;
    for k in 0..size {
        let mut atcg = [0usize; 4];
        for i in 0..KEYLEN {
            atcg[(k >> (i * 2)) & 0x03] += 1;
        }
        if atcg.iter().any(|&n| n >= KEYLEN - 4) {
            continue; // low complexity
        }
        if atcg[2] + atcg[3] >= KEYLEN - 2 {
            continue; // GC heavy
        }
        if k >> 12 == 0xff {
            continue; // starts with GGGG
        }
        let val = counts[k];
        total += u64::from(val);
        for t in (0..TOPNUM).rev() {
            if val < counts[topkeys[t]] {
                if t < TOPNUM - 1 {
                    for m in (t + 2..TOPNUM).rev() {
                        topkeys[m] = topkeys[m - 1];
                    }
                    topkeys[t + 1] = k;
                }
                break;
            } else if t == 0 {
                for m in (1..TOPNUM).rev() {
                    topkeys[m] = topkeys[m - 1];
                }
                topkeys[0] = k;
            }
        }
    }

    const FOLD_THRESHOLD: u64 = 20;
    for &key in topkeys.iter() {
        if key == 0 {
            continue;
        }
        let count = u64::from(counts[key]);
        if count < 10 || count * (size as u64) < total * FOLD_THRESHOLD {
            break;
        }
        let seq = int2seq(key as u32, KEYLEN);
        let diff = seq.windows(2).filter(|w| w[0] != w[1]).count();
        if diff < 3 {
            continue;
        }
        if let Some(detected) = adapter_with_seed(key as i32, reads, shift_tail) {
            return detected;
        }
    }

    Detected::None
}

#[cfg(test)]
#[path = "../tests/detect.rs"]
mod tests;
