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
pub fn check_known_adapters(reads: &[&[u8]]) -> Option<(Vec<u8>, String)> {
    use crate::reads::known_adapters::KNOWN_ADAPTERS;

    const MAX_CHECK_READS: usize = 100_000;
    const MAX_CHECK_BASES: usize = MAX_CHECK_READS * 1000;
    const MAX_HIT: usize = 1000;
    const MATCH_REQ: usize = 8;
    const ALLOW_ONE_MISMATCH_FOR_EACH: usize = 16;

    let mut counts = vec![0usize; KNOWN_ADAPTERS.len()];
    let mut mismatches = vec![0usize; KNOWN_ADAPTERS.len()];

    let mut checked_reads = 0usize;
    let mut checked_bases = 0usize;
    let mut cur_max = 0usize;

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
        for (ai, &(adapter, _)) in KNOWN_ADAPTERS.iter().enumerate() {
            let alen = adapter.len();
            if alen >= rlen {
                continue;
            }
            // Once a front-runner exists, stop scoring the stragglers.
            if cur_max > 20 && counts[ai] < cur_max / 10 {
                continue;
            }
            for pos in 0..rlen.saturating_sub(MATCH_REQ) {
                let cmplen = (rlen - pos).min(alen);
                let allowed = cmplen / ALLOW_ONE_MISMATCH_FOR_EACH;
                let mut mismatch = 0;
                let mut matched = true;
                for i in 0..cmplen {
                    if adapter[i] != r[i + pos] {
                        mismatch += 1;
                        if mismatch > allowed {
                            matched = false;
                            break;
                        }
                    }
                }
                if matched {
                    counts[ai] += 1;
                    if cur_max < counts[ai] {
                        cur_max = counts[ai];
                    }
                    mismatches[ai] += mismatch;
                    break;
                }
            }
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
