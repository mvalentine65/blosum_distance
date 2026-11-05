//! The read loops of SAPPHYRE's flanks mode: placing a gene's reads on its merged
//! row, building consensus outward from an edge, and checking reads recruited
//! to an edge against what the edge holds. Each function gives exactly what the
//! Python it replaces gives, ties included.
use ahash::AHashMap;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

/// Complement as Python's `str.maketrans("ACGTN", "TGCAN")`: other bytes kept.
#[inline]
fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        x => x,
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| comp(b)).collect()
}

/// Upper-case A, C, G, T only, as the Python membership tests in "ACGT".
#[inline]
fn base(b: u8) -> u64 {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

const BASES: [u8; 4] = *b"ACGT";

/// Packed k-mer of every all-ACGT window, None for the others.
fn windows(seq: &[u8], k: usize) -> Vec<Option<u64>> {
    if seq.len() < k {
        return Vec::new();
    }
    let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let mut out = Vec::with_capacity(seq.len() - k + 1);
    let (mut kmer, mut run) = (0u64, 0usize);
    for (i, &b) in seq.iter().enumerate() {
        let c = base(b);
        if c == 4 {
            run = 0;
        } else {
            kmer = ((kmer << 2) | c) & mask;
            run += 1;
        }
        if i + 1 >= k {
            out.push(if run >= k { Some(kmer) } else { None });
        }
    }
    out
}

/// Start of every k-mer seen once in the row.
fn row_index(seq: &[u8], k: usize) -> AHashMap<u64, i64> {
    let mut seen: AHashMap<u64, (u32, i64)> = AHashMap::new();
    for (i, w) in windows(seq, k).into_iter().enumerate() {
        if let Some(km) = w {
            let e = seen.entry(km).or_insert((0, i as i64));
            e.0 += 1;
        }
    }
    seen.into_iter().filter(|(_, (n, _))| *n == 1).map(|(km, (_, p))| (km, p)).collect()
}

fn place(read: &[u8], index: &AHashMap<u64, i64>, k: usize) -> Vec<i64> {
    let mut pos = vec![-1i64; read.len()];
    for (i, w) in windows(read, k).into_iter().enumerate() {
        if let Some(p) = w.and_then(|km| index.get(&km)) {
            for j in 0..k {
                pos[i + j] = p + j as i64;
            }
        }
    }
    pos
}

fn lookup(seq: &[u8], index: &AHashMap<u64, i64>, k: usize) -> Option<i64> {
    windows(seq, k).first().copied().flatten().and_then(|km| index.get(&km).copied())
}

/// Both end k-mers on the row, the read's length apart.
fn inside(read: &[u8], index: &AHashMap<u64, i64>, k: usize) -> Option<i64> {
    if read.len() < k {
        return None;
    }
    let a = lookup(&read[..k], index, k)?;
    let b = lookup(&read[read.len() - k..], index, k)?;
    if b - a == (read.len() - k) as i64 {
        Some(a)
    } else {
        None
    }
}

fn mismatch_tail(piece: &[u8], continuation: &[u8], min_piece: usize, tail_match: f64) -> bool {
    let n = piece.len().min(continuation.len());
    if n < min_piece {
        return false;
    }
    let same = piece.iter().zip(continuation).filter(|(a, b)| a == b).count();
    same as f64 / n as f64 >= tail_match
}

struct Ordered<K> {
    keys: Vec<K>,
    at: AHashMap<K, usize>,
    vals: Vec<Vec<Vec<u8>>>,
}

impl<K: std::hash::Hash + Eq + Clone> Ordered<K> {
    fn new() -> Self {
        Ordered { keys: Vec::new(), at: AHashMap::new(), vals: Vec::new() }
    }
    fn push(&mut self, key: K, piece: Vec<u8>) {
        let i = *self.at.entry(key.clone()).or_insert_with(|| {
            self.keys.push(key);
            self.vals.push(Vec::new());
            self.vals.len() - 1
        });
        self.vals[i].push(piece);
    }
}

type Edges = Vec<((String, i64), Vec<String>)>;
type Whole = Vec<(i64, Vec<String>)>;

/// flanks.fit_junction: the intron that, written after row base `p`, gives back the read. The read holds
/// seq[..=p] + piece + seq[q..]; the row is written as seq[..=p] + intron + seq[p + 1..]. Junction bases both
/// exons share move `q`: with q <= p the read holds seq[q..=p] again after the piece; with q > p + 1 the row holds
/// bases the read lacks, removed from the piece's end when it ends with them. None when it does not.
fn fit_junction(piece: &[u8], seq: &[u8], p: i64, q: i64) -> Option<Vec<u8>> {
    if q == p + 1 {
        return Some(piece.to_vec());
    }
    if q <= p {
        let mut out = piece.to_vec();
        out.extend_from_slice(&seq[q as usize..=p as usize]);
        return Some(out);
    }
    let extra = &seq[p as usize + 1..q as usize];
    if piece.len() > extra.len() && piece.ends_with(extra) {
        return Some(piece[..piece.len() - extra.len()].to_vec());
    }
    None
}

fn strings(v: Vec<Vec<u8>>) -> Vec<String> {
    v.into_iter().map(|p| String::from_utf8_lossy(&p).into_owned()).collect()
}

/// flanks.first_pass: overhangs per growing edge, held outward, and whole introns,
/// from the reads behind one row; exits inside the row kept only where their reads
/// are `exit_share` or more of the reads crossing that boundary.
#[pyfunction]
#[pyo3(signature = (seq, reads, k, min_piece, tail_match, exit_share))]
pub fn flank_first_pass(
    py: Python<'_>,
    seq: &str,
    reads: Vec<String>,
    k: usize,
    min_piece: usize,
    tail_match: f64,
    exit_share: f64,
) -> (Edges, Whole) {
    py.detach(|| {
        let seq = seq.as_bytes();
        let index = row_index(seq, k);
        let mut edges: Ordered<(u8, i64)> = Ordered::new();
        let mut whole: Ordered<i64> = Ordered::new();
        let mut across = vec![0i64; seq.len() + 1];
        let n_seq = seq.len() as i64;
        for read in &reads {
            let fwd = read.as_bytes();
            let bwd = revcomp(fwd);
            if let Some(a) = inside(fwd, &index, k).or_else(|| inside(&bwd, &index, k)) {
                across[a as usize] += 1;
                across[a as usize + fwd.len() - 1] -= 1;
                continue;
            }
            let fw = place(fwd, &index, k);
            let bw = place(&bwd, &index, k);
            let miss = |p: &Vec<i64>| p.iter().filter(|&&x| x == -1).count();
            let (r, pos): (&[u8], &Vec<i64>) = if miss(&bw) < miss(&fw) { (&bwd[..], &bw) } else { (fwd, &fw) };
            if miss(pos) == pos.len() {
                continue;
            }
            for i in 0..pos.len().saturating_sub(1) {
                if pos[i] >= 0 && pos[i + 1] == pos[i] + 1 {
                    across[pos[i] as usize] += 1;
                    across[pos[i] as usize + 1] -= 1;
                }
            }
            // uncovered runs of the read, as flanks._exits
            let mut i = 0usize;
            while i < pos.len() {
                if pos[i] >= 0 {
                    i += 1;
                    continue;
                }
                let mut j = i;
                while j < pos.len() && pos[j] < 0 {
                    j += 1;
                }
                let (a, b) = (i, j);
                i = j;
                if b - a < min_piece {
                    continue;
                }
                let piece = &r[a..b];
                if a == 0 && b < r.len() {
                    let first = pos[b];
                    if first <= 3 {
                        edges.push((b'l', 0), revcomp(piece));
                    } else {
                        let p_rev: Vec<u8> = piece.iter().rev().copied().collect();
                        let row_rev: Vec<u8> = seq[..first as usize].iter().rev().copied().collect();
                        if !mismatch_tail(&p_rev, &row_rev, min_piece, tail_match) {
                            edges.push((b'i', first), revcomp(piece));
                        }
                    }
                } else if b == r.len() && a > 0 {
                    let last = pos[a - 1];
                    if last >= n_seq - 4 {
                        edges.push((b'r', 0), piece.to_vec());
                    } else if !mismatch_tail(piece, &seq[last as usize + 1..], min_piece, tail_match) {
                        edges.push((b'o', last), piece.to_vec());
                    }
                } else if 0 < a && b < r.len() && (0..=3).contains(&(pos[b] - pos[a - 1])) {
                    if let Some(fit) = fit_junction(piece, seq, pos[a - 1], pos[b]) {
                        if fit.len() >= min_piece {
                            whole.push(pos[a - 1], fit);
                        }
                    }
                }
            }
        }
        for t in 1..across.len() {
            across[t] += across[t - 1];
        }
        let enough = |leaving: usize, p: i64| leaving as f64 >= exit_share * (leaving as i64 + across[p as usize]) as f64;
        let mut out_edges = Edges::new();
        for (key, vals) in edges.keys.into_iter().zip(edges.vals) {
            let (side, p) = key;
            let keep = match side {
                b'l' | b'r' => true,
                b'o' => enough(vals.len(), p),
                _ => enough(vals.len(), p - 1),
            };
            if keep {
                let name = match side {
                    b'l' => "left",
                    b'r' => "right",
                    b'o' => "out",
                    _ => "in",
                };
                out_edges.push(((name.to_string(), p), strings(vals)));
            }
        }
        let mut out_whole = Whole::new();
        for (p, vals) in whole.keys.into_iter().zip(whole.vals) {
            if enough(vals.len(), p) {
                out_whole.push((p, strings(vals)));
            }
        }
        (out_edges, out_whole)
    })
}

fn extend_core(overhangs: &[&[u8]], follow: bool, depth: usize, agree: f64) -> Vec<u8> {
    if overhangs.len() < depth {
        return Vec::new();
    }
    let width = overhangs.iter().map(|o| o.len()).max().unwrap_or(0);
    let mut live: Vec<&[u8]> = overhangs.to_vec();
    let mut out = Vec::new();
    for c in 0..width {
        let mut counts = [0usize; 4];
        let mut first = [usize::MAX; 4];
        for (ri, o) in live.iter().enumerate() {
            if o.len() > c {
                let b = base(o[c]) as usize;
                if b < 4 {
                    counts[b] += 1;
                    if first[b] == usize::MAX {
                        first[b] = ri;
                    }
                }
            }
        }
        let n: usize = counts.iter().sum();
        if n < depth {
            break;
        }
        let votes = *counts.iter().max().unwrap();
        let pick = (0..4).filter(|&b| counts[b] == votes).min_by_key(|&b| first[b]).unwrap();
        if (votes as f64) / (n as f64) < agree {
            if !follow || votes < depth {
                break;
            }
            live.retain(|o| o.len() > c && o[c] == BASES[pick]);
        }
        out.push(BASES[pick]);
    }
    out
}

/// flanks.extend: majority consensus of overhangs read outward from an edge,
/// stopping where reads split under `agree`, or with `follow` keeping only the
/// reads that carry the majority base. Ties go to the base seen first.
#[pyfunction]
#[pyo3(signature = (overhangs, follow, depth, agree))]
pub fn flank_extend(overhangs: Vec<String>, follow: bool, depth: usize, agree: f64) -> String {
    let views: Vec<&[u8]> = overhangs.iter().map(|o| o.as_bytes()).collect();
    String::from_utf8(extend_core(&views, follow, depth, agree)).unwrap()
}

/// Bases past the end of `hold` that read r offers at `start`, when the overlap
/// matches at `min_agree` or better over positions both hold ACGT.
fn overlap_tail<'a>(s: &[u8], r: &'a [u8], start: i64, min_agree: f64) -> Option<&'a [u8]> {
    let cut = s.len() as i64 - start;
    let cut = if cut < 0 { (r.len() as i64 + cut).max(0) } else { cut } as usize;
    if cut >= r.len() {
        return None;
    }
    let (mut same, mut both) = (0usize, 0usize);
    let s_from = start.max(0) as usize;
    let r_from = (-start).max(0) as usize;
    if s_from < s.len() && r_from < r.len() {
        for (&x, &y) in s[s_from..].iter().zip(&r[r_from..]) {
            if base(x) < 4 && base(y) < 4 {
                both += 1;
                same += (x == y) as usize;
            }
        }
    }
    if both > 0 && same as f64 / both as f64 >= min_agree {
        Some(&r[cut..])
    } else {
        None
    }
}

/// The growth past an edge each recruited read offers: for read r placed at
/// `start` against the edge's held sequence, its bases past the end, kept when
/// the overlap matches at `min_agree` or better over positions both hold ACGT.
#[pyfunction]
#[pyo3(signature = (hold, reads, starts, min_agree))]
pub fn flank_overlap_tails(hold: &str, reads: Vec<String>, starts: Vec<i64>, min_agree: f64) -> Vec<String> {
    reads
        .iter()
        .zip(&starts)
        .filter_map(|(r, &st)| overlap_tail(hold.as_bytes(), r.as_bytes(), st, min_agree))
        .map(|t| String::from_utf8_lossy(t).into_owned())
        .collect()
}

/// Seeds found in reads, for flanks.grow: each blob record holding a seed on
/// either strand comes back as (strand, seed index, offset), the strand read
/// forward when it holds a seed, else its reverse complement, and the offset of
/// the first seed on that strand, as the Python search over both strands gives.
#[pyclass]
pub struct SeedHits {
    // forward seeds, for the strand and offset of a hit read
    seeds: AHashMap<u64, usize>,
    // the smaller of each seed and its reverse complement, for the one-lookup scan of every read
    canon: ahash::AHashSet<u64>,
    filter: Vec<u64>,
    k: usize,
}

// 2^22 bits, 512 KB: small enough to stay in L2 while the scan checks one k-mer per base
const FILTER_BITS: u32 = 22;

// case-insensitive A, C, G, T; 4 for anything else
const ANY_CASE: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
};

fn rc_kmer(kmer: u64, k: usize) -> u64 {
    let mut out = 0u64;
    let mut x = kmer;
    for _ in 0..k {
        out = (out << 2) | (3 - (x & 3));
        x >>= 2;
    }
    out
}

#[inline]
fn slot(kmer: u64) -> usize {
    (kmer.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> (64 - FILTER_BITS)) as usize
}

impl SeedHits {
    fn build(seeds: &[Vec<u8>], k: usize) -> Self {
        let mut map = AHashMap::with_capacity(seeds.len());
        let mut canon = ahash::AHashSet::with_capacity(seeds.len());
        let mut filter = vec![0u64; 1 << (FILTER_BITS - 6)];
        for (idx, seed) in seeds.iter().enumerate() {
            if let Some(km) = windows(seed, k).first().copied().flatten() {
                map.entry(km).or_insert(idx);
                let c = km.min(rc_kmer(km, k));
                canon.insert(c);
                let s = slot(c);
                filter[s >> 6] |= 1 << (s & 63);
            }
        }
        SeedHits { seeds: map, canon, filter, k }
    }

    /// (strand, seed index, offset) for every record of one blob holding a seed.
    fn hits(&self, batch: &[u8]) -> Result<Vec<(Vec<u8>, usize, usize)>, usize> {
        let mut found = Vec::new();
        let mut pos = 0usize;
        while pos < batch.len() {
            let header_end = match memchr::memchr(b'\n', &batch[pos..]) {
                Some(i) => pos + i,
                None => break,
            };
            if batch[pos] != b'>' {
                return Err(pos);
            }
            pos = header_end + 1;
            let seq_end = memchr::memchr(b'\n', &batch[pos..]).map(|i| pos + i).unwrap_or(batch.len());
            if self.any_seed(&batch[pos..seq_end]) {
                let up = batch[pos..seq_end].to_ascii_uppercase();
                if let Some((idx, off)) = self.first_seed(&up) {
                    found.push((up, idx, off));
                } else {
                    let rc = revcomp(&up);
                    if let Some((idx, off)) = self.first_seed(&rc) {
                        found.push((rc, idx, off));
                    }
                }
            }
            pos = seq_end + 1;
        }
        Ok(found)
    }

    #[inline]
    fn known(&self, canonical: u64) -> bool {
        let s = slot(canonical);
        self.filter[s >> 6] >> (s & 63) & 1 == 1 && self.canon.contains(&canonical)
    }

    /// Any seed on either strand, in one pass and without copying the read.
    fn any_seed(&self, seq: &[u8]) -> bool {
        let mask = if self.k == 32 { u64::MAX } else { (1u64 << (2 * self.k)) - 1 };
        let top = 2 * (self.k as u32 - 1);
        let (mut fwd, mut rev, mut run) = (0u64, 0u64, 0usize);
        for &b in seq {
            let c = ANY_CASE[b as usize] as u64;
            if c == 4 {
                run = 0;
                continue;
            }
            fwd = ((fwd << 2) | c) & mask;
            rev = (rev >> 2) | ((3 - c) << top);
            run += 1;
            if run >= self.k && self.known(fwd.min(rev)) {
                return true;
            }
        }
        false
    }

    fn first_seed(&self, seq: &[u8]) -> Option<(usize, usize)> {
        windows(seq, self.k)
            .into_iter()
            .enumerate()
            .find_map(|(i, w)| w.and_then(|km| self.seeds.get(&km).map(|&idx| (idx, i))))
    }
}

#[pymethods]
impl SeedHits {
    #[new]
    fn new(seeds: Vec<String>, k: usize) -> PyResult<Self> {
        if k == 0 || k > 32 {
            return Err(pyo3::exceptions::PyValueError::new_err("k must be between 1 and 32"));
        }
        let bytes: Vec<Vec<u8>> = seeds.iter().map(|q| q.as_bytes().to_vec()).collect();
        for (seed, q) in seeds.iter().zip(&bytes) {
            if q.len() != k || windows(q, k).first().copied().flatten().is_none() {
                return Err(pyo3::exceptions::PyValueError::new_err(format!("seed {seed:?} is not {k} bases of ACGT")));
            }
        }
        Ok(SeedHits::build(&bytes, k))
    }

    fn scan(&self, py: Python<'_>, batch: &[u8]) -> PyResult<Vec<(Py<PyBytes>, usize, usize)>> {
        let found = py.detach(|| self.hits(batch));
        let found = found.map_err(|pos| {
            pyo3::exceptions::PyValueError::new_err(format!("ntbatch record at byte {pos} does not start with '>'"))
        })?;
        Ok(found.into_iter().map(|(s, i, o)| (PyBytes::new(py, &s).unbind(), i, o)).collect())
    }
}

/// flanks.grow: extend every edge of a sample over all its reads, `rounds` read
/// passes for all edges at once. Edge e grows from `holds[e]` (its row side held
/// outward plus any growth so far) on row `edge_rows[e]`; `batch(i)` returns the
/// sample's i-th ntbatch blob. Returns each edge's growth past `bases[e]`.
#[pyfunction]
#[pyo3(signature = (rows, edge_rows, holds, bases, follows, rounds, batch, n_batches, seed, k, max_hits, read_agree, depth, agree, max_grow))]
#[allow(clippy::too_many_arguments)]
pub fn flank_grow(
    py: Python<'_>,
    rows: Vec<String>,
    edge_rows: Vec<usize>,
    holds: Vec<String>,
    bases: Vec<usize>,
    follows: Vec<bool>,
    rounds: usize,
    batch: &Bound<'_, PyAny>,
    n_batches: usize,
    seed: usize,
    k: usize,
    max_hits: usize,
    read_agree: f64,
    depth: usize,
    agree: f64,
    max_grow: usize,
) -> PyResult<Vec<String>> {
    // k-mers seen once in each row or its reverse complement: growth that reaches one has run back into the row
    let row_kmers: Vec<ahash::AHashSet<u64>> = rows
        .iter()
        .map(|r| {
            let r = r.as_bytes();
            let mut set: ahash::AHashSet<u64> = row_index(r, k).into_keys().collect();
            set.extend(row_index(&revcomp(r), k).into_keys());
            set
        })
        .collect();
    let mut hold: Vec<Vec<u8>> = holds.into_iter().map(|h| h.into_bytes()).collect();
    let mut live: Vec<bool> = vec![true; hold.len()];
    for _ in 0..rounds {
        // seeds from each live edge's outer bases; a seed two edges share, or one edge twice, is dropped
        let mut owner: AHashMap<Vec<u8>, (usize, usize, usize)> = AHashMap::new(); // seed -> (count, edge, position)
        let mut order: Vec<Vec<u8>> = Vec::new();
        for (e, h) in hold.iter().enumerate() {
            if !live[e] {
                continue;
            }
            for back in [seed, seed + 10] {
                if h.len() < back {
                    continue;
                }
                let q = &h[h.len() - back..h.len() - back + seed];
                if q.iter().all(|&b| base(b) < 4) {
                    let entry = owner.entry(q.to_vec()).or_insert_with(|| {
                        order.push(q.to_vec());
                        (0, e, 0)
                    });
                    *entry = (entry.0 + 1, e, h.len() - back);
                }
            }
        }
        let kept: Vec<Vec<u8>> = order.into_iter().filter(|q| owner[q].0 == 1).collect();
        if kept.is_empty() {
            break;
        }
        let finder = SeedHits::build(&kept, seed);
        // per edge, recruited strands in the order first seen, each with its latest start
        let mut hits: Vec<Vec<(Vec<u8>, i64)>> = vec![Vec::new(); hold.len()];
        let mut seen: Vec<AHashMap<Vec<u8>, usize>> = vec![AHashMap::new(); hold.len()];
        for i in 0..n_batches {
            let blob = batch.call1((i,))?;
            let blob = blob.cast::<PyBytes>()?;
            let bytes: &[u8] = blob.as_bytes();
            let found = py.detach(|| finder.hits(bytes)).map_err(|pos| {
                pyo3::exceptions::PyValueError::new_err(format!("ntbatch record at byte {pos} does not start with '>'"))
            })?;
            for (strand, idx, off) in found {
                let (_, e, j) = owner[&kept[idx]];
                let start = j as i64 - off as i64;
                match seen[e].get(&strand) {
                    Some(&at) => hits[e][at].1 = start,
                    None => {
                        seen[e].insert(strand.clone(), hits[e].len());
                        hits[e].push((strand, start));
                    }
                }
            }
        }
        for e in 0..hold.len() {
            if !live[e] {
                continue;
            }
            if hits[e].len() > max_hits {
                live[e] = false;
                continue;
            }
            let tails: Vec<&[u8]> = hits[e]
                .iter()
                .filter_map(|(r, st)| overlap_tail(&hold[e], r, *st, read_agree))
                .collect();
            let mut new = extend_core(&tails, follows[e], depth, agree);
            // growth that runs back into the row: the k-mer may start in bases earlier rounds grew,
            // so their last k - 1 are checked with this round's
            let from = bases[e].max(hold[e].len().saturating_sub(k - 1));
            let old = hold[e].len() - from;
            let mut both = hold[e][from..].to_vec();
            both.extend_from_slice(&new);
            let cut = windows(&both, k)
                .into_iter()
                .position(|w| w.is_some_and(|km| row_kmers[edge_rows[e]].contains(&km)));
            if let Some(c) = cut {
                if c < old {
                    hold[e].truncate(from + c);
                    live[e] = false;
                    continue;
                }
                new.truncate(c - old);
            }
            if new.is_empty() {
                live[e] = false;
                continue;
            }
            hold[e].extend_from_slice(&new);
            if cut.is_some() || hold[e].len() - bases[e] >= max_grow {
                live[e] = false;
            }
        }
    }
    Ok(hold
        .into_iter()
        .zip(bases)
        .map(|(h, b)| String::from_utf8_lossy(&h[b..]).into_owned())
        .collect())
}
