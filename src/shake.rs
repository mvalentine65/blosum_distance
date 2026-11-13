//! Read shaking for SAPPHYRE align: the reads of one sample are joined into contigs
//! by sequence overlap, so reads that share sequence can share columns. Gives
//! exactly what aligner._build_contigs gives, ties included.
use ahash::AHashMap;
use pyo3::prelude::*;
use std::hash::Hash;

/// A counter that keeps first-insertion order, so most_common breaks ties as
/// Python's Counter does.
struct Ordered<K> {
    keys: Vec<K>,
    counts: Vec<u32>,
    at: AHashMap<K, usize>,
}

impl<K: Copy + Eq + Hash> Ordered<K> {
    fn new() -> Self {
        Ordered { keys: Vec::new(), counts: Vec::new(), at: AHashMap::new() }
    }

    fn add(&mut self, k: K) {
        match self.at.get(&k) {
            Some(&i) => self.counts[i] += 1,
            None => {
                self.at.insert(k, self.keys.len());
                self.keys.push(k);
                self.counts.push(1);
            }
        }
    }

    /// Keys by count, highest first; equal counts keep insertion order.
    fn most_common(&self, n: usize) -> Vec<K> {
        let mut idx: Vec<usize> = (0..self.keys.len()).collect();
        idx.sort_by(|&a, &b| self.counts[b].cmp(&self.counts[a]));
        idx.into_iter().take(n).map(|i| self.keys[i]).collect()
    }
}

/// Seed of `s` at `i` packed into a u64 (seed length at most 8).
#[inline]
fn seed(s: &[u8], i: usize, k: usize) -> u64 {
    s[i..i + k].iter().fold(0u64, |acc, &b| (acc << 8) | b as u64)
}

fn seeds_of(s: &[u8], k: usize) -> Vec<u64> {
    if s.len() < k {
        return Vec::new();
    }
    (0..=s.len() - k).map(|i| seed(s, i, k)).collect()
}

/// Overlap length and differences of `v` placed so that u[i] faces v[i + d].
fn mismatches(u: &[u8], v: &[u8], d: i64) -> (i64, i64) {
    let lo = (-d).max(0);
    let hi = (u.len() as i64).min(v.len() as i64 - d);
    let mut mm = 0;
    for i in lo..hi {
        mm += (u[i as usize] != v[(i + d) as usize]) as i64;
    }
    (hi - lo, mm)
}

/// Python's max(2, int(rate * n)).
#[inline]
fn allowed(rate: f64, n: i64) -> i64 {
    2.max((rate * n as f64) as i64)
}

struct Links {
    parent: Vec<usize>,
    off: Vec<i64>,
}

impl Links {
    fn find(&mut self, mut x: usize) -> usize {
        let mut path = Vec::new();
        while self.parent[x] != x {
            path.push(x);
            x = self.parent[x];
        }
        let mut acc = 0;
        for &y in path.iter().rev() {
            acc += self.off[y];
            self.off[y] = acc;
            self.parent[y] = x;
        }
        x
    }

    fn position(&mut self, x: usize) -> (usize, i64) {
        let root = self.find(x);
        (root, if x != root { self.off[x] } else { 0 })
    }
}

/// Contig sequences, and (contig, offset) for every read in input order.
type Shaken = (Vec<String>, Vec<(usize, i64)>);

/// aligner._build_contigs: shake reads into contigs. Returns the contig sequences and,
/// for every read in input order, (contig, offset of the read in it).
#[pyfunction]
#[pyo3(signature = (reads, seed_len, min_overlap, mismatch_rate, edge_mismatch_rate, partners, seed_sample))]
#[allow(clippy::too_many_arguments)]
pub fn shake_contigs(
    py: Python<'_>,
    reads: Vec<String>,
    seed_len: usize,
    min_overlap: i64,
    mismatch_rate: f64,
    edge_mismatch_rate: f64,
    partners: usize,
    seed_sample: usize,
) -> PyResult<Shaken> {
    if seed_len == 0 || seed_len > 8 {
        return Err(pyo3::exceptions::PyValueError::new_err("seed_len must be 1 to 8"));
    }
    let s = seed_len;
    Ok(py.detach(|| {
        // Distinct reads, longest first; equal lengths keep first appearance.
        let mut first: AHashMap<&[u8], usize> = AHashMap::new();
        let mut distinct: Vec<&[u8]> = Vec::new();
        let mut dup: Vec<u32> = Vec::new();
        let mut of_read = Vec::with_capacity(reads.len());
        for r in &reads {
            let b = r.as_bytes();
            let d = *first.entry(b).or_insert_with(|| {
                distinct.push(b);
                dup.push(0);
                distinct.len() - 1
            });
            dup[d] += 1;
            of_read.push(d);
        }
        let mut order: Vec<usize> = (0..distinct.len()).collect();
        order.sort_by(|&a, &b| distinct[b].len().cmp(&distinct[a].len()));
        let n = order.len();
        let seqs_r: Vec<&[u8]> = order.iter().map(|&d| distinct[d]).collect();
        let weight: Vec<u32> = order.iter().map(|&d| dup[d]).collect();
        let mut rank = vec![0usize; n];
        for (r, &d) in order.iter().enumerate() {
            rank[d] = r;
        }

        let mut seeds: AHashMap<u64, Vec<(usize, i64)>> = AHashMap::new();
        for (r, u) in seqs_r.iter().enumerate() {
            for (i, k) in seeds_of(u, s).into_iter().enumerate() {
                seeds.entry(k).or_default().push((r, i as i64));
            }
        }

        // Read level: each read links to the reads sharing the most seeds with it.
        let mut links = Links { parent: (0..n).collect(), off: vec![0; n] };
        for a in 0..n {
            let u = seqs_r[a];
            let mut hits: Ordered<(usize, i64)> = Ordered::new();
            if u.len() >= s {
                for i in (0..=u.len() - s).step_by(3) {
                    if let Some(list) = seeds.get(&seed(u, i, s)) {
                        for &(b, j) in list.iter().take(seed_sample) {
                            if b != a {
                                hits.add((b, j - i as i64));
                            }
                        }
                    }
                }
            }
            let mut made = 0;
            for (b, d) in hits.most_common(usize::MAX) {
                if made >= partners {
                    break;
                }
                let v = seqs_r[b];
                let (olap, mm) = mismatches(u, v, d);
                if olap < min_overlap.min(u.len() as i64).min(v.len() as i64) || mm > allowed(edge_mismatch_rate, olap) {
                    continue;
                }
                let (ra, pa) = links.position(a);
                let (rb, pb) = links.position(b);
                // seed u[i] == v[i + d], so v starts d residues before u
                if ra == rb {
                    made += (pb == pa - d) as usize;
                    continue;
                }
                links.parent[rb] = ra;
                links.off[rb] = pa - d - pb;
                made += 1;
            }
        }

        // Contigs from the linked groups, consensus weighted by duplicate counts.
        let mut group_of: AHashMap<usize, usize> = AHashMap::new();
        let mut groups: Vec<Vec<(usize, i64)>> = Vec::new();
        for r in 0..n {
            let (root, pos) = links.position(r);
            let g = *group_of.entry(root).or_insert_with(|| {
                groups.push(Vec::new());
                groups.len() - 1
            });
            groups[g].push((r, pos));
        }
        let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(groups.len());
        let mut place: Vec<(usize, i64)> = vec![(0, 0); n];
        let mut members: Vec<Vec<usize>> = Vec::with_capacity(groups.len());
        for g in &groups {
            let lo = g.iter().map(|&(_, p)| p).min().unwrap();
            let length = g.iter().map(|&(r, p)| p - lo + seqs_r[r].len() as i64).max().unwrap() as usize;
            let mut counts: Vec<Vec<(u8, u32)>> = vec![Vec::new(); length];
            for &(r, p) in g {
                let w = weight[r];
                for (i, &ch) in seqs_r[r].iter().enumerate() {
                    let col = &mut counts[(p - lo) as usize + i];
                    match col.iter_mut().find(|e| e.0 == ch) {
                        Some(e) => e.1 += w,
                        None => col.push((ch, w)),
                    }
                }
            }
            let ci = seqs.len();
            seqs.push(
                counts
                    .iter()
                    .map(|col| {
                        let mut best = col[0];
                        for &e in &col[1..] {
                            if e.1 > best.1 {
                                best = e;
                            }
                        }
                        best.0
                    })
                    .collect(),
            );
            for &(r, p) in g {
                place[r] = (ci, p - lo);
            }
            members.push(g.iter().map(|&(r, _)| r).collect());
        }

        // Contig level: smallest first, a contig that fits a longer root contig moves into
        // its coordinates; passes repeat until nothing moves. A contig that found no fit is
        // checked again only when it or a candidate has grown since.
        let nc = seqs.len();
        let mut parent: Vec<usize> = (0..nc).collect();
        let mut owner: AHashMap<u64, Vec<usize>> = AHashMap::new(); // every contig that ever held a seed
        let mut seed_at: AHashMap<u64, Vec<(usize, i64)>> = AHashMap::new(); // live contigs: (contig, position)
        fn index(
            b: usize,
            t: &[u8],
            s: usize,
            seed_at: &mut AHashMap<u64, Vec<(usize, i64)>>,
            owner: &mut AHashMap<u64, Vec<usize>>,
        ) {
            for (i, k) in seeds_of(t, s).into_iter().enumerate() {
                seed_at.entry(k).or_default().push((b, i as i64));
                owner.entry(k).or_default().push(b); // repeats removed when candidates are gathered
            }
        }
        fn unindex(b: usize, t: &[u8], s: usize, seed_at: &mut AHashMap<u64, Vec<(usize, i64)>>) {
            for k in seeds_of(t, s) {
                if let Some(list) = seed_at.get_mut(&k) {
                    list.retain(|e| e.0 != b);
                }
            }
        }
        for (b, t) in seqs.iter().enumerate() {
            index(b, t, s, &mut seed_at, &mut owner);
        }
        let mut stamp: i64 = 0;
        let mut grown = vec![0i64; nc];
        let mut checked = vec![-1i64; nc];
        let mut seen = vec![usize::MAX; nc];
        let mut epoch = 0usize; // one per candidate gathering, so a contig checked again starts clean
        let mut changed = true;
        while changed {
            changed = false;
            let mut roots: Vec<usize> = (0..nc).filter(|&i| parent[i] == i).collect();
            roots.sort_by_key(|&i| seqs[i].len());
            for a in roots {
                if parent[a] != a {
                    continue;
                }
                let u = seqs[a].clone();
                let lu = u.len() as i64;
                let useeds = seeds_of(&u, s);
                let mut cand: Vec<usize> = Vec::new();
                epoch += 1;
                for k in &useeds {
                    for &b in owner.get(k).into_iter().flatten() {
                        if seen[b] != epoch {
                            seen[b] = epoch;
                            cand.push(b);
                        }
                    }
                }
                cand.sort_unstable();
                if grown[a] <= checked[a] && cand.iter().all(|&b| grown[b] <= checked[a]) {
                    continue;
                }
                checked[a] = stamp;
                let mut offs: AHashMap<usize, Ordered<i64>> = AHashMap::new();
                let mut ok: AHashMap<usize, bool> = AHashMap::new();
                for (i, k) in useeds.iter().enumerate() {
                    if let Some(list) = seed_at.get(k) {
                        for &(b, j) in list {
                            let good = *ok
                                .entry(b)
                                .or_insert_with(|| b != a && parent[b] == b && seqs[b].len() as i64 >= lu);
                            if good {
                                offs.entry(b).or_insert_with(Ordered::new).add(j - i as i64);
                            }
                        }
                    }
                }
                // best: (contig, offset, differences, overlap); equal fits: lowest contig wins
                let mut best: Option<(usize, i64, i64, i64)> = None;
                for &b in &cand {
                    if b == a || parent[b] != b || (seqs[b].len() as i64) < lu {
                        continue;
                    }
                    let Some(o) = offs.get(&b) else { continue };
                    let cons = &seqs[b];
                    let mut got: Option<(i64, i64, i64)> = None;
                    for d in o.most_common(3) {
                        let (m, mm) = mismatches(&u, cons, d);
                        if m < min_overlap.min(lu) {
                            continue;
                        }
                        if mm <= allowed(mismatch_rate, m) && got.is_none_or(|g| (mm, -m) < (g.1, -g.2)) {
                            got = Some((d, mm, m));
                        }
                    }
                    if let Some((d, mm, m)) = got {
                        if best.is_none_or(|x| (mm, -m) < (x.2, -x.3)) {
                            best = Some((b, d, mm, m));
                        }
                    }
                }
                let Some((b, mut d, _, _)) = best else { continue };
                stamp += 1;
                let lb = seqs[b].len() as i64;
                let left = (-d).max(0);
                let right = (d + lu - lb).max(0);
                if left > 0 || right > 0 {
                    let old = std::mem::take(&mut seqs[b]);
                    unindex(b, &old, s, &mut seed_at);
                    let mut t = u[..left as usize].to_vec();
                    t.extend_from_slice(&old);
                    t.extend_from_slice(&u[(lu - right) as usize..]);
                    index(b, &t, s, &mut seed_at, &mut owner);
                    seqs[b] = t;
                    grown[b] = stamp;
                    for &r in &members[b] {
                        place[r].1 += left;
                    }
                    d += left;
                }
                let moved = std::mem::take(&mut members[a]);
                for &r in &moved {
                    place[r] = (b, place[r].1 + d);
                }
                members[b].extend(moved);
                parent[a] = b;
                unindex(a, &u, s, &mut seed_at);
                changed = true;
            }
        }
        let mut new = vec![usize::MAX; nc];
        let mut out = Vec::new();
        for i in 0..nc {
            if parent[i] == i {
                new[i] = out.len();
                out.push(String::from_utf8_lossy(&seqs[i]).into_owned());
            }
        }
        let at = of_read
            .iter()
            .map(|&d| {
                let (c, o) = place[rank[d]];
                (new[c], o)
            })
            .collect();
        (out, at)
    }))
}
