# Rust Audit — SAPPHYRE Performance Findings

Audited every Rust file in `src/`. Findings are ranked by expected wallclock impact at N=100k candidates / M=3000 columns. Honest scaling bugs only — micro-optimizations dropped.

Files reviewed:
- `src/lib.rs` (551)
- `src/exon_dp.rs` (1951)
- `src/column_cull.rs` (1674)
- `src/consensus.rs` (412)
- `src/dedupe.rs` (259)
- `src/aligner.rs` (199)
- `src/translate.rs` (158)
- `src/flexcull.rs` (40)
- `src/identity.rs` (16)
- `src/overlap.rs` (18)
- `src/interval_tree.rs` (59)

---

## HIGH IMPACT

### 1. `cands.iter().find()` per MXE slot — `src/exon_dp.rs:1460-1461`

```rust
let left_cand = cands.iter().find(|c| c.nf == *left_node);
let right_cand = cands.iter().find(|c| c.nf == *right_node);
```

**Pattern.** Inside the MXE region scan loop over modular slots, two linear scans of `cands` per slot, keyed by `c.nf`.

**Scaling.** Outer loop runs once per modular slot. Per gene: K_slots × 2 × N candidates. For a fat gene with N=100k cands and K_slots=20, that's 4M string-compare ops per gene just here. Scales linearly with both — i.e. O(K × N).

**Fix.** Build `cand_by_nf: HashMap<&str, &Cand>` once at the top of `process_gene` (line ~990). The outer per-cluster loop already does `cluster_node_fields.contains(&c.nf)` filtering, which can use the same map. O(1) lookups instead of O(N) scans.

**Confidence: HIGH.** Textbook "indexable repeated find."

---

### 2. Per-cluster filter scans all candidates — `src/exon_dp.rs:1180-1185`

```rust
let mut aa_subset: Vec<usize> = cands
    .iter()
    .enumerate()
    .filter(|(_, c)| cluster_node_fields.contains(&c.nf))
    .map(|(i, _)| i)
    .collect();
```

**Pattern.** Outer loop iterates over `all_clusters` (K clusters per gene). For each, scans every candidate in `cands` to filter by HashSet membership.

**Scaling.** O(K × N) per gene. K is typically 5–20 clusters; N=100k for fat genes → 1–2M HashSet contains per gene. Each contains is ~30 ns, so ~50 ms per gene before the actual gap work even starts.

**Fix.** Build `cands_by_nf: HashMap<String, Vec<usize>>` once. Then per cluster, iterate `node_tokens` and gather candidate indices via `.get()`. O(N) total across all clusters instead of O(K × N).

**Confidence: HIGH.** Same `cand_by_nf` index that fixes finding 1 also fixes this.

---

### 3. Same identical filter pattern in `flank_scan_gene` — `src/exon_dp.rs:568-574`

```rust
let mut aa_subset: Vec<&FlankNode> = aa_nodes
    .iter()
    .filter(|n| {
        let nf = n.header.split('|').nth(3).unwrap_or("");
        cluster_set.contains(nf)
    })
    .collect();
```

**Pattern.** Same as #2, plus a fresh `header.split('|').nth(3)` per node per cluster.

**Scaling.** K_clusters × N_nodes per gene, with an inner `String::split + nth` allocation each time. For K=10, N=100k: 1M splits per gene. The split itself is small (~50ns), so ~50 ms per gene.

**Fix.** Pre-extract `nf` for each `FlankNode` once into the struct (or build `aa_nodes_by_nf: HashMap<&str, Vec<usize>>` once before the cluster loop).

**Confidence: HIGH.**

---

### 4. Per-record bool-mask write touches every column — `src/column_cull.rs:1015-1022, 1117-1124, 1228-1240`

Three near-identical hot loops:

```rust
// edge mask write (1015-1022)
for (h, _) in &normalized {
    if h.ends_with(ref_suffix) { continue; }
    if let Some(rv) = removed_columns.get_mut(h) {
        for i in 0..orig_aln_len {
            if edge_mask[i] { rv[i] = true; }
        }
    }
}
```

**Pattern.** For each non-ref record (N records), iterate every alignment column (M) and conditionally mark the per-record `removed_columns` Vec<bool>. Repeats for `gap_cull_mask` (1117) and `global_empty_mask` (1228).

**Scaling.** Three passes of O(N × M) per gene. For N=100k, M=3000: each pass is 300M iterations. Three passes = 900M tight-loop iterations per gene. Even at ~1 ns per iteration that's roughly 1 second/gene; across many genes this is the dominant cost in `cull_columns` when N is large.

**Fix.** Precompute the active indices once: `let edge_indices: Vec<usize> = (0..orig_aln_len).filter(|&i| edge_mask[i]).collect();`. Then per record: `for &i in &edge_indices { rv[i] = true; }`. Edge masks are typically sparse (only the leading/trailing edges flip, so K << M). Same fix applies to all three sites. 5–30× speedup on these loops.

**Confidence: HIGH.**

---

### 5. PSSM stored as `HashMap<usize, ...>` keyed by column index — `src/column_cull.rs:81-119, 124-148`

```rust
fn build_blosum_pssm(...) -> HashMap<usize, (HashMap<u8, f64>, f64)> { ... }

fn score_fragment(seq: &[u8], pssm: &HashMap<usize, (...)>, start: usize, end: usize) -> ... {
    for i in start..end {
        ...
        if let Some((scores, col_max)) = pssm.get(&i) { ... }
    }
}
```

**Pattern.** `pssm` is a HashMap keyed by `col: usize`. Each call to `score_fragment` iterates over a column range and does `pssm.get(&i)` per column — a hash lookup keyed by an integer that's already a dense index. Plus the inner `scores: HashMap<u8, f64>` keyed by amino acid byte (only 20 possible values).

**Scaling.** `mask_stop_blosum` calls `score_fragment` twice per stop codon per record. For N=100k records × ~M=3000 cols visited (worst case) × 2 calls = 600M HashMap lookups per gene. ~30ns each = ~18s per gene.

**Fix.** Replace the outer `HashMap<usize, ...>` with `Vec<Option<(...)>>` indexed by column — direct index access is single-digit ns. Replace inner `HashMap<u8, f64>` with `[f64; 20]` indexed by the same `idx(c)` function already in `blosum62`. Both are keyed by tiny dense ranges; HashMap is the wrong structure.

**Confidence: HIGH.**

---

### 6. Per-record sequence index mapping rebuilt with HashMap-style sparse pattern — `src/column_cull.rs:1131-1141, 1247-1267`

```rust
let after_gap_cull: Vec<(String, Vec<u8>)> = if current_to_orig.len() < orig_aln_len {
    after_structural
        .iter()
        .map(|(h, s)| {
            let new_seq: Vec<u8> = current_to_orig.iter().map(|&i| s[i]).collect();
            (h.clone(), new_seq)
        })
        .collect()
}
```

**Pattern.** Per record: rebuild a Vec by gathering bytes at `current_to_orig` indices. Plus the `(h.clone(), new_seq)` clones the header String once per record (also true at 1247-1267). For N=100k records, that's 100k String::clone calls — header strings are typically 50-100 bytes each, so ~10 MB of unnecessary heap traffic.

**Scaling.** O(N × K) where K = surviving columns. The byte-gather work itself is necessary, but the header.clone() for each record is avoidable. Across multiple stages (`after_structural` → `after_gap_cull` → `mask_retained_introns` → `mask_stop_blosum` → `masked_records` → `trimmed_records`), every record gets cloned 5+ times.

**Fix.** Switch the pipeline to take `&[(String, Vec<u8>)]` and produce `Vec<Vec<u8>>` (sequences only) at intermediate stages, joining headers back at the end. Or move headers out of records once and pass index arrays. Avoids ~5× redundant String::clone per record.

**Confidence: MEDIUM-HIGH.** Real but constant-factor; not a complexity bug.

---

## MEDIUM IMPACT

### 7. `cluster_node_fields.iter().any()` instead of indexed lookup — `src/exon_dp.rs:1280-1297`

```rust
let foreign: Vec<(usize, usize)> = sibling_intervals
    .iter()
    .filter(|(sscaf, ss, se)| {
        *sscaf == ga.scaffold
            && *ss >= rs
            && *se <= re
            && !cluster_node_fields.iter().any(|nf| {
                if let Some(entry) = gff.get(nf.as_str()) {
                    entry.scaffold == *sscaf && entry.start == *ss && entry.end == *se
                } else { false }
            })
    })
    ...
```

**Pattern.** Inside the gap loop, for each sibling interval, linearly scan `cluster_node_fields` (a HashSet but used as a Vec via `.iter()`), doing a `gff.get` + 3 field comparisons per element.

**Scaling.** Per MXE gap: S × C ops where S ≈ 5–50 sibling intervals, C ≈ 5–50 cluster nodes. Per gene: K_mxe × N_gaps × S × C ≈ 5 × 30 × 20 × 20 = 60k ops. Across 100k genes ≈ 6B ops, with each op being a HashMap get. Real but bounded.

**Fix.** Hoist outside the gap loop: pre-build `cluster_intervals: HashSet<(scaffold, usize, usize)>` for this cluster once. Then the filter check is O(1).

**Confidence: HIGH that the pattern is suboptimal, MEDIUM on real impact.** Constant factor matters — `gff.get` is a string hash, not free.

---

### 8. Ref consensus stored as `HashMap<usize, Vec<u8>>` keyed by column — `src/exon_dp.rs:990-1019, 497-528`

```rust
let mut rc: HashMap<usize, Vec<u8>> = HashMap::new();
...
for (i, &bp) in sb.iter().enumerate() {
    let ent = rc.entry(i).or_default();
    ent.push(if i >= s && i <= e { bp } else { b'-' });
}
```

**Pattern.** `rc` (and `ref_consensus` in `flank_scan_gene` at lines 510-515) is a HashMap keyed by column index. Per ref × per col: a `.entry(i).or_default()` hash lookup. This is the exact same anti-pattern as finding 5 — dense integer key in a HashMap.

**Scaling.** R refs × M cols hash operations per gene. R=50, M=3000 → 150k hash ops per gene before the work even starts. ~5 ms/gene wasted per call site, two call sites = 10 ms/gene.

**Fix.** Replace with `Vec<Vec<u8>>` indexed by column. `build_is_ref_gap` already iterates by column index — make that explicit.

**Confidence: HIGH.**

---

### 9. `delete_empty_columns_pairs` clones buffer per record — `src/lib.rs:180-191`

```rust
let mut buf = Vec::with_capacity(kept_indices.len());
for (header, seq) in records {
    let bytes = seq.as_bytes();
    buf.clear();
    for &i in &kept_indices {
        if i < bytes.len() {
            buf.push(bytes[i]);
        }
    }
    out.push((header, String::from_utf8(buf.clone()).unwrap()));
}
```

**Pattern.** `buf.clone()` clones the byte vector every record so `buf` can be reused; `String::from_utf8(...).unwrap()` then validates UTF-8 of ASCII alignment data.

**Scaling.** N records × K cols. N=100k records, K=3000 cols → 300M bytes cloned + 300M bytes UTF-8 validated per call. Real wallclock cost on hot paths.

**Fix.** `let buf_out = std::mem::replace(&mut buf, Vec::with_capacity(kept_indices.len()));` then `unsafe { String::from_utf8_unchecked(buf_out) }`. Saves both the clone and the validation.

**Confidence: HIGH.**

---

### 10. `blosum62_distance` rebuilds 27-element HashSet per call — `src/lib.rs:64-98, 101-141`

```rust
fn blosum62_distance(one: String, two: String) -> f64 {
    ...
    let allowed: HashSet<u8> = HashSet::from([
        65, 84, 67, 71, 73, 68, 82, 80, 87, 77, 69, 81, 83, 72, 86, 76, 75, 70, 89, 78, 88, 90, 74,
        66, 79, 85, 42,
    ]);
    for i in 0..length {
        ...
        if !(allowed.contains(&char1)) { panic!(...) }
        if !(allowed.contains(&char2)) { panic!(...) }
        ...
    }
}
```

**Pattern.** Allocates a 27-element HashSet on every call, then does `contains()` checks per character.

**Scaling.** If called from per-record Python in tight loops (ref-to-ref blosum is N²/2 calls), the allocation alone is the dominant cost — `HashSet::from` for 27 elements is ~500 ns. For pairwise on 1k refs that's 500k allocations = ~250 ms wasted per gene. Larger ref sets are quadratic in this cost.

**Fix.** Replace with `static ALLOWED: [bool; 256]` lookup table built at compile time, or a `match` arm in a `#[inline]` function. O(1) lookup with zero allocation.

**Confidence: HIGH.** Same fix needed for `blosum62_candidate_to_reference`.

---

## LOW IMPACT (dropped — not theatrical, but small)

- **`bio_revcomp` UTF-8 validation roundtrip** (`src/lib.rs:22-24`) — Per-call O(L) validation pass; ~10 ns/byte on 500 bp sequences = 5 µs/call. Real but small. Use `String::from_utf8_unchecked` since `bio::alphabets::dna::revcomp` preserves ASCII.
- **MXE scan ORF projection over `known_intervals`** (`src/exon_dp.rs:1546-1577`) — O(orfs × known) per slot. Counts: 100 × 20 = 2k per slot. Fine.
- **Per-block `(bs..be).filter(|&i| ...).count()` in `mask_retained_introns`** (`src/column_cull.rs:472, 491, 494`) — Block sizes are small, called per (record × block). Tolerable.
- **`build_is_ref_gap` inner filter** (`src/exon_dp.rs:182-193`) — O(R × M) once per gene. Necessary work.

---

## Disputed Findings From The Earlier Pass

The first audit's "Triple-nested `.contains()` on Vec at exon_dp.rs:1140-1142" claim was **wrong**. `children` there is `&Vec<HashSet<String>>` (constructed at `src/exon_dp.rs:1121` as `iso_sets: HashMap<String, Vec<HashSet<String>>>`), so `child.contains(*node)` is O(1) HashSet lookup, not O(child_len). The total cost is O(B × C) where B=base_fields and C=children, both small (<50). Not a finding.

---

## Files With Nothing To Flag

- `src/consensus.rs` — `_dumb_consensus*` are O(N × L) by necessity; inner alphabet loop is bounded at 27. Clean.
- `src/dedupe.rs` — Open-addressing hash table with proper grow strategy. Clean.
- `src/aligner.rs` — Already uses `String::from_utf8_unchecked` on the hot return path. Clean.
- `src/translate.rs` — Static codon tables, single-pass `chunks_exact(3)`. Clean.
- `src/flexcull.rs`, `src/identity.rs`, `src/overlap.rs`, `src/interval_tree.rs` — All small, O(N) at worst, HashSet membership checks where appropriate.

---

## Ranking By Expected Wallclock Impact

| # | File:Line | Impact | Confidence |
|---|---|---|---|
| 1 | `exon_dp.rs:1460-1461` — `cands.iter().find()` × 2 per slot | HIGH | HIGH |
| 4 | `column_cull.rs:1015,1117,1228` — N×M bool-mask writes ×3 | HIGH | HIGH |
| 5 | `column_cull.rs:81-148` — PSSM as `HashMap<usize,...>` | HIGH | HIGH |
| 2 | `exon_dp.rs:1180-1185` — per-cluster filter scans all cands | HIGH | HIGH |
| 3 | `exon_dp.rs:568-574` — same in flank_scan_gene + per-iter split | HIGH | HIGH |
| 6 | `column_cull.rs:1131,1247` — repeated `header.clone()` per stage | MEDIUM-HIGH | MEDIUM-HIGH |
| 8 | `exon_dp.rs:990, 510` — ref_consensus as `HashMap<usize,Vec<u8>>` | MEDIUM | HIGH |
| 9 | `lib.rs:180-191` — `delete_empty_columns_pairs` clone+validate per record | MEDIUM | HIGH |
| 7 | `exon_dp.rs:1280-1297` — `cluster_node_fields.iter().any()` in gap loop | MEDIUM | MEDIUM |
| 10 | `lib.rs:64,101` — blosum62 allocs HashSet per call | MEDIUM | HIGH |
