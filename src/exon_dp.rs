use flate2::read::GzDecoder;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use rocksdb::{Options, DB};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;
use std::io::Read;
use std::path::Path;

/// Fraction of refs that must be gap/space for a column to be masked out.
const REF_GAP_THR: f64 = 0.67;
/// Minimum effective (masked) gap columns to trigger a scan.
const MINIMUM_GAP_AA: usize = 10;
/// Maximum effective (masked) gap columns; larger gaps are skipped.
const MAX_GAP_AA: usize = 250;
/// Genomic search region cap for flank scans (bp).
const FLANK_BP: usize = 15000;

fn codon_to_aa(c1: u8, c2: u8, c3: u8) -> u8 {
    match (c1, c2, c3) {
        (b'T', b'T', b'T') | (b'T', b'T', b'C') => b'F',
        (b'T', b'T', b'A') | (b'T', b'T', b'G') => b'L',
        (b'C', b'T', b'T') | (b'C', b'T', b'C') | (b'C', b'T', b'A') | (b'C', b'T', b'G') => b'L',
        (b'A', b'T', b'T') | (b'A', b'T', b'C') | (b'A', b'T', b'A') => b'I',
        (b'A', b'T', b'G') => b'M',
        (b'G', b'T', b'T') | (b'G', b'T', b'C') | (b'G', b'T', b'A') | (b'G', b'T', b'G') => b'V',
        (b'T', b'C', b'T') | (b'T', b'C', b'C') | (b'T', b'C', b'A') | (b'T', b'C', b'G') => b'S',
        (b'C', b'C', b'T') | (b'C', b'C', b'C') | (b'C', b'C', b'A') | (b'C', b'C', b'G') => b'P',
        (b'A', b'C', b'T') | (b'A', b'C', b'C') | (b'A', b'C', b'A') | (b'A', b'C', b'G') => b'T',
        (b'G', b'C', b'T') | (b'G', b'C', b'C') | (b'G', b'C', b'A') | (b'G', b'C', b'G') => b'A',
        (b'T', b'A', b'T') | (b'T', b'A', b'C') => b'Y',
        (b'T', b'A', b'A') | (b'T', b'A', b'G') => b'*',
        (b'C', b'A', b'T') | (b'C', b'A', b'C') => b'H',
        (b'C', b'A', b'A') | (b'C', b'A', b'G') => b'Q',
        (b'A', b'A', b'T') | (b'A', b'A', b'C') => b'N',
        (b'A', b'A', b'A') | (b'A', b'A', b'G') => b'K',
        (b'G', b'A', b'T') | (b'G', b'A', b'C') => b'D',
        (b'G', b'A', b'A') | (b'G', b'A', b'G') => b'E',
        (b'T', b'G', b'T') | (b'T', b'G', b'C') => b'C',
        (b'T', b'G', b'A') => b'*',
        (b'T', b'G', b'G') => b'W',
        (b'C', b'G', b'T') | (b'C', b'G', b'C') | (b'C', b'G', b'A') | (b'C', b'G', b'G') => b'R',
        (b'A', b'G', b'T') | (b'A', b'G', b'C') => b'S',
        (b'A', b'G', b'A') | (b'A', b'G', b'G') => b'R',
        (b'G', b'G', b'T') | (b'G', b'G', b'C') | (b'G', b'G', b'A') | (b'G', b'G', b'G') => b'G',
        _ => 0,
    }
}

fn bio_revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            other => other,
        })
        .collect()
}

fn find_index_pair(seq: &[u8], gap: u8) -> (usize, usize) {
    let s = seq.iter().position(|&c| c != gap).unwrap_or(0);
    let e = seq.iter().rposition(|&c| c != gap).unwrap_or(0);
    (s, e)
}

fn parse_node_field(h: &str) -> &str {
    h.split('|').nth(3).unwrap_or("")
}

// ---------------------------------------------------------------------------
// Flank scan: ORF-based flank extension
// ---------------------------------------------------------------------------

struct FlankNode {
    header: String,
    frame: i32,
    start: usize,
    end: usize,
}

struct FlankCandidate {
    aa_seq: String,
    nt_seq: String,
    frame: i32,
    scaffold: String,
    hit_start: usize,
    hit_end: usize,
    strand: String,
    node_name: String,
    is_leading: bool,
    gap_start: usize,
    gap_end: usize,
    cluster_key: String,
}

/// Build sorted exon interval lists per scaffold from GFF entries.
fn build_scaffold_intervals(
    gff: &HashMap<String, GffEntry>,
) -> HashMap<String, Vec<(usize, usize, String)>> {
    let mut intervals: HashMap<String, Vec<(usize, usize, String)>> = HashMap::new();
    for (name, entry) in gff {
        intervals
            .entry(entry.scaffold.clone())
            .or_default()
            .push((entry.start, entry.end, name.clone()));
    }
    for v in intervals.values_mut() {
        v.sort();
    }
    intervals
}

/// Find the end of the nearest exon LEFT of `position` on this scaffold.
/// Intervals are sorted by start, so binary search for the insertion point.
fn find_left_bound(
    scaffold: &str,
    position: usize,
    intervals: &HashMap<String, Vec<(usize, usize, String)>>,
) -> usize {
    if let Some(ivs) = intervals.get(scaffold) {
        // Find the last interval whose start < position
        let idx = ivs.partition_point(|&(s, _, _)| s < position);
        // Walk backwards from idx to find the best end < position
        let mut best = 0usize;
        for i in (0..idx).rev() {
            if ivs[i].1 < position {
                best = best.max(ivs[i].1);
                break; // sorted by start, so earlier intervals have smaller ends in practice
            }
        }
        best
    } else {
        0
    }
}

/// Find the start of the nearest exon RIGHT of `position` on this scaffold.
/// Intervals are sorted by start, so binary search for first start > position.
fn find_right_bound(
    scaffold: &str,
    position: usize,
    intervals: &HashMap<String, Vec<(usize, usize, String)>>,
    scaffold_len: usize,
) -> usize {
    if let Some(ivs) = intervals.get(scaffold) {
        let idx = ivs.partition_point(|&(s, _, _)| s <= position);
        if idx < ivs.len() {
            return ivs[idx].0;
        }
    }
    scaffold_len
}

/// Count reference gap columns in [gap_start..gap_end).
fn flank_count_ref_gaps(
    gap_start: usize,
    gap_end: usize,
    ref_consensus: &HashMap<usize, Vec<u8>>,
    ref_gap_thresh: f64,
) -> (HashSet<usize>, Vec<usize>, usize) {
    let mut ref_gaps = HashSet::new();
    let mut insert_at = Vec::new();
    let mut longest_consecutive: usize = 0;
    let mut consecutive_non_gap: usize = 0;
    for (i, x) in (gap_start..gap_end).enumerate() {
        let col = ref_consensus.get(&x);
        let is_gap = if let Some(col) = col {
            let gap_count = col.iter().filter(|&&c| c == b'-' || c == b' ').count();
            (gap_count as f64 / col.len() as f64) >= ref_gap_thresh
        } else {
            true
        };
        if is_gap {
            ref_gaps.insert(x);
            insert_at.push(i);
            longest_consecutive = longest_consecutive.max(consecutive_non_gap);
            consecutive_non_gap = 0;
        } else {
            consecutive_non_gap += 1;
        }
    }
    longest_consecutive = longest_consecutive.max(consecutive_non_gap);
    (ref_gaps, insert_at, longest_consecutive)
}

/// Check whether a gap qualifies for scanning. Returns
/// `Some((effective_gap, scaled_min_aa))` if it does, `None` otherwise.
fn qualify_gap(
    gap_start: usize,
    gap_end: usize,
    ref_consensus: &HashMap<usize, Vec<u8>>,
) -> Option<(usize, usize)> {
    let (ref_gaps, _, _) =
        flank_count_ref_gaps(gap_start, gap_end, ref_consensus, REF_GAP_THR);
    let total_cols = gap_end - gap_start;
    let effective_gap = total_cols - ref_gaps.len();
    if effective_gap < MINIMUM_GAP_AA || effective_gap >= MAX_GAP_AA {
        return None;
    }
    let scaled_min_aa = MINIMUM_GAP_AA.max(effective_gap / 2);
    Some((effective_gap, scaled_min_aa))
}

/// Translate nucleotide bytes to amino acid string (standard code).
fn flank_translate(nt: &[u8]) -> Vec<u8> {
    let mut aa = Vec::with_capacity(nt.len() / 3);
    let mut i = 0;
    while i + 2 < nt.len() {
        let c = codon_to_aa(
            nt[i].to_ascii_uppercase(),
            nt[i + 1].to_ascii_uppercase(),
            nt[i + 2].to_ascii_uppercase(),
        );
        aa.push(if c == 0 { b'X' } else { c });
        i += 3;
    }
    aa
}

/// Trim sequence to multiple of 3.
fn trim_to_codon(seq: &[u8]) -> &[u8] {
    let r = seq.len() % 3;
    if r > 0 { &seq[..seq.len() - r] } else { seq }
}

/// Extract ORFs from a nucleotide sequence across all 3 reading frames.
/// Returns (aa_seq, nt_start, nt_end, frame) for each ORF >= min_aa residues.
/// ORFs are delimited by stop codons or sequence boundaries.
///
/// Leading flanks (require_start=true):
///   - Trim 5' to first M (ATG). No M = discard.
///
/// Results are deduplicated by AA containment (shorter contained in longer removed).
fn extract_flank_orfs(
    seq: &[u8],
    min_aa: usize,
    require_start: bool,
) -> Vec<(Vec<u8>, usize, usize, usize)> {
    let mut raw: Vec<(Vec<u8>, usize, usize, usize)> = Vec::new();
    let seq_len = seq.len();

    for frame in 0..3usize {
        if frame >= seq_len {
            continue;
        }
        let trimmed = trim_to_codon(&seq[frame..]);
        let protein = flank_translate(trimmed);

        // Collect stop-to-stop segments: (aa_start_inclusive, aa_end_exclusive)
        let mut segments: Vec<(usize, usize)> = Vec::new();
        let mut seg_start: usize = 0;
        for (i, &aa) in protein.iter().enumerate() {
            if aa == b'*' {
                if i > seg_start {
                    segments.push((seg_start, i));
                }
                seg_start = i + 1;
            }
        }
        if protein.len() > seg_start {
            segments.push((seg_start, protein.len()));
        }

        for (seg_start_aa, seg_end_aa) in segments {
            let mut aa_from = seg_start_aa;
            let aa_to = seg_end_aa;

            // --- Leading flank: trim 5' to first M ---
            if require_start {
                let slice = &protein[aa_from..aa_to];
                if let Some(m_pos) = slice.iter().position(|&a| a == b'M') {
                    aa_from += m_pos;
                } else {
                    continue;
                }
            }

            if aa_to > aa_from && aa_to - aa_from >= min_aa {
                let nt_start = frame + aa_from * 3;
                let nt_end = frame + aa_to * 3;
                raw.push((
                    protein[aa_from..aa_to].to_vec(),
                    nt_start,
                    nt_end,
                    frame,
                ));
            }
        }
    }

    if raw.is_empty() {
        return raw;
    }

    // Sort by AA length descending for containment check
    raw.sort_by(|a, b| b.0.len().cmp(&a.0.len()));

    let mut kept: Vec<(Vec<u8>, usize, usize, usize)> = Vec::new();
    for entry in raw {
        let dominated = kept.iter().any(|k| {
            let needle = &entry.0;
            let haystack = &k.0;
            if needle.len() > haystack.len() {
                return false;
            }
            haystack
                .windows(needle.len())
                .any(|w| w == needle.as_slice())
        });
        if !dominated {
            kept.push(entry);
        }
    }

    kept
}

/// Extract candidate ORFs from a genomic flank region.
/// Extracts candidate ORFs from a genomic flank region.
fn flank_extract_orfs(
    gap_start: usize,
    gap_end: usize,
    node: &FlankNode,
    is_leading: bool,
    log: &mut Vec<String>,
    ref_consensus: &HashMap<usize, Vec<u8>>,
    gff: &HashMap<String, GffEntry>,
    genome: &HashMap<String, Vec<u8>>,
    scaffold_intervals: &HashMap<String, Vec<(usize, usize, String)>>,
    cluster_key: &str,
    results: &mut Vec<FlankCandidate>,
) {
    let direction = if is_leading { "Left leading" } else { "Right trailing" };

    let Some((effective_gap, scaled_min_aa)) = qualify_gap(gap_start, gap_end, ref_consensus) else {
        return;
    };
    let total_cols = gap_end - gap_start;

    log.push(format!(
        "{} gap of {} ({}/{} cols effective/total)",
        direction, node.header, effective_gap, total_cols
    ));

    let node_field = parse_node_field(&node.header).to_string();
    let node_name = if node_field.contains("&&") {
        if is_leading {
            node_field.split("&&").next().unwrap_or(&node_field)
        } else {
            node_field.split("&&").last().unwrap_or(&node_field)
        }
        .split('_')
        .next()
        .unwrap_or(&node_field)
    } else {
        node_field.split('_').next().unwrap_or(&node_field)
    };

    let Some(gff_entry) = gff.get(&node_field).or_else(|| gff.get(node_name)) else {
        log.push(format!("No GFF entry for {}\n", node_field));
        return;
    };

    let Some(scaffold_seq) = genome.get(&gff_entry.scaffold) else {
        log.push(format!("Scaffold {} not loaded\n", gff_entry.scaffold));
        return;
    };
    let scaffold_len = scaffold_seq.len();

    let (seq, rs, re) = if is_leading {
        if gff_entry.strand == "+" {
            let left = find_left_bound(&gff_entry.scaffold, gff_entry.start, scaffold_intervals);
            let s = left.max(gff_entry.start.saturating_sub(1).saturating_sub(FLANK_BP));
            let e = gff_entry.start.saturating_sub(1);
            (scaffold_seq[s..e].to_vec(), s, e)
        } else {
            let right = find_right_bound(
                &gff_entry.scaffold,
                gff_entry.end,
                scaffold_intervals,
                scaffold_len,
            );
            let s = gff_entry.end;
            let e = right.min(gff_entry.end + FLANK_BP);
            (bio_revcomp(&scaffold_seq[s..e]), s, e)
        }
    } else {
        if gff_entry.strand == "+" {
            let right = find_right_bound(
                &gff_entry.scaffold,
                gff_entry.end,
                scaffold_intervals,
                scaffold_len,
            );
            let s = gff_entry.end;
            let e = right.min(gff_entry.end + FLANK_BP);
            (scaffold_seq[s..e].to_vec(), s, e)
        } else {
            let left = find_left_bound(&gff_entry.scaffold, gff_entry.start, scaffold_intervals);
            let s = left.max(gff_entry.start.saturating_sub(1).saturating_sub(FLANK_BP));
            let e = gff_entry.start.saturating_sub(1);
            (bio_revcomp(&scaffold_seq[s..e]), s, e)
        }
    };

    if seq.is_empty() {
        log.push("Empty flank region\n".to_string());
        return;
    }

    log.push(format!(
        "  Node: {}:{}-{} ({}) Search: {}:{}-{} len={}",
        gff_entry.scaffold,
        gff_entry.start,
        gff_entry.end,
        gff_entry.strand,
        gff_entry.scaffold,
        rs + 1,
        re,
        seq.len()
    ));

    // Extract ORFs from all 3 frames
    let orfs = extract_flank_orfs(&seq, scaled_min_aa, is_leading);
    if orfs.is_empty() {
        log.push(format!(
            "No ORFs >= {} AA found in flank region\n",
            scaled_min_aa
        ));
        return;
    }

    log.push(format!(
        "  {} candidate ORFs (min_aa={}, gap={}/{} cols effective/total)",
        orfs.len(), scaled_min_aa, effective_gap, total_cols
    ));

    let strand = &gff_entry.strand;

    for (aa_bytes, nt_start, nt_end, frame) in &orfs {
        let aa_str = String::from_utf8_lossy(aa_bytes).to_string();
        let nt_slice = &seq[*nt_start..*nt_end];
        let nt_str = String::from_utf8_lossy(nt_slice).to_string();

        // Compute genomic coordinates
        let (hit_start, hit_end) = if strand == "+" {
            (rs + nt_start + 1, rs + nt_end)
        } else {
            let slen = seq.len();
            (rs + slen - nt_end + 1, rs + slen - nt_start)
        };

        let mut frame_val = (*frame as i32) + 1;
        if node.frame < 0 {
            frame_val = -frame_val;
        }

        log.push(format!(
            "  ORF: frame={} aa_len={} genomic={}:{}-{}({})\n  aa={}\n  nt={}",
            frame_val,
            aa_str.len(),
            gff_entry.scaffold,
            hit_start,
            hit_end,
            strand,
            aa_str,
            nt_str
        ));

        results.push(FlankCandidate {
            aa_seq: aa_str,
            nt_seq: nt_str,
            frame: frame_val,
            scaffold: gff_entry.scaffold.clone(),
            hit_start,
            hit_end,
            strand: strand.clone(),
            node_name: node_field.to_string(),
            is_leading,
            gap_start,
            gap_end,
            cluster_key: cluster_key.to_string(),
        });
    }

    log.push(String::new());
}

/// Find the start and end of non-gap content in a sequence.
fn flank_find_index_pair(seq: &str) -> (usize, usize) {
    let bytes = seq.as_bytes();
    let mut start = 0;
    let mut end = bytes.len();
    for (i, &b) in bytes.iter().enumerate() {
        if b != b'-' {
            start = i;
            break;
        }
    }
    for i in (0..bytes.len()).rev() {
        if bytes[i] != b'-' {
            end = i + 1;
            break;
        }
    }
    (start, end)
}

/// Run flank scan on a single gene.
fn flank_scan_gene(
    gene: &str,
    aa_entries: &[(String, String)],
    clusters_map: &HashMap<String, HashMap<String, (Vec<String>, String)>>,
    gff: &HashMap<String, GffEntry>,
    genome: &HashMap<String, Vec<u8>>,
    scaffold_intervals: &HashMap<String, Vec<(usize, usize, String)>>,
) -> (Vec<String>, Vec<FlankCandidate>) {
    let mut log: Vec<String> = Vec::new();
    let mut all_candidates: Vec<FlankCandidate> = Vec::new();

    log.push(format!("=== Flank ORF scan: {} ===", gene));

    // Build ref consensus and node list from AA entries
    let mut ref_consensus: HashMap<usize, Vec<u8>> = HashMap::new();
    let mut ref_starts: Vec<usize> = Vec::new();
    let mut ref_ends: Vec<usize> = Vec::new();
    let mut ref_count: usize = 0;
    let mut aa_nodes: Vec<FlankNode> = Vec::new();

    for (header, seq) in aa_entries {
        if header.ends_with('.') {
            ref_count += 1;
            let (start, end) = flank_find_index_pair(seq);
            ref_starts.push(start);
            ref_ends.push(end);
            for (i, b) in seq.bytes().enumerate() {
                let col = ref_consensus.entry(i).or_default();
                if i >= start && i <= end {
                    col.push(b);
                } else {
                    col.push(b'-');
                }
            }
            continue;
        }
        let fields: Vec<&str> = header.split('|').collect();
        let frame: i32 = fields.get(4).and_then(|s| s.parse().ok()).unwrap_or(1);
        let (start, end) = flank_find_index_pair(seq);
        aa_nodes.push(FlankNode {
            header: header.clone(),
            frame,
            start,
            end,
        });
    }

    if aa_nodes.is_empty() || ref_count == 0 || ref_starts.is_empty() {
        return (log, all_candidates);
    }

    // Ref median start/end
    let median_idx = ref_starts.len() / 2;
    let ref_median_start = *ref_starts.select_nth_unstable(median_idx).1;
    let ref_median_end = *ref_ends.select_nth_unstable(median_idx).1;

    // Get cluster sets for this gene
    let Some(clusters) = clusters_map.get(gene) else {
        return (log, all_candidates);
    };

    // Sort by cluster key for deterministic iteration order: `clusters`
    // is a HashMap whose iteration order is hash-randomized per process.
    // Without this sort, `all_candidates` ends up in a run-varying order
    // and downstream Python tie-breaking produces different survivors.
    let mut cluster_sets: Vec<(&str, HashSet<String>, &str)> = clusters
        .iter()
        .map(|(key, (nodes, iso_type))| {
            (
                key.as_str(),
                nodes.iter().cloned().collect(),
                iso_type.as_str(),
            )
        })
        .collect();
    cluster_sets.sort_unstable_by_key(|(ck, _, _)| *ck);

    for (cluster_key, cluster_set, iso_type) in &cluster_sets {
        // Only search flanks for base clusters, not isoforms (e.g. "6" not "6_1", "6_2").
        // Isoforms share edge nodes with the base cluster so the flank search is identical.
        if cluster_key.contains('_') {
            continue;
        }
        let mut aa_subset: Vec<&FlankNode> = aa_nodes
            .iter()
            .filter(|n| {
                let nf = n.header.split('|').nth(3).unwrap_or("");
                cluster_set.contains(nf)
            })
            .collect();

        if aa_subset.is_empty() {
            continue;
        }

        if aa_subset[0].frame < 0 {
            aa_subset.sort_by(|a, b| b.start.cmp(&a.start));
        } else {
            aa_subset.sort_by(|a, b| a.start.cmp(&b.start));
        }

        let (first_node, last_node) = if aa_subset[0].frame < 0 {
            (
                aa_subset[aa_subset.len() - 1],
                aa_subset[0],
            )
        } else {
            (aa_subset[0], aa_subset[aa_subset.len() - 1])
        };

        // Log window decisions for this cluster
        log.push(format!(
            "Cluster {}: {} nodes, ref_median={}-{}, first_node={}(cols {}-{}), last_node={}(cols {}-{})",
            cluster_key,
            aa_subset.len(),
            ref_median_start,
            ref_median_end,
            parse_node_field(&first_node.header),
            first_node.start,
            first_node.end,
            parse_node_field(&last_node.header),
            last_node.start,
            last_node.end,
        ));
        for node in &aa_subset {
            log.push(format!(
                "  node: {} cols {}-{} frame={}",
                parse_node_field(&node.header),
                node.start,
                node.end,
                node.frame,
            ));
        }

        // Leading gap: from ref_median_start to first_node.start
        let has_leading = first_node.start > ref_median_start;
        let has_trailing = ref_median_end > last_node.end;
        log.push(format!(
            "  Leading gap: {} (ref_median_start={} < first_node.start={}? {})",
            if has_leading { format!("cols {}-{} ({} cols)", ref_median_start, first_node.start, first_node.start - ref_median_start) } else { "NONE".to_string() },
            ref_median_start,
            first_node.start,
            has_leading,
        ));
        log.push(format!(
            "  Trailing gap: {} (last_node.end={} < ref_median_end={}? {})",
            if has_trailing { format!("cols {}-{} ({} cols)", last_node.end, ref_median_end, ref_median_end - last_node.end) } else { "NONE".to_string() },
            last_node.end,
            ref_median_end,
            has_trailing,
        ));

        // N-terminal isoforms: skip leading (left) flank scan.
        // C-terminal isoforms: skip trailing (right) flank scan.
        // Base clusters and other types: scan both.
        let skip_leading = *iso_type == "N-terminal";
        let skip_trailing = *iso_type == "C-terminal";

        if has_leading && !skip_leading {
            flank_extract_orfs(
                ref_median_start,
                first_node.start,
                first_node,
                true,
                &mut log,
                &ref_consensus,
                gff,
                genome,
                scaffold_intervals,
                cluster_key,
                &mut all_candidates,
            );
        }

        // Trailing gap: from last_node.end to ref_median_end
        if has_trailing && !skip_trailing {
            flank_extract_orfs(
                last_node.end,
                ref_median_end,
                last_node,
                false,
                &mut log,
                &ref_consensus,
                gff,
                genome,
                scaffold_intervals,
                cluster_key,
                &mut all_candidates,
            );
        }
    }

    (log, all_candidates)
}

// ---------------------------------------------------------------------------
// End flank scan
// ---------------------------------------------------------------------------

struct GffEntry {
    scaffold: String,
    start: usize,
    end: usize,
    strand: String,
}

fn collect_needed_scaffolds(
    gene_clusters: &HashMap<String, HashMap<String, (Vec<String>, String)>>,
    gff: &HashMap<String, GffEntry>,
) -> HashSet<String> {
    let mut needed = HashSet::new();
    for (_gene, clusters) in gene_clusters {
        for (_ck, (node_tokens, _iso_type)) in clusters {
            for token in node_tokens {
                if let Some(entry) = gff.get(token.as_str()) {
                    needed.insert(entry.scaffold.clone());
                }
            }
        }
    }
    needed
}

fn load_genome_from_rocksdb(
    rocksdb_path: &str,
    needed: &HashSet<String>,
) -> HashMap<String, Vec<u8>> {
    let mut opts = Options::default();
    opts.set_error_if_exists(false);
    let db = match DB::open_for_read_only(&opts, rocksdb_path, false) {
        Ok(db) => db,
        Err(e) => {
            eprintln!("Failed to open RocksDB at {}: {}", rocksdb_path, e);
            return HashMap::new();
        }
    };

    // Try indexed lookups first (new format: scaffold_index key)
    if let Ok(Some(idx_raw)) = db.get(b"scaffold_index") {
        let idx_str = String::from_utf8_lossy(&idx_raw);
        // Parse index: each line is "name\tbatch_index\tbyte_offset\tseq_length"
        // Group needed scaffolds by batch so we only load each batch once.
        let mut batch_requests: HashMap<String, Vec<(String, usize, usize)>> = HashMap::new();
        for line in idx_str.lines() {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 4 {
                continue;
            }
            let name = parts[0];
            if !needed.contains(name) {
                continue;
            }
            let batch_id = parts[1].to_string();
            let offset: usize = parts[2].parse().unwrap_or(0);
            let length: usize = parts[3].parse().unwrap_or(0);
            batch_requests
                .entry(batch_id)
                .or_default()
                .push((name.to_string(), offset, length));
        }

        let mut scaffolds = HashMap::with_capacity(needed.len());
        for (batch_id, entries) in &batch_requests {
            let key = format!("parentbatch:{}", batch_id);
            let data = match db.get(key.as_bytes()) {
                Ok(Some(v)) => v,
                _ => continue,
            };
            for (name, offset, length) in entries {
                let end = (*offset + *length).min(data.len());
                let start = (*offset).min(end);
                scaffolds.insert(name.clone(), data[start..end].to_vec());
            }
        }
        return scaffolds;
    }

    // Fallback: legacy batched format (scan all batches)
    let recipe = match db.get(b"getall:parentbatches") {
        Ok(Some(v)) => String::from_utf8_lossy(&v).to_string(),
        _ => {
            eprintln!("No parentbatches recipe in RocksDB");
            return HashMap::new();
        }
    };

    let batch_keys: Vec<&str> = recipe.split(',').collect();
    let mut scaffolds = HashMap::new();

    for bk in &batch_keys {
        let key = format!("parentbatch:{}", bk);
        let data = match db.get(key.as_bytes()) {
            Ok(Some(v)) => String::from_utf8_lossy(&v).to_string(),
            _ => continue,
        };
        let mut pos = 0;
        let bytes = data.as_bytes();
        let dlen = bytes.len();
        while pos < dlen {
            let nl = match data[pos..].find('\n') {
                Some(p) => pos + p,
                None => break,
            };
            let name = &data[pos + 1..nl];
            pos = nl + 1;
            let nl2 = match data[pos..].find('\n') {
                Some(p) => pos + p,
                None => {
                    if needed.contains(name) {
                        scaffolds.insert(name.to_string(), data[pos..].as_bytes().to_vec());
                    }
                    break;
                }
            };
            if needed.contains(name) {
                scaffolds.insert(name.to_string(), data[pos..nl2].as_bytes().to_vec());
            }
            pos = nl2 + 1;
        }
    }
    scaffolds
}

fn read_gff(path: &str) -> HashMap<String, GffEntry> {
    let mut nodes = HashMap::new();
    let data = match fs::read_to_string(path) {
        Ok(data) => data,
        Err(_) => return nodes,
    };
    for line in data.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 {
            continue;
        }
        let scaffold = parts[0].to_string();
        let start: usize = parts[3].parse().unwrap_or(0);
        let end: usize = parts[4].parse().unwrap_or(0);
        let strand = parts[6].to_string();
        let attrs = parts[8];
        let mut name = None;
        for attr in attrs.split(';') {
            if let Some(v) = attr.strip_prefix("Name=") {
                name = Some(v.to_string());
                break;
            }
        }
        if let Some(n) = name {
            nodes.insert(n, GffEntry { scaffold, start, end, strand });
        }
    }
    nodes
}

fn read_clusters(path: &str) -> HashMap<String, HashMap<String, (Vec<String>, String)>> {
    let mut gene_clusters: HashMap<String, HashMap<String, (Vec<String>, String)>> = HashMap::new();
    let data = match fs::read_to_string(path) {
        Ok(data) => data,
        Err(_) => return gene_clusters,
    };
    let mut first = true;
    for line in data.lines() {
        if first {
            first = false;
            continue;
        }
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 7 {
            continue;
        }
        let gene = fields[0].trim_end_matches(".gz").to_string();
        let cluster_key = fields[1].to_string();
        let nodes: Vec<String> = fields[6]
            .split(';')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect();
        let isoform_type = if fields.len() > 8 {
            fields[8].trim().to_string()
        } else {
            String::new()
        };
        gene_clusters
            .entry(gene)
            .or_default()
            .insert(cluster_key, (nodes, isoform_type));
    }
    gene_clusters
}

fn read_aa_fasta(path: &str) -> Vec<(String, String)> {
    let mut entries = Vec::new();
    let data = if path.ends_with(".gz") {
        let file = match fs::File::open(path) {
            Ok(f) => f,
            Err(_) => return entries,
        };
        let mut decoder = GzDecoder::new(file);
        let mut s = String::new();
        if decoder.read_to_string(&mut s).is_err() {
            return entries;
        }
        s
    } else {
        match fs::read_to_string(path) {
            Ok(data) => data,
            Err(_) => return entries,
        }
    };
    let mut header = String::new();
    let mut seq = String::new();
    for line in data.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if !header.is_empty() {
                entries.push((header.clone(), seq.clone()));
            }
            header = rest.to_string();
            seq.clear();
        } else {
            seq.push_str(line.trim());
        }
    }
    if !header.is_empty() {
        entries.push((header, seq));
    }
    entries
}

struct RecoveredExon {
    gene: String,
    header: String,
    aa_seq: String,
    nt_seq: String,
    region: String,
    score: f64,
    cluster_key: String,
    node_a_name: String,
    node_b_name: String,
}

struct GapCandidate {
    aa_seq: String,
    nt_seq: String,
    frame: i32,
    scaffold: String,
    hit_start: usize,
    hit_end: usize,
    strand: String,
    gap_start: usize,
    gap_end: usize,
    cluster_key: String,
    node_a_name: String,
    node_b_name: String,
    /// For MXE region scan candidates: the base cluster's modular node
    /// at this slot.  Empty for normal gap candidates.
    mxe_slot_node: String,
}

struct GeneResult {
    log: Vec<String>,
    gaps: usize,
    hits: usize,
    gff_lines: Vec<String>,
    recovered: Vec<RecoveredExon>,
    gap_candidates: Vec<GapCandidate>,
}

fn process_gene(
    aa_file: &str,
    entries: &[(String, String)],
    clusters_map: &HashMap<String, HashMap<String, (Vec<String>, String)>>,
    gff: &HashMap<String, GffEntry>,
    genome: &HashMap<String, Vec<u8>>,
) -> GeneResult {
    let gene = aa_file.trim_end_matches(".gz");
    let mut flog: Vec<String> = Vec::new();
    let gff_lines: Vec<String> = Vec::new();
    let recovered: Vec<RecoveredExon> = Vec::new();
    let mut gap_candidates: Vec<GapCandidate> = Vec::new();
    let (mut tgaps, thits) = (0usize, 0usize);

    let Some(clusters) = clusters_map.get(gene) else {
        return GeneResult {
            log: flog, gaps: 0, hits: 0,
            gff_lines, recovered, gap_candidates,
        };
    };

    #[allow(dead_code)]
    struct Cand {
        header: String,
        seq: Vec<u8>,
        frame: i32,
        start: usize,
        end: usize,
        nf: String,
    }

    let mut rc: HashMap<usize, Vec<u8>> = HashMap::new();
    let mut rcount = 0usize;
    let mut cands: Vec<Cand> = Vec::new();

    for (header, seq) in entries {
        let sb = seq.as_bytes();
        if header.ends_with('.') {
            rcount += 1;
            let (s, e) = find_index_pair(sb, b'-');
            for (i, &bp) in sb.iter().enumerate() {
                let ent = rc.entry(i).or_default();
                ent.push(if i >= s && i <= e { bp } else { b'-' });
            }
        } else {
            let parts: Vec<&str> = header.split('|').collect();
            if parts.len() > 4 {
                if let Ok(frame) = parts[4].parse::<i32>() {
                    let (s, e) = find_index_pair(sb, b'-');
                    cands.push(Cand {
                        header: header.clone(),
                        seq: sb.to_vec(),
                        frame,
                        start: s,
                        end: e,
                        nf: parse_node_field(header).to_string(),
                    });
                }
            }
        }
    }

    if rcount == 0 || cands.is_empty() {
        return GeneResult {
            log: flog, gaps: 0, hits: 0,
            gff_lines, recovered, gap_candidates,
        };
    }

    flog.push(format!(
        "=== {}: {} refs, {} candidates ===",
        gene, rcount, cands.len()
    ));

    // Separate clusters into base and MXE, skip N-terminal
    let mut base_clusters: Vec<(&String, &Vec<String>, &String)> = Vec::new();
    let mut mxe_clusters: Vec<(&String, &Vec<String>, &String)> = Vec::new();
    for (ck, (node_tokens, iso_type)) in clusters.iter() {
        if iso_type == "N-terminal" {
            continue;
        }
        if iso_type == "MXE" {
            mxe_clusters.push((ck, node_tokens, iso_type));
        } else {
            base_clusters.push((ck, node_tokens, iso_type));
        }
    }
    base_clusters.sort_by_key(|(ck, _, _)| *ck);
    mxe_clusters.sort_by_key(|(ck, _, _)| *ck);

    let all_clusters: Vec<_> = base_clusters
        .iter()
        .chain(mxe_clusters.iter())
        .cloned()
        .collect();

    // Build a map of known GFF intervals for each MXE cluster group.
    // Key = base cluster key (e.g. "86"), Value = vec of (scaffold, start, end)
    // for every node across all sibling isoforms (including the base).
    // Used to narrow MXE gap search regions so the DP does not scan
    // genomic intervals already occupied by known cassette variants
    // from other isoforms.
    let mxe_sibling_gff: HashMap<String, Vec<(String, usize, usize)>> = {
        // Collect all MXE cluster keys grouped by base key
        let mut groups: HashMap<String, Vec<&Vec<String>>> = HashMap::new();
        for (ck, node_tokens, iso_type) in &all_clusters {
            if *iso_type != "MXE" {
                continue;
            }
            let base_key = ck
                .rsplit_once('_')
                .map(|(k, _)| k.to_string())
                .unwrap_or_else(|| (*ck).clone());
            groups.entry(base_key).or_default().push(node_tokens);
        }
        let mut result: HashMap<String, Vec<(String, usize, usize)>> = HashMap::new();
        for (base_key, token_lists) in &groups {
            let mut intervals = Vec::new();
            for tokens in token_lists {
                for token in *tokens {
                    if let Some(entry) = gff.get(token.as_str()) {
                        intervals.push((
                            entry.scaffold.clone(),
                            entry.start,
                            entry.end,
                        ));
                    }
                }
            }
            // Also include nodes from the base cluster itself, which
            // may not be typed MXE (e.g. "C-terminal" with MXE children).
            if let Some((base_tokens, _)) = clusters.get(base_key.as_str()) {
                for token in base_tokens {
                    if let Some(entry) = gff.get(token.as_str()) {
                        intervals.push((
                            entry.scaffold.clone(),
                            entry.start,
                            entry.end,
                        ));
                    }
                }
            }
            result.insert(base_key.clone(), intervals);
        }
        result
    };

    // For each MXE base cluster key, compute which of its nodes are
    // modular (i.e. get swapped out in at least one isoform).
    // A base node is modular if any isoform child does NOT contain it.
    // Note: the base cluster itself may not be typed "MXE" (e.g. it
    // could be "C-terminal" if the cluster has both MXE and C-terminal
    // splicing).  We identify base clusters by checking whether any
    // child key (with '_') is typed MXE.
    // BTreeMap: sorted iteration order is free and stable across runs.
    // Downstream loops push to `gap_candidates` in this order; a HashMap
    // here would randomize that and break downstream tie-breaking.
    let mxe_base_modular: BTreeMap<String, HashSet<String>> = {
        let mut result: BTreeMap<String, HashSet<String>> = BTreeMap::new();
        // Group MXE isoform node-sets by base key
        let mut iso_sets: HashMap<String, Vec<HashSet<String>>> = HashMap::new();
        for (ck, (node_tokens, iso_type)) in clusters.iter() {
            if iso_type != "MXE" {
                continue;
            }
            let fields: HashSet<String> = node_tokens.iter().cloned().collect();
            if let Some(bk) = ck.rsplit_once('_').map(|(k, _)| k.to_string()) {
                iso_sets.entry(bk).or_default().push(fields);
            }
        }
        // For each base key that has MXE children, look up the base
        // cluster's node list (regardless of its own iso_type) and
        // find which nodes are swapped out in at least one child.
        for (bk, children) in &iso_sets {
            let base_tokens = match clusters.get(bk.as_str()) {
                Some((tokens, _)) => tokens,
                None => continue,
            };
            let base_fields: HashSet<String> = base_tokens.iter().cloned().collect();
            let modular: HashSet<String> = base_fields
                .iter()
                .filter(|node| children.iter().any(|child| !child.contains(*node)))
                .cloned()
                .collect();
            result.insert(bk.clone(), modular);
        }
        result
    };

    for (ck, node_tokens, iso_type) in &all_clusters {
        let cluster_node_fields: HashSet<String> = node_tokens.iter().cloned().collect();

        // For MXE clusters, identify modular nodes.
        // Isoform clusters: nodes present in the isoform but not the base.
        // Base clusters: nodes in the base that are replaced in at
        //   least one isoform child.
        let modular_node_fields: HashSet<String> = if *iso_type == "MXE" {
            if let Some(base_key) = ck.rsplit_once('_').map(|(k, _)| k.to_string()) {
                if let Some((base_tokens, _)) = clusters.get(&base_key) {
                    let base_fields: HashSet<String> = base_tokens.iter().cloned().collect();
                    cluster_node_fields
                        .difference(&base_fields)
                        .cloned()
                        .collect()
                } else {
                    HashSet::new()
                }
            } else {
                // Base cluster: use precomputed modular set
                mxe_base_modular
                    .get(*ck)
                    .cloned()
                    .unwrap_or_default()
            }
        } else {
            HashSet::new()
        };

        // Filter candidates by cluster membership
        let mut aa_subset: Vec<usize> = cands
            .iter()
            .enumerate()
            .filter(|(_, c)| cluster_node_fields.contains(&c.nf))
            .map(|(i, _)| i)
            .collect();

        if aa_subset.len() < 2 {
            continue;
        }

        let is_negative = cands[aa_subset[0]].frame < 0;
        if is_negative {
            aa_subset.sort_by(|&a, &b| cands[b].start.cmp(&cands[a].start));
        } else {
            aa_subset.sort_by(|&a, &b| cands[a].start.cmp(&cands[b].start));
        }

        for wi in 1..aa_subset.len() {
            let na_idx = aa_subset[wi - 1];
            let nb_idx = aa_subset[wi];
            let na = &cands[na_idx];
            let nb = &cands[nb_idx];

            let (left_end, right_start) = if is_negative {
                (nb.end, na.start)
            } else {
                (na.end, nb.start)
            };

            let gap_aa = if right_start > left_end {
                right_start - left_end
            } else {
                0
            };
            if gap_aa == 0 {
                continue;
            }

            let gap_start = left_end;
            let gap_end = right_start;

            let Some((effective_gap, gap_min_aa)) = qualify_gap(gap_start, gap_end, &rc) else {
                continue;
            };
            let total_gap_cols = gap_end - gap_start;

            tgaps += 1;

            let node_a_name = na.nf.clone();
            let node_b_name = nb.nf.clone();

            // MXE isoforms: skip gaps between two shared (non-modular)
            // nodes.  These gaps duplicate work done by the base
            // cluster.  Only gaps touching a modular node may contain
            // undiscovered cassette variants in the flanking intron.
            // The base cluster itself must process shared-shared gaps
            // because no other cluster will.
            if *iso_type == "MXE" && ck.contains('_') && !modular_node_fields.is_empty() {
                let a_mod = modular_node_fields.contains(&na.nf);
                let b_mod = modular_node_fields.contains(&nb.nf);
                if !a_mod && !b_mod {
                    continue;
                }
            }

            let Some(ga) = gff.get(&node_a_name) else { continue; };
            let Some(gb) = gff.get(&node_b_name) else { continue; };

            if ga.scaffold != gb.scaffold {
                continue;
            }
            let strand = &ga.strand;

            let (region_start, region_end) = if ga.start < gb.start {
                (ga.end + 1, if gb.start > 0 { gb.start - 1 } else { 0 })
            } else {
                (gb.end + 1, if ga.start > 0 { ga.start - 1 } else { 0 })
            };

            // MXE gap narrowing.  Two complementary mechanisms:
            //
            // (A) All MXE clusters: narrow around pre-existing GFF
            //     nodes from sibling isoforms so the DP does not
            //     scan over known cassette exons and produce
            //     overlapping hits.
            //
            // (B) Isoform clusters only: further narrow around
            //     DP-recovered exons from the base cluster (which
            //     runs first) so each isoform searches only the
            //     sub-region belonging to its own cassette module.
            let (mut rs, mut re) = (region_start, region_end);
            if *iso_type == "MXE" {
                let base_key = ck
                    .rsplit_once('_')
                    .map(|(k, _)| k.to_string())
                    .unwrap_or_else(|| (*ck).clone());

                // (A) Narrow around sibling GFF nodes
                if let Some(sibling_intervals) = mxe_sibling_gff.get(&base_key) {
                    let foreign: Vec<(usize, usize)> = sibling_intervals
                        .iter()
                        .filter(|(sscaf, ss, se)| {
                            *sscaf == ga.scaffold
                                && *ss >= rs
                                && *se <= re
                                && !cluster_node_fields.iter().any(|nf| {
                                    if let Some(entry) = gff.get(nf.as_str()) {
                                        entry.scaffold == *sscaf
                                            && entry.start == *ss
                                            && entry.end == *se
                                    } else {
                                        false
                                    }
                                })
                        })
                        .map(|(_, ss, se)| (*ss, *se))
                        .collect();
                    if !foreign.is_empty() {
                        let a_mod = modular_node_fields.contains(&na.nf);
                        let b_mod = modular_node_fields.contains(&nb.nf);
                        let leftmost_foreign = foreign.iter().map(|(s, _)| *s).min().unwrap();
                        let rightmost_foreign_end =
                            foreign.iter().map(|(_, e)| *e).max().unwrap();
                        let mod_is_left = if a_mod && !b_mod {
                            ga.start < gb.start
                        } else if b_mod && !a_mod {
                            gb.start < ga.start
                        } else {
                            let left_size = leftmost_foreign.saturating_sub(rs);
                            let right_size = re.saturating_sub(rightmost_foreign_end);
                            left_size >= right_size
                        };
                        if mod_is_left {
                            re = re.min(leftmost_foreign - 1);
                        } else {
                            rs = rs.max(rightmost_foreign_end + 1);
                        }
                        if rs >= re {
                            continue;
                        }
                    }
                }

            }

            let Some(scaffold_seq) = genome.get(&ga.scaffold) else { continue; };

            if re < rs || rs == 0 || re > scaffold_seq.len() {
                continue;
            }

            flog.push(format!(
                "  Gap {}-{}: {}:{}-{} ({})",
                node_a_name, node_b_name, ga.scaffold, rs, re, strand
            ));

            let sr: Vec<u8> = if strand == "+" {
                scaffold_seq[rs - 1..re].to_vec()
            } else {
                bio_revcomp(&scaffold_seq[rs - 1..re])
            };

            // Log the entire NT gap region in FASTA format
            flog.push(format!(
                "  >gap_region {}:{}-{}({}) len={}",
                ga.scaffold, rs, re, strand, sr.len()
            ));
            flog.push(format!("  {}", String::from_utf8_lossy(&sr)));

            if sr.len() < 6 {
                continue;
            }

            // Extract ORFs from the gap region.  No M requirement for
            // internal exons. min_aa scaled by effective gap (ref-presence cols).
            let orfs = extract_flank_orfs(&sr, gap_min_aa, false);

            flog.push(format!(
                "  {} candidate ORFs (min_aa={}, gap={}/{} eff/total cols)",
                orfs.len(), gap_min_aa, effective_gap, total_gap_cols
            ));

            if orfs.is_empty() {
                continue;
            }

            for (aa_bytes, nt_start, nt_end, frame) in &orfs {
                let aa_str = String::from_utf8_lossy(aa_bytes).to_string();
                let nt_slice = &sr[*nt_start..*nt_end];
                let nt_str = String::from_utf8_lossy(nt_slice).to_string();

                let (hss, hse) = if strand == "+" {
                    (rs + nt_start, rs + nt_end - 1)
                } else {
                    let rl = re - rs + 1;
                    (rs + (rl - nt_end), rs + (rl - nt_start) - 1)
                };

                let fv = if strand == "+" {
                    (*frame as i32) + 1
                } else {
                    -((*frame as i32) + 1)
                };

                flog.push(format!(
                    "    ORF: frame={} aa_len={} genomic={}:{}-{}({})\n      aa={}",
                    fv, aa_str.len(), ga.scaffold, hss, hse, strand, aa_str
                ));

                gap_candidates.push(GapCandidate {
                    aa_seq: aa_str,
                    nt_seq: nt_str,
                    frame: fv,
                    scaffold: ga.scaffold.clone(),
                    hit_start: hss,
                    hit_end: hse,
                    strand: strand.clone(),
                    gap_start,
                    gap_end,
                    cluster_key: ck.to_string(),
                    node_a_name: node_a_name.clone(),
                    node_b_name: node_b_name.clone(),
                    mxe_slot_node: String::new(),
                });
            }
        }
    }

    // MXE region scan: for each modular slot, simulate a gap between
    // the two constitutive flanking nodes and extract ORFs from the
    // full intervening genomic region.
    {
        let mut mxe_scanned: HashSet<(String, String, String)> = HashSet::new();

        for (base_key, modular_set) in &mxe_base_modular {
            let Some((base_tokens, _)) = clusters.get(base_key.as_str()) else {
                continue;
            };
            if modular_set.is_empty() {
                continue;
            }

            // Known sibling cassette intervals for this MXE group.
            // Used to trim (not reject) overlapping ORFs so novel cassettes
            // flanking or between known siblings can still surface.
            let known_intervals: Vec<(String, usize, usize)> = mxe_sibling_gff
                .get(base_key)
                .cloned()
                .unwrap_or_default();

            for (slot_idx, token) in base_tokens.iter().enumerate() {
                if !modular_set.contains(token) {
                    continue;
                }
                if slot_idx == 0 || slot_idx >= base_tokens.len() - 1 {
                    continue;
                }

                let left_node = &base_tokens[slot_idx - 1];
                let right_node = &base_tokens[slot_idx + 1];
                if modular_set.contains(left_node) || modular_set.contains(right_node) {
                    continue;
                }

                let Some(left_gff) = gff.get(left_node.as_str()) else { continue; };
                let Some(right_gff) = gff.get(right_node.as_str()) else { continue; };
                if left_gff.scaffold != right_gff.scaffold {
                    continue;
                }
                let scaffold = &left_gff.scaffold;
                let strand = &left_gff.strand;

                let scan_key = (scaffold.clone(), left_node.clone(), right_node.clone());
                if mxe_scanned.contains(&scan_key) {
                    continue;
                }
                mxe_scanned.insert(scan_key);

                // MSA gap window between the two constitutive flanks
                let left_cand = cands.iter().find(|c| c.nf == *left_node);
                let right_cand = cands.iter().find(|c| c.nf == *right_node);
                let (gap_start, gap_end) = match (left_cand, right_cand) {
                    (Some(lc), Some(rc)) => {
                        if strand == "-" { (rc.end, lc.start) } else { (lc.end, rc.start) }
                    }
                    _ => continue,
                };
                if gap_end <= gap_start {
                    continue;
                }

                let Some((effective_gap, gap_min_aa)) = qualify_gap(gap_start, gap_end, &rc) else {
                    continue;
                };

                // Genomic region between the two constitutive flanks
                let (rs, re) = if left_gff.start < right_gff.start {
                    (left_gff.end + 1, right_gff.start.saturating_sub(1))
                } else {
                    (right_gff.end + 1, left_gff.start.saturating_sub(1))
                };
                if rs >= re {
                    continue;
                }

                let Some(scaffold_seq) = genome.get(scaffold) else { continue; };
                if re > scaffold_seq.len() {
                    continue;
                }

                let total_gap_cols = gap_end - gap_start;
                flog.push(format!(
                    "  MXE_SCAN slot={} flanks={}..{} region={}:{}-{}({}) gap_cols={}-{} eff={}/{} min_aa={}",
                    token, left_node, right_node, scaffold,
                    rs, re, strand, gap_start, gap_end, effective_gap, total_gap_cols, gap_min_aa,
                ));

                let sr: Vec<u8> = if strand == "+" {
                    scaffold_seq[rs - 1..re].to_vec()
                } else {
                    bio_revcomp(&scaffold_seq[rs - 1..re])
                };
                if sr.len() < 6 {
                    continue;
                }

                flog.push(format!(
                    "  >mxe_region {}:{}-{}({}) len={}",
                    scaffold, rs, re, strand, sr.len()
                ));
                flog.push(format!("  {}", String::from_utf8_lossy(&sr)));

                let orfs = extract_flank_orfs(&sr, gap_min_aa, false);
                let mut kept = 0usize;

                flog.push(format!(
                    "  {} candidate ORFs (min_aa={}, gap={}/{} eff/total cols)",
                    orfs.len(), gap_min_aa, effective_gap, total_gap_cols
                ));

                for (aa_bytes, nt_start, nt_end, frame) in &orfs {
                    let (orf_hss, orf_hse) = if strand == "+" {
                        (rs + nt_start, rs + nt_end - 1)
                    } else {
                        let rl = re - rs + 1;
                        (rs + (rl - nt_end), rs + (rl - nt_start) - 1)
                    };

                    let fv = if strand == "+" {
                        (*frame as i32) + 1
                    } else {
                        -((*frame as i32) + 1)
                    };

                    let aa_len = aa_bytes.len();

                    // Project each overlapping known interval onto this ORF's
                    // AA index space, rounding inward (ceil on both sides) so
                    // fragments never contain partial-codon overlap with a known.
                    // For + strand the ORF's AA index i covers genomic
                    // [orf_hss + i*3, orf_hss + i*3 + 2]. For - strand the AA
                    // string runs 5' -> 3' on the reverse complement, so AA
                    // index i covers genomic [orf_hse - i*3 - 2, orf_hse - i*3],
                    // and knowns project with the orientation flipped.
                    let mut proj: Vec<(usize, usize)> = Vec::new();
                    for (ks_scaf, ks, ke) in &known_intervals {
                        if ks_scaf != scaffold {
                            continue;
                        }
                        let ks = *ks;
                        let ke = *ke;
                        if ke < orf_hss || ks > orf_hse {
                            continue;
                        }
                        let (aa_ks, aa_ke) = if strand == "+" {
                            // AA i on + strand covers genomic [orf_hss + 3i, orf_hss + 3i + 2].
                            // First overlapping AA: smallest i with (orf_hss + 3i + 2) >= ks
                            //   => i >= ceil((ks - orf_hss - 2) / 3)
                            // First non-overlapping AA past the known: smallest i with
                            // (orf_hss + 3i) > ke  =>  i >= ceil((ke + 1 - orf_hss) / 3)
                            let a = ((ks.saturating_sub(orf_hss + 2) + 2) / 3).min(aa_len);
                            let b = (((ke + 1).saturating_sub(orf_hss) + 2) / 3).min(aa_len);
                            (a, b)
                        } else {
                            // AA i on - strand covers genomic [orf_hse - 3i - 2, orf_hse - 3i].
                            // First overlapping AA: smallest i with (orf_hse - 3i - 2) <= ke
                            //   => i >= ceil((orf_hse - ke - 2) / 3)
                            // First non-overlapping AA past the known: smallest i with
                            // (orf_hse - 3i) < ks  =>  i >= ceil((orf_hse + 1 - ks) / 3)
                            let a = ((orf_hse.saturating_sub(ke + 2) + 2) / 3).min(aa_len);
                            let b = (((orf_hse + 1).saturating_sub(ks) + 2) / 3).min(aa_len);
                            (a, b)
                        };
                        if aa_ke > aa_ks {
                            proj.push((aa_ks, aa_ke));
                        }
                    }
                    proj.sort_by_key(|&(s, _)| s);

                    // Walk left-to-right, emitting fragments between knowns.
                    let mut fragments: Vec<(usize, usize)> = Vec::new();
                    let mut cursor = 0usize;
                    for (s, e) in &proj {
                        if *s > cursor {
                            fragments.push((cursor, *s));
                        }
                        if *e > cursor {
                            cursor = *e;
                        }
                    }
                    if cursor < aa_len {
                        fragments.push((cursor, aa_len));
                    }

                    for (aa_from, aa_to) in fragments {
                        let frag_len = aa_to - aa_from;
                        if frag_len < gap_min_aa {
                            continue;
                        }

                        let frag_aa = &aa_bytes[aa_from..aa_to];
                        let frag_nt_start = nt_start + aa_from * 3;
                        let frag_nt_end = nt_start + aa_to * 3;

                        let (hss, hse) = if strand == "+" {
                            (rs + frag_nt_start, rs + frag_nt_end - 1)
                        } else {
                            let rl = re - rs + 1;
                            (rs + (rl - frag_nt_end), rs + (rl - frag_nt_start) - 1)
                        };

                        let aa_str = String::from_utf8_lossy(frag_aa).to_string();
                        let nt_slice = &sr[frag_nt_start..frag_nt_end];
                        let nt_str = String::from_utf8_lossy(nt_slice).to_string();

                        flog.push(format!(
                            "    MXE_ORF: frame={} aa_len={} genomic={}:{}-{}({})\n      aa={}",
                            fv, frag_len, scaffold, hss, hse, strand, aa_str,
                        ));

                        gap_candidates.push(GapCandidate {
                            aa_seq: aa_str,
                            nt_seq: nt_str,
                            frame: fv,
                            scaffold: scaffold.clone(),
                            hit_start: hss,
                            hit_end: hse,
                            strand: strand.clone(),
                            gap_start,
                            gap_end,
                            cluster_key: base_key.clone(),
                            node_a_name: left_node.clone(),
                            node_b_name: right_node.clone(),
                            mxe_slot_node: token.clone(),
                        });
                        kept += 1;
                    }
                }

                flog.push(format!(
                    "  MXE_SCAN_RESULT: {} ORFs extracted, {} new candidates",
                    orfs.len(), kept,
                ));
            }
        }
    }

    GeneResult {
        log: flog,
        gaps: tgaps,
        hits: thits,
        gff_lines,
        recovered,
        gap_candidates,
    }
}

fn find_gff_path(folder: &Path, sub_dir: &str) -> Option<String> {
    let candidates = [
        folder.join("outlier").join(sub_dir),
        folder.join("outlier").join("resolve"),
        folder.join(sub_dir),
    ];
    for dir in &candidates {
        if dir.exists() {
            if let Ok(entries) = fs::read_dir(dir) {
                let mut gffs: Vec<String> = entries
                    .flatten()
                    .filter_map(|e| {
                        let name = e.file_name().to_string_lossy().to_string();
                        if name.ends_with("_coords.gff") {
                            Some(e.path().to_string_lossy().to_string())
                        } else {
                            None
                        }
                    })
                    .collect();
                gffs.sort();
                if let Some(path) = gffs.first() {
                    return Some(path.clone());
                }
            }
        }
    }
    None
}

fn find_input_folder(folder: &Path, sub_dir: &str) -> Option<std::path::PathBuf> {
    let p1 = folder.join("outlier").join(sub_dir);
    if p1.exists() {
        return Some(p1);
    }
    let p2 = folder.join(sub_dir);
    if p2.exists() {
        return Some(p2);
    }
    None
}

#[pyfunction]
#[pyo3(signature = (folder, sub_dir, taxa_path))]
pub fn exon_dp(
    py: Python<'_>,
    folder: String,
    sub_dir: String,
    taxa_path: String,
) -> PyResult<Py<PyAny>> {
    // `folder` is the per-orthoset root (<taxa>/<orthoset>); subdirs like
    // outlier/<sub_dir> live under it. `taxa_path` is the parent <taxa> folder
    // and holds the shared sequences RocksDB.
    let folder_path = Path::new(&folder);
    let taxa_path_p = Path::new(&taxa_path);

    let input_folder = find_input_folder(folder_path, &sub_dir).ok_or_else(|| {
        pyo3::exceptions::PyFileNotFoundError::new_err(format!(
            "Input folder not found for sub_dir={}",
            sub_dir
        ))
    })?;

    let aa_dir = input_folder.join("aa");

    let gff_path = find_gff_path(folder_path, &sub_dir);

    let mut gff_nodes = match &gff_path {
        Some(path) => read_gff(path),
        None => HashMap::new(),
    };

    let cluster_csv = input_folder.join("resolve_clusters.csv");
    let gene_clusters = if cluster_csv.exists() {
        read_clusters(&cluster_csv.to_string_lossy())
    } else {
        HashMap::new()
    };

    // Determine which scaffolds are actually needed before hitting RocksDB.
    // If no cluster nodes map to any GFF entry we have nothing to process.
    let needed_scaffolds = collect_needed_scaffolds(&gene_clusters, &gff_nodes);

    if needed_scaffolds.is_empty() {
        let output = PyDict::new(py);
        output.set_item("results", PyList::empty(py))?;
        output.set_item("gff_nodes", PyDict::new(py))?;
        output.set_item("genes", PyList::empty(py))?;
        output.set_item("genome_count", 0usize)?;
        return Ok(output.into());
    }

    let rocksdb_path = taxa_path_p.join("sequences").join("nt");
    let genome = load_genome_from_rocksdb(&rocksdb_path.to_string_lossy(), &needed_scaffolds);

    if genome.is_empty() {
        return Err(pyo3::exceptions::PyRuntimeError::new_err(format!(
            "No parent sequences found in RocksDB at {:?}",
            rocksdb_path
        )));
    }

    let mut gene_files: Vec<String> = Vec::new();
    if let Ok(entries) = fs::read_dir(&aa_dir) {
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().to_string();
            if name.ends_with(".fa") || name.contains(".fa.") {
                gene_files.push(name);
            }
        }
    }
    gene_files.sort();

    let mut results_list: Vec<GeneResult> = Vec::new();
    let mut output_genes: Vec<String> = Vec::new();
    let mut gene_entries_cache: Vec<Vec<(String, String)>> = Vec::new();

    for gene_file in &gene_files {
        let gene_key = gene_file.trim_end_matches(".gz");
        if !gene_clusters.contains_key(gene_key) {
            continue;
        }

        let fpath = aa_dir.join(gene_file);
        let entries = read_aa_fasta(&fpath.to_string_lossy());

        if entries.is_empty() {
            continue;
        }

        output_genes.push(gene_file.clone());

        results_list.push(process_gene(
            gene_key,
            &entries,
            &gene_clusters,
            &gff_nodes,
            &genome,
        ));
        gene_entries_cache.push(entries);
    }

    // -----------------------------------------------------------------------
    // Flank scan: update GFF with recovered exons, then scan flanks
    // -----------------------------------------------------------------------
    // Register recovered DP nodes into gff_nodes so flank scan can see them
    for result in &results_list {
        for rec in &result.recovered {
            let dp_name = parse_node_field(&rec.header).to_string();
            // Parse region string: "scaffold:start-end(strand)"
            if let Some(colon) = rec.region.find(':') {
                if let Some(paren) = rec.region.find('(') {
                    let scaffold = &rec.region[..colon];
                    let coords = &rec.region[colon + 1..paren];
                    let strand = &rec.region[paren + 1..paren + 2];
                    let parts: Vec<&str> = coords.split('-').collect();
                    if parts.len() == 2 {
                        if let (Ok(s), Ok(e)) = (parts[0].parse(), parts[1].parse()) {
                            gff_nodes.insert(dp_name, GffEntry {
                                scaffold: scaffold.to_string(),
                                start: s, end: e,
                                strand: strand.to_string(),
                            });
                        }
                    }
                }
            }
        }
    }

    let scaffold_intervals = build_scaffold_intervals(&gff_nodes);

    let mut flank_log_all: Vec<String> = Vec::new();
    let mut flank_candidates: Vec<(String, Vec<FlankCandidate>)> = Vec::new();

    for (gi, gene_file) in output_genes.iter().enumerate() {
        let gene_key = gene_file.trim_end_matches(".gz");
        let entries = &gene_entries_cache[gi];

        let (log, candidates) = flank_scan_gene(
            gene_key,
            entries,
            &gene_clusters,
            &gff_nodes,
            &genome,
            &scaffold_intervals,
        );

        if !log.is_empty() {
            flank_log_all.extend(log);
        }
        if !candidates.is_empty() {
            flank_candidates.push((gene_key.to_string(), candidates));
        }
    }

    // -----------------------------------------------------------------------
    // Build Python output
    // -----------------------------------------------------------------------
    let output = PyDict::new(py);

    let py_results = PyList::empty(py);
    for result in &results_list {
        let rd = PyDict::new(py);
        rd.set_item("log", &result.log)?;
        rd.set_item("gaps", result.gaps)?;
        rd.set_item("hits", result.hits)?;
        rd.set_item("gff_lines", &result.gff_lines)?;

        let py_recovered = PyList::empty(py);
        for rec in &result.recovered {
            let rec_dict = PyDict::new(py);
            rec_dict.set_item("gene", &rec.gene)?;
            rec_dict.set_item("header", &rec.header)?;
            rec_dict.set_item("aa_seq", &rec.aa_seq)?;
            rec_dict.set_item("nt_seq", &rec.nt_seq)?;
            rec_dict.set_item("region", &rec.region)?;
            rec_dict.set_item("score", rec.score)?;
            rec_dict.set_item("cluster_key", &rec.cluster_key)?;
            rec_dict.set_item("node_a_name", &rec.node_a_name)?;
            rec_dict.set_item("node_b_name", &rec.node_b_name)?;
            py_recovered.append(rec_dict)?;
        }
        rd.set_item("recovered", py_recovered)?;

        let py_gap_cands = PyList::empty(py);
        for gc in &result.gap_candidates {
            let gcd = PyDict::new(py);
            gcd.set_item("aa_seq", &gc.aa_seq)?;
            gcd.set_item("nt_seq", &gc.nt_seq)?;
            gcd.set_item("frame", gc.frame)?;
            gcd.set_item("scaffold", &gc.scaffold)?;
            gcd.set_item("hit_start", gc.hit_start)?;
            gcd.set_item("hit_end", gc.hit_end)?;
            gcd.set_item("strand", &gc.strand)?;
            gcd.set_item("gap_start", gc.gap_start)?;
            gcd.set_item("gap_end", gc.gap_end)?;
            gcd.set_item("cluster_key", &gc.cluster_key)?;
            gcd.set_item("node_a_name", &gc.node_a_name)?;
            gcd.set_item("node_b_name", &gc.node_b_name)?;
            if !gc.mxe_slot_node.is_empty() {
                gcd.set_item("mxe_slot_node", &gc.mxe_slot_node)?;
            }
            py_gap_cands.append(gcd)?;
        }
        rd.set_item("gap_candidates", py_gap_cands)?;

        py_results.append(rd)?;
    }
    output.set_item("results", py_results)?;

    // Flank candidate results
    let py_flank = PyList::empty(py);
    for (gene, candidates) in &flank_candidates {
        let gd = PyDict::new(py);
        gd.set_item("gene", gene)?;
        let py_cands = PyList::empty(py);
        for cand in candidates {
            let cd = PyDict::new(py);
            cd.set_item("aa_seq", &cand.aa_seq)?;
            cd.set_item("nt_seq", &cand.nt_seq)?;
            cd.set_item("frame", cand.frame)?;
            cd.set_item("scaffold", &cand.scaffold)?;
            cd.set_item("hit_start", cand.hit_start)?;
            cd.set_item("hit_end", cand.hit_end)?;
            cd.set_item("strand", &cand.strand)?;
            cd.set_item("node_name", &cand.node_name)?;
            cd.set_item("is_leading", cand.is_leading)?;
            cd.set_item("gap_start", cand.gap_start)?;
            cd.set_item("gap_end", cand.gap_end)?;
            cd.set_item("cluster_key", &cand.cluster_key)?;
            py_cands.append(cd)?;
        }
        gd.set_item("candidates", py_cands)?;
        py_flank.append(gd)?;
    }
    output.set_item("flank_candidates", py_flank)?;
    output.set_item("flank_log", flank_log_all)?;

    let py_gff = PyDict::new(py);
    for (name, entry) in &gff_nodes {
        let ed = PyDict::new(py);
        ed.set_item("scaffold", &entry.scaffold)?;
        ed.set_item("start", entry.start)?;
        ed.set_item("end", entry.end)?;
        ed.set_item("strand", &entry.strand)?;
        py_gff.set_item(name, ed)?;
    }
    output.set_item("gff_nodes", py_gff)?;
    output.set_item("genes", &output_genes)?;
    output.set_item("genome_count", genome.len())?;

    Ok(output.into())
}
