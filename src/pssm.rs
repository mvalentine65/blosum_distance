//! PSSM helpers for ExonFinder candidate scoring.
//!
//! `compute_pssm_for_window` is the Rust port of exonfinder's
//! `_compute_ref_data_cols` + `_precompute_pssm`.  Returns the
//! `(data_cols, pssm, full_max)` triple Python's window cache stores.
//!
//! Behaviour matches the Python helpers: a column enters `data_cols`
//! when `>= min_ref_occ` of the refs are non-gap/non-dot (`*` counts
//! toward the threshold but is filtered when tallying AA counts, since
//! it isn't in the BLOSUM62 standard set).  A column is dropped from the
//! PSSM if its max score is `<= 0`.

use pyo3::prelude::*;
use std::collections::{HashMap, HashSet};

/// The 20 standard AAs in the same order Python's `_BLOSUM62` enumerates them.
const STANDARD_AAS: [u8; 20] = *b"ARNDCQEGHILKMFPSTWYV";

/// Map an uppercase ASCII byte to the dense `STANDARD_AAS` index, or `None`
/// if the byte isn't one of the 20 standard AAs.
#[inline]
fn aa_index(c: u8) -> Option<usize> {
    match c {
        b'A' => Some(0),
        b'R' => Some(1),
        b'N' => Some(2),
        b'D' => Some(3),
        b'C' => Some(4),
        b'Q' => Some(5),
        b'E' => Some(6),
        b'G' => Some(7),
        b'H' => Some(8),
        b'I' => Some(9),
        b'L' => Some(10),
        b'K' => Some(11),
        b'M' => Some(12),
        b'F' => Some(13),
        b'P' => Some(14),
        b'S' => Some(15),
        b'T' => Some(16),
        b'W' => Some(17),
        b'Y' => Some(18),
        b'V' => Some(19),
        _ => None,
    }
}

/// Combined `_compute_ref_data_cols` + `_precompute_pssm` port.
///
/// `data_cols` holds alignment columns in `[win_start, win_end)` that pass
/// the `>= min_ref_occ` ref-occupancy gate.  `pssm` maps each kept column
/// to its BLOSUM62-weighted score per query AA plus the column max.
/// `full_max` is the sum of those column maxima.
#[pyfunction]
#[pyo3(signature = (ref_seqs, win_start, win_end, min_ref_occ = 0.30))]
pub fn compute_pssm_for_window(
    ref_seqs: Vec<String>,
    win_start: usize,
    win_end: usize,
    min_ref_occ: f64,
) -> PyResult<(
    HashSet<usize>,
    HashMap<usize, (HashMap<String, f64>, f64)>,
    f64,
)> {
    let n_refs = ref_seqs.len();
    if n_refs == 0 || win_end <= win_start {
        return Ok((HashSet::new(), HashMap::new(), 0.0));
    }
    let n_refs_f = n_refs as f64;

    // Pre-convert refs to byte slices once for fast indexing.
    let refs_bytes: Vec<&[u8]> = ref_seqs.iter().map(|s| s.as_bytes()).collect();

    // -----------------------------------------------------------------
    // _compute_ref_data_cols
    // -----------------------------------------------------------------
    let mut data_cols: HashSet<usize> = HashSet::new();
    for col in win_start..win_end {
        let mut count = 0usize;
        for r in &refs_bytes {
            if col < r.len() {
                let c = r[col];
                if c != b'-' && c != b'.' {
                    count += 1;
                }
            }
        }
        if (count as f64) / n_refs_f >= min_ref_occ {
            data_cols.insert(col);
        }
    }

    // -----------------------------------------------------------------
    // _precompute_pssm
    // -----------------------------------------------------------------
    let mut pssm: HashMap<usize, (HashMap<String, f64>, f64)> = HashMap::new();
    let mut full_max = 0.0_f64;

    for &col in &data_cols {
        // Tally per-AA counts on a flat 20-slot array (faster than HashMap).
        let mut counts = [0i32; 20];
        let mut total = 0i32;
        for r in &refs_bytes {
            if col >= r.len() {
                continue;
            }
            let c = r[col];
            if c == b'-' || c == b'.' {
                continue;
            }
            let cu = c.to_ascii_uppercase();
            if let Some(idx) = aa_index(cu) {
                counts[idx] += 1;
                total += 1;
            }
        }
        if total == 0 {
            continue;
        }
        let total_f = total as f64;

        // BLOSUM-weighted average per query AA across the ref distribution.
        let mut scores: HashMap<String, f64> = HashMap::with_capacity(20);
        let mut col_max = f64::NEG_INFINITY;
        for &q in &STANDARD_AAS {
            let mut s = 0.0_f64;
            for (i, &cnt) in counts.iter().enumerate() {
                if cnt > 0 {
                    s += (bio::scores::blosum62(q, STANDARD_AAS[i]) as f64) * (cnt as f64);
                }
            }
            let normalized = s / total_f;
            scores.insert((q as char).to_string(), normalized);
            if normalized > col_max {
                col_max = normalized;
            }
        }
        if col_max > 0.0 {
            pssm.insert(col, (scores, col_max));
            full_max += col_max;
        }
    }

    Ok((data_cols, pssm, full_max))
}
