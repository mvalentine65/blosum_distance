use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedStr;
use pyo3::types::PyDict;
use std::collections::HashSet;

#[pyfunction]
pub fn join_with_exclusions(string: &str, column_cull: HashSet<usize>) -> String {
    let chars = string.as_bytes();
    let mut result = Vec::with_capacity(string.len());
    for i in 0..string.len() {
        match column_cull.contains(&(i * 3)) {
            true => result.push(b'-'),
            false => result.push(chars[i]),
        }
    }
    String::from_utf8(result).unwrap()
}

/// Mask excluded codons with `---` across a gene's NT sequences.
///
/// A codon is excluded when its start position is in that sequence's
/// `exclusions` entry or in `shared_exclusion`, which is the same for every
/// sequence and so crosses from Python once. Sequences must be ASCII and a whole
/// number of codons.
#[pyfunction]
pub fn join_triplets_with_exclusions_many(
    sequences: Vec<PyBackedStr>,
    exclusions: Vec<HashSet<usize>>,
    shared_exclusion: HashSet<usize>,
) -> PyResult<Vec<String>> {
    if sequences.len() != exclusions.len() {
        return Err(PyValueError::new_err(
            "sequences and exclusions must be the same length",
        ));
    }
    sequences
        .iter()
        .zip(exclusions.iter())
        .map(|(sequence, exclusion)| {
            let bytes = sequence.as_bytes();
            if !sequence.is_ascii() || bytes.len() % 3 != 0 {
                return Err(PyValueError::new_err(
                    "sequences must be ASCII and a whole number of codons",
                ));
            }
            if exclusion.is_empty() && shared_exclusion.is_empty() {
                return Ok(sequence.to_string());
            }
            let mut result = Vec::with_capacity(bytes.len());
            for (i, codon) in bytes.chunks_exact(3).enumerate() {
                let position = i * 3;
                if exclusion.contains(&position) || shared_exclusion.contains(&position) {
                    result.extend_from_slice(b"---");
                } else {
                    result.extend_from_slice(codon);
                }
            }
            // Only ASCII bytes and '-' were copied.
            Ok(String::from_utf8(result).unwrap())
        })
        .collect()
}

type ByteSet = [u64; 4];

#[inline]
fn byte_set_has(set: &ByteSet, byte: u8) -> bool {
    set[(byte >> 6) as usize] & (1u64 << (byte & 63)) != 0
}

/// Bits for the single-character ASCII strings in a Python iterable. A longer
/// or non-ASCII string can never equal one ASCII sequence character, so it is
/// dropped without changing any membership answer.
fn byte_set_from(items: &Bound<'_, PyAny>) -> PyResult<ByteSet> {
    let mut set: ByteSet = [0; 4];
    for item in items.try_iter()? {
        let text: PyBackedStr = item?.extract()?;
        let bytes = text.as_bytes();
        if bytes.len() == 1 && bytes[0].is_ascii() {
            set[(bytes[0] >> 6) as usize] |= 1u64 << (bytes[0] & 63);
        }
    }
    Ok(set)
}

/// Per-column reference tables for flexcull's leading/trailing cull, built once
/// per gene.
#[pyclass]
pub struct CullTables {
    allowed: Vec<ByteSet>,
    blosum: Vec<ByteSet>,
    all_dashes: Vec<bool>,
    gap_present: Vec<bool>,
}

#[pymethods]
impl CullTables {
    /// Takes `process_refs`'s outputs. Columns run `0..len(all_dashes_by_index)`;
    /// a column missing from `all_dashes_by_index` or `gap_present_threshold` is
    /// a `ValueError`.
    #[new]
    fn new(
        character_at_each_pos: &Bound<'_, PyDict>,
        blosum_at_each_pos: &Bound<'_, PyDict>,
        all_dashes_by_index: &Bound<'_, PyDict>,
        gap_present_threshold: &Bound<'_, PyDict>,
    ) -> PyResult<Self> {
        let columns = all_dashes_by_index.len();
        let mut tables = CullTables {
            allowed: Vec::with_capacity(columns),
            blosum: Vec::with_capacity(columns),
            all_dashes: Vec::with_capacity(columns),
            gap_present: Vec::with_capacity(columns),
        };
        let missing = |name: &str, column: usize| {
            PyValueError::new_err(format!("{name} has no entry for column {column}"))
        };
        for column in 0..columns {
            tables.all_dashes.push(
                all_dashes_by_index
                    .get_item(column)?
                    .ok_or_else(|| missing("all_dashes_by_index", column))?
                    .extract()?,
            );
            tables.gap_present.push(
                gap_present_threshold
                    .get_item(column)?
                    .ok_or_else(|| missing("gap_present_threshold", column))?
                    .extract()?,
            );
            // Both are defaultdicts in Python: a missing column is empty.
            tables.allowed.push(match character_at_each_pos.get_item(column)? {
                Some(items) => byte_set_from(&items)?,
                None => [0; 4],
            });
            tables.blosum.push(match blosum_at_each_pos.get_item(column)? {
                Some(items) => byte_set_from(&items)?,
                None => [0; 4],
            });
        }
        Ok(tables)
    }

    fn __len__(&self) -> usize {
        self.all_dashes.len()
    }

    /// Leading/trailing cull for each sequence: `(cull_start, cull_end, kick)`.
    /// Sequences must be ASCII and no longer than the tables.
    fn cull_many(
        &self,
        sequences: Vec<PyBackedStr>,
        offset: i64,
        amt_matches: i64,
        mismatches: i64,
        blosum_max_percent: f64,
    ) -> PyResult<Vec<(Option<usize>, Option<usize>, bool)>> {
        sequences
            .iter()
            .map(|sequence| {
                if !sequence.is_ascii() || sequence.len() > self.all_dashes.len() {
                    return Err(PyValueError::new_err(format!(
                        "sequence must be ASCII and at most {} columns",
                        self.all_dashes.len()
                    )));
                }
                Ok(self.cull(
                    sequence.as_bytes(),
                    offset,
                    amt_matches,
                    mismatches,
                    blosum_max_percent,
                ))
            })
            .collect()
    }
}

impl CullTables {
    #[inline]
    fn window_passes(&self, sequence: &[u8], i: i64, step: i64, amt_matches: i64, mismatches: i64,
                     blosum_max_percent: f64) -> bool {
        let n = sequence.len() as i64;
        let first = sequence[i as usize];
        let mut mismatch = mismatches;
        let mut blosum_matches: i64 = byte_set_has(&self.blosum[i as usize], first) as i64;
        let mut checks = amt_matches - 1;
        let mut match_i: i64 = 1;

        while checks > 0 {
            let j = i + step * match_i;
            if step > 0 {
                if i + match_i + 1 >= n {
                    return false;
                }
            } else if j < 0 {
                return false;
            }
            let ju = j as usize;
            let byte = sequence[ju];
            if byte == b'-' {
                let next_gap = if step > 0 {
                    self.gap_present[ju + 1]
                } else {
                    j - 1 < 0 || self.gap_present[ju - 1]
                };
                if self.gap_present[ju] && next_gap {
                    return false;
                }
                match_i += 1;
            } else if !byte_set_has(&self.allowed[ju], byte) {
                mismatch -= 1;
                if mismatch < 0 || byte == b'*' {
                    return false;
                }
                match_i += 1;
                checks -= 1;
            } else {
                if byte_set_has(&self.blosum[ju], byte) {
                    blosum_matches += 1;
                }
                match_i += 1;
                checks -= 1;
            }
        }

        !(blosum_max_percent != -1.0 && blosum_matches as f64 / match_i as f64 > blosum_max_percent)
    }

    /// Scan forward for the first column that opens a passing match window,
    /// then backward for the last. A residue that does not match its column
    /// never opens a window. Kicks when the scans run into `offset`.
    fn cull(&self, sequence: &[u8], offset: i64, amt_matches: i64, mismatches: i64,
            blosum_max_percent: f64) -> (Option<usize>, Option<usize>, bool) {
        let n = sequence.len() as i64;
        let mut cull_start: Option<i64> = None;
        let mut cull_end: Option<i64> = None;
        let mut kick = false;

        for i in 0..n {
            if i == n - offset {
                kick = true;
                break;
            }
            let iu = i as usize;
            let byte = sequence[iu];
            if byte == b'-' || self.all_dashes[iu] || !byte_set_has(&self.allowed[iu], byte) {
                continue;
            }
            if mismatches < 0 {
                continue;
            }
            if self.window_passes(sequence, i, 1, amt_matches, mismatches, blosum_max_percent) {
                cull_start = Some(i);
                break;
            }
        }

        if let (false, Some(start)) = (kick, cull_start) {
            for i in (0..n).rev() {
                if i < start + offset {
                    kick = true;
                    break;
                }
                let iu = i as usize;
                let byte = sequence[iu];
                if byte == b'-' || self.all_dashes[iu] || !byte_set_has(&self.allowed[iu], byte) {
                    continue;
                }
                if mismatches < 0 {
                    continue;
                }
                if self.window_passes(sequence, i, -1, amt_matches, mismatches, blosum_max_percent) {
                    cull_end = Some(i + 1);
                    break;
                }
            }
        }

        (cull_start.map(|v| v as usize), cull_end.map(|v| v as usize), kick)
    }
}
