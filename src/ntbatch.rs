use ahash::AHashSet;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict};

/// Pulls the sequences for a set of node ids out of prepare's `ntbatch` blobs.
///
/// A blob is concatenated FASTA (`>id\nSEQ\n`) holding half a million reads.
/// Diamond wants roughly 1% of them but has to walk every record to find the
/// boundaries, so on a large reads set the Python loop parses ~112M integer
/// headers to keep ~1.4M sequences. The set is built once and reused for every
/// blob.
#[pyclass]
pub struct NtBatchScanner {
    wanted: AHashSet<u64>,
}

/// Parse a decimal node id. Rejects anything `int()` would not have accepted
/// here, so a malformed batch fails loudly rather than silently skipping reads.
fn parse_id(field: &[u8]) -> Option<u64> {
    if field.is_empty() {
        return None;
    }
    let mut value: u64 = 0;
    for &byte in field {
        let digit = byte.wrapping_sub(b'0');
        if digit > 9 {
            return None;
        }
        value = value.checked_mul(10)?.checked_add(digit as u64)?;
    }
    Some(value)
}

#[inline]
fn find_newline(haystack: &[u8], from: usize) -> Option<usize> {
    haystack[from..].iter().position(|&b| b == b'\n').map(|i| from + i)
}

#[pymethods]
impl NtBatchScanner {
    #[new]
    fn new(nodes: Vec<u64>) -> Self {
        let mut wanted = AHashSet::with_capacity(nodes.len());
        wanted.extend(nodes);
        NtBatchScanner { wanted }
    }

    fn __len__(&self) -> usize {
        self.wanted.len()
    }

    /// Scan one blob, inserting `node id -> sequence` into `out` for every
    /// wanted record. Returns how many were inserted.
    ///
    /// Writes straight into the caller's dict so scanning every blob in a
    /// recipe needs no intermediate lists and no merge.
    fn scan_into(
        &self,
        py: Python<'_>,
        batch: &[u8],
        out: &Bound<'_, PyDict>,
    ) -> PyResult<usize> {
        let mut found = 0usize;
        let mut pos = 0usize;

        while pos < batch.len() {
            let header_end = match find_newline(batch, pos) {
                Some(i) => i,
                None => break,
            };
            if batch[pos] != b'>' {
                return Err(PyValueError::new_err(format!(
                    "ntbatch record at byte {pos} does not start with '>'"
                )));
            }
            let id = parse_id(&batch[pos + 1..header_end]).ok_or_else(|| {
                PyValueError::new_err(format!(
                    "ntbatch header at byte {pos} is not a node id"
                ))
            })?;

            pos = header_end + 1;
            let seq_end = match find_newline(batch, pos) {
                Some(i) => i,
                None => break,
            };

            if self.wanted.contains(&id) {
                out.set_item(id, PyBytes::new(py, &batch[pos..seq_end]))?;
                found += 1;
            }
            pos = seq_end + 1;
        }

        Ok(found)
    }
}
