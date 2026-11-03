use ahash::AHashSet;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

/// Finds the reads of prepare's `ntbatch` blobs that hold any of a set of seeds.
///
/// Flanks mode grows every edge of a sample by searching all its reads for
/// short seeds, a few rounds over hundreds of MB. A literal multi-pattern
/// search slows to ~55 MB/s once there are tens of thousands of DNA seeds; here
/// each read's k-mers are packed into integers and looked up, so the cost does
/// not depend on how many seeds there are. Case is ignored, as after upper().
#[pyclass]
pub struct SeedScanner {
    seeds: AHashSet<u64>,
    // one bit per hashed k-mer, checked before the set
    filter: Vec<u64>,
    k: usize,
}

const FILTER_BITS: u32 = 24;

#[inline]
fn code(byte: u8) -> u64 {
    match byte {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

#[inline]
fn slot(kmer: u64) -> usize {
    (kmer.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> (64 - FILTER_BITS)) as usize
}

#[inline]
fn find_newline(haystack: &[u8], from: usize) -> Option<usize> {
    memchr::memchr(b'\n', &haystack[from..]).map(|i| from + i)
}

impl SeedScanner {
    fn holds_seed(&self, seq: &[u8]) -> bool {
        let mask = if self.k == 32 { u64::MAX } else { (1u64 << (2 * self.k)) - 1 };
        let (mut kmer, mut run) = (0u64, 0usize);
        for &byte in seq {
            let c = code(byte);
            if c == 4 {
                run = 0;
                continue;
            }
            kmer = ((kmer << 2) | c) & mask;
            run += 1;
            if run >= self.k {
                let s = slot(kmer);
                if self.filter[s >> 6] >> (s & 63) & 1 == 1 && self.seeds.contains(&kmer) {
                    return true;
                }
            }
        }
        false
    }
}

#[pymethods]
impl SeedScanner {
    /// `seeds` are k-mers of A, C, G and T, all of length `k` (at most 32).
    /// Give both strands if both are wanted.
    #[new]
    fn new(seeds: Vec<String>, k: usize) -> PyResult<Self> {
        if k == 0 || k > 32 {
            return Err(PyValueError::new_err("k must be between 1 and 32"));
        }
        let mut set = AHashSet::with_capacity(seeds.len());
        let mut filter = vec![0u64; 1 << (FILTER_BITS - 6)];
        for seed in &seeds {
            let bytes = seed.as_bytes();
            if bytes.len() != k {
                return Err(PyValueError::new_err(format!("seed {seed:?} is not {k} long")));
            }
            let mut kmer = 0u64;
            for &byte in bytes {
                let c = code(byte);
                if c == 4 {
                    return Err(PyValueError::new_err(format!("seed {seed:?} holds a base other than ACGT")));
                }
                kmer = (kmer << 2) | c;
            }
            set.insert(kmer);
            let s = slot(kmer);
            filter[s >> 6] |= 1 << (s & 63);
        }
        Ok(SeedScanner { seeds: set, filter, k })
    }

    fn __len__(&self) -> usize {
        self.seeds.len()
    }

    /// Scan one blob and return, in order, the upper-cased sequence of every
    /// record that holds a seed.
    fn scan(&self, py: Python<'_>, batch: &[u8]) -> PyResult<Vec<Py<PyBytes>>> {
        let spans = py.detach(|| {
            let mut spans = Vec::new();
            let mut pos = 0usize;
            while pos < batch.len() {
                let header_end = match find_newline(batch, pos) {
                    Some(i) => i,
                    None => break,
                };
                if batch[pos] != b'>' {
                    return Err(pos);
                }
                pos = header_end + 1;
                let seq_end = find_newline(batch, pos).unwrap_or(batch.len());
                if self.holds_seed(&batch[pos..seq_end]) {
                    spans.push((pos, seq_end));
                }
                pos = seq_end + 1;
            }
            Ok(spans)
        });
        let spans = spans.map_err(|pos| {
            PyValueError::new_err(format!("ntbatch record at byte {pos} does not start with '>'"))
        })?;
        Ok(spans
            .into_iter()
            .map(|(a, b)| PyBytes::new(py, &batch[a..b].to_ascii_uppercase()).unbind())
            .collect())
    }
}
