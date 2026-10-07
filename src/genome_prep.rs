//! Genome FASTA for prepare: scaffolds read a block at a time, their residues
//! checked and cut into chunks, their N-free runs followed, their windows laid out.

use memchr::memchr;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyList};

/// residues shown either side of an illegal one
const PREVIEW: usize = 20;

const OK: u8 = 0;
const HAS_N: u8 = 0x40;
const BAD: u8 = 0x80;

/// OK for a IUPAC nucleotide code in either case, HAS_N for `N`, BAD otherwise.
const CLASS: [u8; 256] = {
    let mut t = [BAD; 256];
    let valid = b"ACGTUNRYSWKMBDHVacgtunryswkmbdhv";
    let mut i = 0;
    while i < valid.len() {
        t[valid[i] as usize] = OK;
        i += 1;
    }
    t[b'N' as usize] = HAS_N;
    t
};

/// The classes of a line's bytes, ORed.
#[inline]
fn classes(line: &[u8]) -> u8 {
    let (mut a, mut b, mut c, mut d) = (0u8, 0u8, 0u8, 0u8); // four, so bytes don't wait on each other
    let mut quads = line.chunks_exact(4);
    for q in &mut quads {
        a |= CLASS[q[0] as usize];
        b |= CLASS[q[1] as usize];
        c |= CLASS[q[2] as usize];
        d |= CLASS[q[3] as usize];
    }
    for &x in quads.remainder() {
        a |= CLASS[x as usize];
    }
    a | b | c | d
}

/// What an illegal residue's error reports.
struct Illegal {
    pos: u64,
    seen: [bool; 256],
    preview: Vec<u8>,
    /// residues the preview still lacks after the illegal one
    short: usize,
}

/// Reads a genome's FASTA from the blocks it is fed and reports each scaffold
/// as events: `("H", header)`, `("C", chunk)` for every full chunk of its
/// residues, then `("E", length, saw_n, runs, last_chunk)`, or
/// `("X", length, offset, characters, preview)` when it holds a residue that is
/// no nucleotide code. `runs` are its N-free stretches of `min_len` or more.
#[pyclass]
pub struct GenomeScanner {
    chunk_size: usize,
    min_len: u64,
    /// a header line the block edge cut, without its `>`
    carry: Vec<u8>,
    in_header: bool,
    line_start: bool,
    in_record: bool,
    /// nothing has followed the last header's line
    bare: bool,
    chunk: Vec<u8>,
    /// the end of the chunk before this one
    before: Vec<u8>,
    length: u64,
    saw_n: bool,
    run_start: Option<u64>,
    runs: Vec<(u64, u64)>,
    illegal: Option<Illegal>,
}

impl GenomeScanner {
    fn close_run(&mut self, end: u64) {
        if let Some(start) = self.run_start.take() {
            if end - start >= self.min_len {
                self.runs.push((start, end));
            }
        }
    }

    /// Follow the N-free runs through a line that holds an `N`.
    fn track_runs(&mut self, line: &[u8]) {
        let (start, n) = (self.length, line.len());
        self.saw_n = true;
        let (mut pos, mut last) = (0, 0);
        while pos < n {
            let Some(skip) = line[pos..].iter().position(|&b| b != b'N') else { break };
            let first = pos + skip;
            let end = memchr(b'N', &line[first..]).map_or(n, |k| first + k);
            if first > 0 {
                self.close_run(start + last as u64);
            }
            if self.run_start.is_none() {
                self.run_start = Some(start + first as u64);
            }
            last = end;
            pos = end;
        }
        if last < n {
            self.close_run(start + last as u64);
        }
    }

    fn store(&mut self, py: Python<'_>, mut line: &[u8], events: &Bound<'_, PyList>) -> PyResult<()> {
        while self.chunk.len() + line.len() >= self.chunk_size {
            let take = self.chunk_size - self.chunk.len();
            self.chunk.extend_from_slice(&line[..take]);
            events.append(("C", PyBytes::new(py, &self.chunk)))?;
            self.before.clear();
            self.before.extend_from_slice(&self.chunk[self.chunk.len().saturating_sub(PREVIEW)..]);
            self.chunk.clear();
            line = &line[take..];
        }
        self.chunk.extend_from_slice(line);
        Ok(())
    }

    /// The residues of one line, or of a piece of one.
    fn residues(&mut self, py: Python<'_>, line: &[u8], events: &Bound<'_, PyList>) -> PyResult<()> {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty() {
            return Ok(());
        }
        let class = classes(line);
        if class & BAD == 0 && self.illegal.is_none() {
            if class & HAS_N != 0 {
                self.track_runs(line);
            } else if self.run_start.is_none() {
                self.run_start = Some(self.length);
            }
            self.store(py, line, events)?;
            self.length += line.len() as u64;
            return Ok(());
        }
        // a CR inside the line is dropped; anything else that is no residue is illegal
        let cleaned: Vec<u8> = line.iter().copied().filter(|&b| b != b'\r').collect();
        if self.illegal.is_none() && classes(&cleaned) & BAD == 0 {
            return self.residues(py, &cleaned, events);
        }
        match &mut self.illegal {
            None => {
                let at = cleaned.iter().position(|&b| CLASS[b as usize] == BAD).unwrap_or(0);
                let mut preview: Vec<u8> = self.before.iter().chain(&self.chunk).chain(&cleaned[..at]).copied().collect();
                preview.drain(..preview.len().saturating_sub(PREVIEW));
                let shown = &cleaned[at..cleaned.len().min(at + PREVIEW + 1)];
                preview.extend_from_slice(shown);
                let mut seen = [false; 256];
                for &b in cleaned.iter().filter(|&&b| CLASS[b as usize] == BAD) {
                    seen[b as usize] = true;
                }
                self.illegal = Some(Illegal { pos: self.length + at as u64, seen, preview, short: PREVIEW + 1 - shown.len() });
            }
            Some(bad) => {
                let more = &cleaned[..cleaned.len().min(bad.short)];
                bad.preview.extend_from_slice(more);
                bad.short -= more.len();
                for &b in cleaned.iter().filter(|&&b| CLASS[b as usize] == BAD) {
                    bad.seen[b as usize] = true;
                }
            }
        }
        self.length += cleaned.len() as u64;
        Ok(())
    }

    /// Report the open scaffold, if any, and clear it.
    fn end_scaffold(&mut self, py: Python<'_>, events: &Bound<'_, PyList>) -> PyResult<()> {
        if !self.in_record {
            return Ok(());
        }
        if let Some(bad) = self.illegal.take() {
            let latin = |bytes: &[u8]| bytes.iter().map(|&b| b as char).collect::<String>();
            let characters: Vec<String> = (0..256usize).filter(|&b| bad.seen[b]).map(|b| (b as u8 as char).to_string()).collect();
            events.append(("X", self.length, bad.pos, characters, latin(&bad.preview)))?;
        } else {
            self.close_run(self.length);
            events.append(("E", self.length, self.saw_n, std::mem::take(&mut self.runs), PyBytes::new(py, &self.chunk)))?;
        }
        self.chunk.clear();
        self.before.clear();
        self.runs.clear();
        self.length = 0;
        self.saw_n = false;
        self.run_start = None;
        self.in_record = false;
        Ok(())
    }

    fn header(&mut self, line: &[u8], events: &Bound<'_, PyList>) -> PyResult<()> {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        let text = std::str::from_utf8(line).map_err(|e| PyValueError::new_err(format!("a header is not UTF-8: {e}")))?;
        events.append(("H", text))?;
        self.in_record = true;
        self.bare = true;
        Ok(())
    }
}

#[pymethods]
impl GenomeScanner {
    #[new]
    fn new(chunk_size: usize, min_len: u64) -> PyResult<Self> {
        if chunk_size == 0 {
            return Err(PyValueError::new_err("chunk_size must be positive"));
        }
        Ok(GenomeScanner {
            chunk_size,
            min_len,
            carry: Vec::new(),
            in_header: false,
            line_start: true,
            in_record: false,
            bare: false,
            chunk: Vec::with_capacity(chunk_size),
            before: Vec::with_capacity(PREVIEW),
            length: 0,
            saw_n: false,
            run_start: None,
            runs: Vec::new(),
            illegal: None,
        })
    }

    /// The events of the next block of the file.
    fn feed<'py>(&mut self, py: Python<'py>, block: &[u8]) -> PyResult<Bound<'py, PyList>> {
        let events = PyList::empty(py);
        let (mut pos, n) = (0, block.len());
        if self.in_header {
            let Some(k) = memchr(b'\n', block) else {
                self.carry.extend_from_slice(block);
                return Ok(events);
            };
            self.carry.extend_from_slice(&block[..k]);
            let line = std::mem::take(&mut self.carry);
            self.header(&line, &events)?;
            self.in_header = false;
            self.line_start = true;
            pos = k + 1;
        }
        while pos < n {
            self.bare = false;
            let end = memchr(b'\n', &block[pos..]).map(|k| pos + k);
            if self.line_start && block[pos] == b'>' {
                // the scaffold before is whole once the next one's header begins
                self.end_scaffold(py, &events)?;
                let Some(end) = end else {
                    self.carry.extend_from_slice(&block[pos + 1..]);
                    self.in_header = true;
                    break;
                };
                self.header(&block[pos + 1..end], &events)?;
                pos = end + 1;
                continue;
            }
            let line = &block[pos..end.unwrap_or(n)];
            if self.in_record {
                self.residues(py, line, &events)?;
            } else if line.iter().any(|b| !b.is_ascii_whitespace()) {
                return Err(PyValueError::new_err("sequence data before the first header"));
            }
            self.line_start = end.is_some();
            pos = end.map_or(n, |e| e + 1);
        }
        Ok(events)
    }

    /// The events left once the file is read through.
    fn finish<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        if self.in_header || self.bare {
            return Err(PyValueError::new_err("the file ends on a header"));
        }
        let events = PyList::empty(py);
        self.end_scaffold(py, &events)?;
        Ok(events)
    }
}

/// Windows of `window` nt, `overlap` shared between neighbours, over a
/// scaffold's `runs`; one shorter than `min_len` is left out. Returns their
/// count and, as native arrays, the columns prepare keeps: node, scaffold
/// index, start, end, scaffold length, window length less one, run start, run end.
#[pyfunction]
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn genome_windows<'py>(
    py: Python<'py>,
    runs: Vec<(i64, i64)>,
    parent_len: i64,
    window: i64,
    overlap: i64,
    min_len: i64,
    first_node: i64,
    scaffold_idx: i32,
) -> PyResult<(usize, Vec<Bound<'py, PyBytes>>)> {
    let step = window - overlap;
    if step <= 0 {
        return Err(PyValueError::new_err("overlap must be shorter than the window"));
    }
    let mut cols: [Vec<i64>; 7] = Default::default();
    for &(start, end) in &runs {
        let mut i = 0;
        while i < end - start {
            let child = window.min(end - start - i);
            if child >= min_len {
                let node = first_node + cols[0].len() as i64;
                for (col, value) in cols.iter_mut().zip([node, start + i, start + i + child - 1, parent_len, child - 1, start, end]) {
                    col.push(value);
                }
            }
            i += step;
        }
    }
    let count = cols[0].len();
    let bytes = |values: &[i64]| PyBytes::new(py, &values.iter().flat_map(|v| v.to_ne_bytes()).collect::<Vec<u8>>());
    let ids: Vec<u8> = std::iter::repeat_n(scaffold_idx.to_ne_bytes(), count).flatten().collect();
    let mut out = vec![bytes(&cols[0]), PyBytes::new(py, &ids)];
    out.extend(cols[1..].iter().map(|c| bytes(c)));
    Ok((count, out))
}
