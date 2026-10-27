//! Genome chunks at 2 bits per base ("chunk-v2" pseq store).
//!
//! A chunk is stored as a 12-byte header (u32 length, u32 non-ACGT runs, u32
//! lowercase runs), its bases packed 4 to a byte (A 0, C 1, G 2, T 3, first
//! base in the low bits; 0 where the base is not ACGT), then the non-ACGT runs
//! (u16 start, u16 length, u8 uppercase byte) and the lowercase runs (u16
//! start, u16 length). All little-endian; a run length of 0 means 65536.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

/// Longest chunk a value can hold: run starts and lengths are u16.
pub const MAX_CHUNK: usize = 1 << 16;
const HEADER: usize = 12;
const EX_RUN: usize = 5;
const CASE_RUN: usize = 4;

/// Per byte: its 2-bit code (0 where not ACGT), OTHER for any other letter,
/// LOWER for lowercase, BAD for anything that is not a letter.
const OTHER: u8 = 4;
const LOWER: u8 = 8;
const BAD: u8 = 0x80;
static CLASS: [u8; 256] = {
    let mut t = [BAD; 256];
    let mut c = b'A';
    while c <= b'Z' {
        t[c as usize] = OTHER;
        t[(c | 0x20) as usize] = OTHER | LOWER;
        c += 1;
    }
    let bases = *b"ACGT";
    let mut i = 0;
    while i < 4 {
        t[bases[i] as usize] = i as u8;
        t[(bases[i] | 0x20) as usize] = i as u8 | LOWER;
        i += 1;
    }
    t
};

/// The four bases a packed byte holds, first base in the low bits.
static BASES: [[u8; 4]; 256] = {
    let mut t = [[0u8; 4]; 256];
    let mut byte = 0;
    while byte < 256 {
        let mut i = 0;
        while i < 4 {
            t[byte][i] = b"ACGT"[(byte >> (2 * i)) & 3];
            i += 1;
        }
        byte += 1;
    }
    t
};

fn put_run(out: &mut Vec<u8>, start: usize, len: usize) {
    out.extend_from_slice(&(start as u16).to_le_bytes());
    out.extend_from_slice(&((len % MAX_CHUNK) as u16).to_le_bytes());
}

/// Runs of non-ACGT letters and of lowercase, as `pack` finds them.
struct Runs {
    ex: Vec<u8>,
    case: Vec<u8>,
    n_ex: u32,
    n_case: u32,
    // open runs: start, or usize::MAX when none is open
    ex_start: usize,
    ex_byte: u8,
    case_start: usize,
}

impl Runs {
    fn close_ex(&mut self, end: usize) {
        put_run(&mut self.ex, self.ex_start, end - self.ex_start);
        self.ex.push(self.ex_byte);
        self.n_ex += 1;
        self.ex_start = usize::MAX;
    }

    fn close_case(&mut self, end: usize) {
        put_run(&mut self.case, self.case_start, end - self.case_start);
        self.n_case += 1;
        self.case_start = usize::MAX;
    }

    /// Record base `i` (class `k`, byte `b`).
    fn step(&mut self, i: usize, b: u8, k: u8) {
        if k & OTHER != 0 {
            let up = b & 0xDF;
            if self.ex_start != usize::MAX && up != self.ex_byte {
                self.close_ex(i);
            }
            if self.ex_start == usize::MAX {
                (self.ex_start, self.ex_byte) = (i, up);
            }
        } else if self.ex_start != usize::MAX {
            self.close_ex(i);
        }
        if k & LOWER != 0 {
            if self.case_start == usize::MAX {
                self.case_start = i;
            }
        } else if self.case_start != usize::MAX {
            self.close_case(i);
        }
    }
}

const LANES: u64 = 0x0101_0101_0101_0101;

/// High bit set in each byte of `y` that is zero, and only there.
#[inline]
fn zero_bytes(y: u64) -> u64 {
    !(((y & (0x7F * LANES)) + 0x7F * LANES) | y) & (0x80 * LANES)
}

/// Eight bases packed to 16 bits, when all are ACGT in either case; else None.
#[inline]
fn pack8(x: u64) -> Option<u16> {
    let u = x & (0xDF * LANES);
    let hits = zero_bytes(u ^ (b'A' as u64 * LANES))
        | zero_bytes(u ^ (b'C' as u64 * LANES))
        | zero_bytes(u ^ (b'G' as u64 * LANES))
        | zero_bytes(u ^ (b'T' as u64 * LANES));
    if hits != 0x80 * LANES {
        return None;
    }
    // bits 1-2 of a base: A 0, C 1, T 2, G 3; swapping T and G gives the stored A 0, C 1, G 2, T 3
    let t = (u >> 1) & (3 * LANES);
    let c = t ^ ((t >> 1) & LANES);
    let v = (c | (c >> 6)) & 0x000F_000F_000F_000F;
    let v = (v | (v >> 12)) & 0x0000_00FF_0000_00FF;
    Some((v | (v >> 24)) as u16)
}

/// `chunk` (ASCII letters only, 1 to MAX_CHUNK of them) in its stored form, in one pass.
pub fn pack(chunk: &[u8]) -> Result<Vec<u8>, String> {
    let n = chunk.len();
    if n == 0 || n > MAX_CHUNK {
        return Err(format!("pseq chunk of {n} bases; expected 1 to {MAX_CHUNK}"));
    }
    let mut out = vec![0u8; HEADER + n.div_ceil(4)];
    let mut runs = Runs {
        ex: Vec::new(),
        case: Vec::new(),
        n_ex: 0,
        n_case: 0,
        ex_start: usize::MAX,
        ex_byte: 0,
        case_start: usize::MAX,
    };
    let mut i = 0;
    while i < n {
        // Eight ACGT with no non-ACGT run open: pack them at once; only case edges to record.
        if i % 4 == 0 && i + 8 <= n && runs.ex_start == usize::MAX {
            let x = u64::from_le_bytes(chunk[i..i + 8].try_into().unwrap());
            if let Some(two) = pack8(x) {
                out[HEADER + i / 4..HEADER + i / 4 + 2].copy_from_slice(&two.to_le_bytes());
                let lower = x & (0x20 * LANES);
                if lower == 0 {
                    if runs.case_start != usize::MAX {
                        runs.close_case(i);
                    }
                } else if lower == 0x20 * LANES {
                    if runs.case_start == usize::MAX {
                        runs.case_start = i;
                    }
                } else {
                    for (j, &b) in chunk[i..i + 8].iter().enumerate() {
                        runs.step(i + j, b, CLASS[b as usize]);
                    }
                }
                i += 8;
                continue;
            }
        }
        let b = chunk[i];
        let k = CLASS[b as usize];
        if k & BAD != 0 {
            return Err(format!("pseq chunk holds byte {b:#04x} at {i}; only letters can be packed"));
        }
        out[HEADER + i / 4] |= (k & 3) << (2 * (i % 4));
        runs.step(i, b, k);
        i += 1;
    }
    if runs.ex_start != usize::MAX {
        runs.close_ex(n);
    }
    if runs.case_start != usize::MAX {
        runs.close_case(n);
    }
    out[0..4].copy_from_slice(&(n as u32).to_le_bytes());
    out[4..8].copy_from_slice(&runs.n_ex.to_le_bytes());
    out[8..12].copy_from_slice(&runs.n_case.to_le_bytes());
    out.extend_from_slice(&runs.ex);
    out.extend_from_slice(&runs.case);
    Ok(out)
}

fn u32_at(v: &[u8], p: usize) -> usize {
    u32::from_le_bytes(v[p..p + 4].try_into().unwrap()) as usize
}

fn run_at(v: &[u8], p: usize) -> (usize, usize) {
    let s = u16::from_le_bytes([v[p], v[p + 1]]) as usize;
    let n = u16::from_le_bytes([v[p + 2], v[p + 3]]) as usize;
    (s, if n == 0 { MAX_CHUNK } else { n })
}

/// Bases a stored value holds, after checking its layout.
pub fn stored_len(value: &[u8]) -> Result<usize, String> {
    if value.len() < HEADER {
        return Err(format!("pseq value of {} bytes is shorter than its header", value.len()));
    }
    let (n, n_ex, n_case) = (u32_at(value, 0), u32_at(value, 4), u32_at(value, 8));
    if n == 0 || n > MAX_CHUNK {
        return Err(format!("pseq value claims {n} bases; expected 1 to {MAX_CHUNK}"));
    }
    let want = HEADER + n.div_ceil(4) + EX_RUN * n_ex + CASE_RUN * n_case;
    if value.len() != want {
        return Err(format!("pseq value is {} bytes; its header implies {want}", value.len()));
    }
    Ok(n)
}

/// Append the bases `value` stores to `out`; `upper` leaves out the lowercase runs.
pub fn unpack_into(value: &[u8], upper: bool, out: &mut Vec<u8>) -> Result<(), String> {
    let n = stored_len(value)?;
    let (n_ex, n_case) = (u32_at(value, 4), u32_at(value, 8));
    let base = out.len();
    let mut p = HEADER + n.div_ceil(4);
    out.reserve(n.div_ceil(4) * 4);
    for &b in &value[HEADER..p] {
        out.extend_from_slice(&BASES[b as usize]);
    }
    out.truncate(base + n);
    let seq = &mut out[base..];
    for _ in 0..n_ex {
        let (s, len) = run_at(value, p);
        let b = value[p + 4];
        if s + len > n || !b.is_ascii_uppercase() {
            return Err(format!("pseq value has a bad non-ACGT run at {s}+{len} ({b:#04x})"));
        }
        seq[s..s + len].fill(b);
        p += EX_RUN;
    }
    for _ in 0..n_case {
        let (s, len) = run_at(value, p);
        if s + len > n {
            return Err(format!("pseq value has a lowercase run past its end at {s}+{len}"));
        }
        if !upper {
            for c in &mut seq[s..s + len] {
                *c |= 0x20;
            }
        }
        p += CASE_RUN;
    }
    Ok(())
}

/// A genome chunk in its 2-bit stored form.
#[pyfunction]
pub fn pack_pseq<'py>(py: Python<'py>, chunk: &[u8]) -> PyResult<Bound<'py, PyBytes>> {
    let out = pack(chunk).map_err(PyValueError::new_err)?;
    Ok(PyBytes::new(py, &out))
}

/// The bases a 2-bit stored chunk holds; `upper` returns them uppercased.
#[pyfunction]
#[pyo3(signature = (value, upper = false))]
pub fn unpack_pseq<'py>(py: Python<'py>, value: &[u8], upper: bool) -> PyResult<Bound<'py, PyBytes>> {
    let mut out = Vec::new();
    unpack_into(value, upper, &mut out).map_err(PyValueError::new_err)?;
    Ok(PyBytes::new(py, &out))
}

/// Several 2-bit stored chunks, unpacked and joined in order.
#[pyfunction]
#[pyo3(signature = (values, upper = false))]
pub fn unpack_pseq_many<'py>(
    py: Python<'py>,
    values: Vec<Bound<'py, PyBytes>>,
    upper: bool,
) -> PyResult<Bound<'py, PyBytes>> {
    let mut total = 0;
    for v in &values {
        total += stored_len(v.as_bytes()).map_err(PyValueError::new_err)?;
    }
    let mut out = Vec::with_capacity(total + 3);
    for v in &values {
        unpack_into(v.as_bytes(), upper, &mut out).map_err(PyValueError::new_err)?;
    }
    Ok(PyBytes::new(py, &out))
}

#[cfg(test)]
#[path = "tests/pseq.rs"]
mod tests;
