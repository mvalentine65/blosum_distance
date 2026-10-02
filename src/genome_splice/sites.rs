//! Base codes, codons and splice-site scores.
//!
//! The donor PSSM is exonerate 2.4.0's (splice.c; Senapathy et al. 1990 primate
//! frequencies); the acceptor's is pooled from 12 arthropod genomes. Both are
//! log-odds scaled as exonerate does, rounded to integers.

pub const HMM_AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
const CODE: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"; // TCAG order

pub const RES_X: usize = 20;
pub const RES_STOP: usize = 21;

// A C G T
const SS5: [[i32; 4]; 9] = [
    [28, 40, 17, 14], [59, 14, 13, 14], [8, 5, 81, 6],
    [0, 0, 100, 0], [0, 0, 0, 100],
    [54, 2, 42, 2], [74, 8, 11, 8], [5, 6, 85, 4], [16, 18, 21, 45],
];
const SS5_AFTER: i64 = 3;
const SS3: [[i32; 4]; 15] = [
    [24, 16, 11, 48], [24, 16, 11, 50], [22, 16, 11, 52], [20, 16, 10, 54],
    [20, 18, 9, 53], [19, 19, 10, 51], [21, 19, 11, 49], [21, 17, 12, 50],
    [10, 14, 6, 70], [9, 14, 4, 73], [27, 13, 21, 39], [4, 58, 0, 38],
    [100, 0, 0, 0], [0, 0, 100, 0],
    [31, 14, 41, 14],
];
const SS3_LAST: i64 = 13;

/// A0 C1 G2 T3, anything else 4.
#[inline]
pub fn base(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

#[inline]
fn tcag(c: u8) -> i32 {
    match c {
        b'T' | b't' => 0,
        b'C' | b'c' => 1,
        b'A' | b'a' => 2,
        b'G' | b'g' => 3,
        _ => -1,
    }
}

/// Standard code; 'X' for a codon with an N.
#[inline]
pub fn translate(c: &[u8]) -> u8 {
    let (a, b, d) = (tcag(c[0]), tcag(c[1]), tcag(c[2]));
    if a < 0 || b < 0 || d < 0 {
        return b'X';
    }
    CODE[(a * 16 + b * 4 + d) as usize]
}

/// 0..19 HMMER order, X 20, stop 21.
pub fn residue(aa: u8) -> usize {
    if aa == b'*' {
        return RES_STOP;
    }
    HMM_AA.iter().position(|&x| x == aa).unwrap_or(RES_X)
}

/// ACGT codon index -> residue.
pub fn cod64() -> [usize; 64] {
    const ACGT: &[u8; 4] = b"ACGT";
    let mut out = [0usize; 64];
    for i in 0..64 {
        let c = [ACGT[i >> 4], ACGT[(i >> 2) & 3], ACGT[i & 3]];
        out[i] = residue(translate(&c));
    }
    out
}

fn lod(f: i32) -> f64 {
    ((1.0 + f as f64) / 26.0).ln() * 1.5
}

fn rnd(x: f64) -> f64 {
    if x < 0.0 { (x - 0.5).ceil() } else { (x + 0.5).floor() }
}

/// Donor score at p (the GT starts at p) and acceptor score at p (the AG ends
/// at p). Columns off the sequence add nothing; N adds 0.
pub fn splice_scores(b: &[u8], ss5: &mut [f64], ss3: &mut [f64]) {
    let d = b.len() as i64;
    for p in 0..d {
        let (mut s5, mut s3) = (0.0, 0.0);
        for (i, row) in SS5.iter().enumerate() {
            let j = p - SS5_AFTER + i as i64;
            if j >= 0 && j < d && b[j as usize] < 4 {
                s5 += lod(row[b[j as usize] as usize]);
            }
        }
        for (i, row) in SS3.iter().enumerate() {
            let j = p - SS3_LAST + i as i64;
            if j >= 0 && j < d && b[j as usize] < 4 {
                s3 += lod(row[b[j as usize] as usize]);
            }
        }
        ss5[p as usize] = rnd(s5);
        ss3[p as usize] = rnd(s3);
    }
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C',
            b'a' => b't', b't' => b'a', b'c' => b'g', b'g' => b'c',
            _ => b'N',
        })
        .collect()
}
