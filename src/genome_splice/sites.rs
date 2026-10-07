//! Base codes, codons and splice-site scores.
//!
//! The donor PSSM is the primate frequencies of Senapathy et al. 1990 (Methods
//! in Enzymology 183:252-278), as tabulated in exonerate 2.4.0; the acceptor's
//! is pooled from 12 arthropod genomes, or learned from the run's own junctions
//! (learn_acceptor). Both are log-odds scaled as exonerate does, rounded to
//! integers.

pub const HMM_AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
/// A genetic code: the residue of each codon in TCAG order, as NCBI writes its tables.
pub type Code = &'static [u8; 64];
pub const STANDARD: Code = crate::translate::TABLE_1;

/// The NCBI table with this id.
pub fn code(table: u8) -> Option<Code> {
    crate::translate::genetic_code_table(table)
}

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
pub const SS3_LAST: i64 = 13;

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

/// 'X' for a codon with an N.
#[inline]
pub fn translate(code: Code, c: &[u8]) -> u8 {
    let (a, b, d) = (tcag(c[0]), tcag(c[1]), tcag(c[2]));
    if a < 0 || b < 0 || d < 0 {
        return b'X';
    }
    code[(a * 16 + b * 4 + d) as usize]
}

/// 0..19 HMMER order, X 20, stop 21.
pub fn residue(aa: u8) -> usize {
    if aa == b'*' {
        return RES_STOP;
    }
    HMM_AA.iter().position(|&x| x == aa).unwrap_or(RES_X)
}

/// ACGT codon index -> residue.
pub fn cod64(code: Code) -> [usize; 64] {
    const ACGT: &[u8; 4] = b"ACGT";
    let mut out = [0usize; 64];
    for i in 0..64 {
        let c = [ACGT[i >> 4], ACGT[(i >> 2) & 3], ACGT[i & 3]];
        out[i] = residue(translate(code, &c));
    }
    out
}

fn lod(f: i32) -> f64 {
    ((1.0 + f as f64) / 26.0).ln() * 1.5
}

/// Acceptor log-odds, one row per SS3 position.
pub type AccTable = [[f64; 4]; 15];

/// SS3 against equal base frequencies.
pub fn acceptor_default() -> AccTable {
    let mut t = [[0.0; 4]; 15];
    for (i, row) in SS3.iter().enumerate() {
        for b in 0..4 { t[i][b] = lod(row[b]); }
    }
    t
}

/// SS3 positions that share one base distribution when learned; the AG itself is not learned
const TIED: [&[usize]; 5] = [&[0, 1, 2, 3, 4, 5, 6, 7], &[8, 9], &[10], &[11], &[14]];
/// observations the prior counts for in a learned distribution
const LEARN_PRIOR: f64 = 50.0;

/// Acceptor log-odds of true sites against the false AGs around them, each window holding its
/// true G at flank + SS3_LAST (base codes). True sites are blended with SS3, false ones with
/// equal frequencies; true sites keep their mean score under SS3. None without a window.
pub fn learn_acceptor(wins: &[Vec<u8>], flank: usize) -> Option<AccTable> {
    if wins.is_empty() { return None; }
    let last = SS3_LAST as usize;
    let (mut tr, mut de) = ([[0f64; 4]; 15], [[0f64; 4]; 15]);
    for w in wins {
        for j in last..w.len().saturating_sub(1) {
            if w[j - 1] != 0 || w[j] != 2 { continue; }
            let c = if j == flank + last { &mut tr } else { &mut de };
            for (i, row) in c.iter_mut().enumerate() {
                let b = w[j - last + i];
                if b < 4 { row[b as usize] += 1.0; }
            }
        }
    }
    let base = acceptor_default();
    let mut t = base;
    for g in TIED {
        let pct = |c: &[[f64; 4]; 15]| -> ([f64; 4], f64) {
            let mut s = [0f64; 4];
            for &i in g.iter() { for b in 0..4 { s[b] += c[i][b]; } }
            let n: f64 = s.iter().sum();
            (s.map(|x| if n > 0.0 { 100.0 * x / n } else { 0.0 }), n)
        };
        let ((ft, nt), (fd, nd)) = (pct(&tr), pct(&de));
        for &i in g.iter() {
            for b in 0..4 {
                let yes = (ft[b] * nt + LEARN_PRIOR * SS3[i][b] as f64) / (nt + LEARN_PRIOR);
                let no = (fd[b] * nd + LEARN_PRIOR * 25.0) / (nd + LEARN_PRIOR);
                t[i][b] = ((1.0 + yes) / (1.0 + no)).ln() * 1.5;
            }
        }
    }
    let mean = |m: &AccTable| -> f64 {
        wins.iter().map(|w| (0..15).map(|i| { let b = w[flank + i]; if b < 4 { m[i][b as usize] } else { 0.0 } }).sum::<f64>()).sum::<f64>() / wins.len() as f64
    };
    t[last][2] += mean(&base) - mean(&t);
    Some(t)
}

fn rnd(x: f64) -> f64 {
    if x < 0.0 { (x - 0.5).ceil() } else { (x + 0.5).floor() }
}

/// Donor score at p (the GT starts at p) and acceptor score at p (the AG ends
/// at p) under acc. Columns off the sequence add nothing; N adds 0.
pub fn splice_scores(b: &[u8], acc: &AccTable, ss5: &mut [f64], ss3: &mut [f64]) {
    let d = b.len() as i64;
    for p in 0..d {
        let (mut s5, mut s3) = (0.0, 0.0);
        for (i, row) in SS5.iter().enumerate() {
            let j = p - SS5_AFTER + i as i64;
            if j >= 0 && j < d && b[j as usize] < 4 {
                s5 += lod(row[b[j as usize] as usize]);
            }
        }
        for (i, row) in acc.iter().enumerate() {
            let j = p - SS3_LAST + i as i64;
            if j >= 0 && j < d && b[j as usize] < 4 {
                s3 += row[b[j as usize] as usize];
            }
        }
        ss5[p as usize] = rnd(s5);
        ss3[p as usize] = rnd(s3);
    }
}

/// Donor score of the intron after scaffold base don_g and acceptor score of the
/// one before acc_g (1-based, the gene read on its strand), under acc.
pub fn site_scores(sq: &[u8], plus: bool, don_g: i64, acc_g: i64, acc: &AccTable) -> (f64, f64) {
    // n bases as the gene reads them from scaffold base g; N off the scaffold
    let read = |g: i64, n: i64| -> Vec<u8> {
        (0..n).map(|i| {
            let p = if plus { g + i } else { g - i };
            if p < 1 || p > sq.len() as i64 { return 4; }
            let b = base(sq[(p - 1) as usize]);
            if plus || b > 3 { b } else { 3 - b }
        }).collect()
    };
    let dir = if plus { 1 } else { -1 };
    let score = |t: &[u8]| {
        let (mut ss5, mut ss3) = (vec![0f64; t.len()], vec![0f64; t.len()]);
        splice_scores(t, acc, &mut ss5, &mut ss3);
        (ss5, ss3)
    };
    // the donor's G follows SS5_AFTER exon bases; the acceptor's G is one before the exon
    let don = score(&read(don_g - dir * (SS5_AFTER - 1), SS5.len() as i64)).0[SS5_AFTER as usize];
    let ac = score(&read(acc_g - dir * (SS3_LAST + 1), SS3.len() as i64)).1[SS3_LAST as usize];
    (don, ac)
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
