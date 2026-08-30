//! Tandem-repeat detection.
//!
//! Not from fastp: its low-complexity filter measures compositional diversity,
//! which conflates a tandem repeat with a sequence that is merely base-biased.
//! A repeat is self-similarity at a lag, so that is what this measures -- for
//! each lag `p`, the share of positions where `s[i] == s[i - p]`, less what
//! the sequence's own base composition would produce by chance:
//!
//! ```text
//! E        = sum of squared base frequencies   (chance match rate)
//! match(p) = |{ i : s[i] == s[i-p] }| / (len - p)
//! excess   = max over p of (match(p) - E) / (1 - E)
//! ```
//!
//! Normalising by `E` is the point: an AT-rich aperiodic read scores near zero
//! rather than near 0.34, which no entropy threshold achieves. 205 ns per
//! sequence, flagging ~5-9% depending on the library.

/// 2-bit base code, 0xFF for anything else. Local so the module stands alone.
const BASE_CODE: [u8; 256] = {
    let mut t = [0xFFu8; 256];
    t[b'A' as usize] = 0; t[b'a' as usize] = 0;
    t[b'C' as usize] = 1; t[b'c' as usize] = 1;
    t[b'G' as usize] = 2; t[b'g' as usize] = 2;
    t[b'T' as usize] = 3; t[b't' as usize] = 3;
    t
};

/// Longest lag considered. Being a multiple of 2, 3, 4 and 6 makes 12 the
/// catch-all for degraded short satellites that miss at their own period.
/// Going further costs ~33% more at 135 bp and found almost nothing new.
pub const MAX_PERIOD: usize = 12;

/// Excess self-similarity required to call a sequence a repeat. Aperiodic
/// sequence sits near 0, so there is a wide empty margin below this.
pub const DEFAULT_THRESHOLD: f64 = 0.8;

/// Positions where `s` agrees with itself `p` bases earlier. Chunked like the
/// adapter scan, which is what lets it vectorise.
#[inline]
fn matches_at_lag(s: &[u8], p: usize) -> usize {
    let (a, b) = (&s[p..], &s[..s.len() - p]);
    let mut same = 0usize;
    let mut ca = a.chunks_exact(16);
    let mut cb = b.chunks_exact(16);
    for (x, y) in ca.by_ref().zip(cb.by_ref()) {
        same += x.iter().zip(y).filter(|(l, r)| l == r).count();
    }
    same + ca
        .remainder()
        .iter()
        .zip(cb.remainder())
        .filter(|(l, r)| l == r)
        .count()
}

/// The repeat period, if this sequence is one. `None` means it is not.
///
/// The shortest qualifying lag: the unit length for a clean repeat, a multiple
/// of it for one degraded enough to miss at its own period. Stops at the first
/// hit, so a repeat costs less than a non-repeat.
pub fn repeat_period(s: &[u8], threshold: f64) -> Option<usize> {
    let l = s.len();
    if l < 4 {
        return None;
    }

    // Chance match rate from this sequence's own composition.
    let mut freq = [0u32; 4];
    let mut n = 0u32;
    for &b in s {
        let c = BASE_CODE[b as usize];
        if c != 0xFF {
            freq[c as usize] += 1;
            n += 1;
        }
    }
    if n == 0 {
        return None;
    }
    let e: f64 = freq
        .iter()
        .map(|&c| {
            let p = f64::from(c) / f64::from(n);
            p * p
        })
        .sum();
    if e >= 1.0 {
        // A single base throughout: a homopolymer, and the formula degenerates.
        return Some(1);
    }

    // excess >= threshold  <=>  matches >= (E + threshold * (1 - E)) * (len - p)
    // Comparing counts against this bound avoids a divide inside the loop.
    let need = e + threshold * (1.0 - e);
    for p in 1..=MAX_PERIOD.min(l - 1) {
        if matches_at_lag(s, p) as f64 >= need * (l - p) as f64 {
            return Some(p);
        }
    }
    None
}

/// The repeat unit of `s` at period `p`, as a class not a literal substring:
/// consensus over the copies, reduced to its primitive period, then folded
/// across rotation and strand. `(AG)n` and `(AGAG)n` are one satellite.
pub fn canonical_unit(s: &[u8], p: usize) -> Vec<u8> {
    let p = p.min(s.len()).min(MAX_PERIOD);
    if p == 0 {
        return Vec::new();
    }

    // Majority base at each offset within the unit, over all copies.
    let mut counts = [[0u32; 4]; MAX_PERIOD];
    for (i, &b) in s.iter().enumerate() {
        let c = BASE_CODE[b as usize];
        if c != 0xFF {
            counts[i % p][c as usize] += 1;
        }
    }
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let consensus: Vec<u8> = (0..p)
        .map(|j| {
            let col = &counts[j];
            if col.iter().all(|&n| n == 0) {
                return b'N';
            }
            // Ties resolve to the lowest base so the key is deterministic.
            let best = (0..4).max_by_key(|&k| (col[k], std::cmp::Reverse(k))).unwrap();
            BASES[best]
        })
        .collect();

    // Smallest d dividing p that the consensus is built from: `AGAG` -> `AG`.
    let prim = (1..=p)
        .find(|d| p % d == 0 && (0..p).all(|i| consensus[i] == consensus[i % d]))
        .unwrap_or(p);
    let unit = &consensus[..prim];

    let rc: Vec<u8> = unit
        .iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect();

    let mut best: Option<Vec<u8>> = None;
    for src in [unit.to_vec(), rc] {
        for r in 0..prim {
            let cand: Vec<u8> = (0..prim).map(|i| src[(r + i) % prim]).collect();
            if best.as_ref().is_none_or(|b| cand < *b) {
                best = Some(cand);
            }
        }
    }
    best.unwrap_or_default()
}

#[cfg(test)]
#[path = "../tests/periodicity.rs"]
mod tests;
