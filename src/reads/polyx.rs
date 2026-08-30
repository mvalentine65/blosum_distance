//! polyG and polyX tail trimming.
//!
//! Ported from fastp 1.3.6 `src/polyx.cpp` (MIT, (c) 2016 OpenGene).
//!
//! On two-colour chemistry (NovaSeq, NextSeq) a dark cycle reads as G, so a
//! read that runs off the end of its cluster finishes in a run of Gs. polyX is
//! the same idea generalised to whichever base dominates the tail.

use crate::reads::read::ReadRec;

const ALLOW_ONE_MISMATCH_FOR_EACH: usize = 8;
const MAX_MISMATCH: usize = 5;

/// fastp's `ATCG_BASES`, and the order its `atcgNumbers` tallies use.
const ATCG_BASES: [u8; 4] = [b'A', b'T', b'C', b'G'];

/// fastp's `POLYX_BASE_IDX`: A=0, T=1, C=2, G=3, N=4, anything else 5.
/// Uppercase only, exactly as fastp has it — reads reaching this point come
/// straight from an Illumina FASTQ.
#[inline]
fn base_idx(b: u8) -> u8 {
    match b {
        b'A' => 0,
        b'T' => 1,
        b'C' => 2,
        b'G' => 3,
        b'N' => 4,
        _ => 5,
    }
}

/// Trim a polyG tail. `compare_req` is the shortest run worth cutting.
///
/// Returns the number of bases removed.
pub fn trim_poly_g(r: &mut ReadRec<'_>, compare_req: usize) -> usize {
    let data = r.bases();
    let rlen = data.len();
    if rlen == 0 {
        return 0;
    }

    let mut mismatch = 0usize;
    let mut first_g_pos = rlen - 1;
    let mut i = 0usize;
    while i < rlen {
        if data[rlen - i - 1] != b'G' {
            mismatch += 1;
        } else {
            first_g_pos = rlen - i - 1;
        }
        let allowed_mismatch = (i + 1) / ALLOW_ONE_MISMATCH_FOR_EACH;
        if mismatch > MAX_MISMATCH || (mismatch > allowed_mismatch && i + 1 >= compare_req) {
            break;
        }
        i += 1;
    }

    if i >= compare_req {
        let removed = rlen - first_g_pos;
        r.resize(first_g_pos);
        removed
    } else {
        0
    }
}

/// Trim a tail dominated by any single base. Returns `(base, bases removed)`.
pub fn trim_poly_x(r: &mut ReadRec<'_>, compare_req: usize) -> Option<(u8, usize)> {
    let data = r.bases();
    let rlen = data.len();
    if rlen == 0 {
        return None;
    }

    let mut atcg = [0usize; 4];
    let mut pos = 0usize;
    while pos < rlen {
        let idx = base_idx(data[rlen - pos - 1]);
        if idx < 4 {
            atcg[idx as usize] += 1;
        } else if idx == 4 {
            // fastp lets an N count toward every base
            for slot in atcg.iter_mut() {
                *slot += 1;
            }
        }

        let cmp = pos + 1;
        let allowed_mismatch = MAX_MISMATCH.min(cmp / ALLOW_ONE_MISMATCH_FOR_EACH);

        let mut need_to_break = true;
        for b in 0..4 {
            if cmp - atcg[b] <= allowed_mismatch {
                need_to_break = false;
            }
        }
        if need_to_break
            && (pos >= ALLOW_ONE_MISMATCH_FOR_EACH || pos + 1 >= compare_req.saturating_sub(1))
        {
            break;
        }
        pos += 1;
    }

    if pos + 1 < compare_req {
        return None;
    }

    let mut poly = 0usize;
    let mut max_count = 0usize;
    for b in 0..4 {
        if atcg[b] > max_count {
            max_count = atcg[b];
            poly = b;
        }
    }
    let poly_base = ATCG_BASES[poly];

    // Walk back to the last position actually holding the poly base. fastp
    // tests `data[rlen-pos-1]` before the bounds check and so reads the string
    // terminator at pos == -1; we check the bound first, which lands on the
    // same outcome (no trim) without the out-of-range read.
    // `pos` reaches rlen when the scan never broke, i.e. the whole read is the
    // run; the walk-back has to start at the last base that exists.
    let mut back = pos.min(rlen - 1) as isize;
    while back >= 0 && data[rlen - back as usize - 1] != poly_base {
        back -= 1;
    }
    if back < 0 {
        return None;
    }
    let removed = back as usize + 1;
    r.resize(rlen - removed);
    Some((poly_base, removed))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// fastp's own `PolyX::test()`, expectations included.
    #[test]
    fn fastp_polyx_test_case() {
        let seq = b"ATTTTAAAAAAAAAATAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAT";
        let mut r = ReadRec::from_seq(seq);
        let trimmed = trim_poly_x(&mut r, 10);
        assert_eq!(r.bases(), b"ATTTT");
        let (base, removed) = trimmed.expect("a polyA tail should be found");
        assert_eq!(base, b'A');
        assert_eq!(removed, 51);
    }

    /// The cut lands on the leftmost G that the mismatch budget can still
    /// reach, not on the first G of the visible run: walking back over the 12
    /// Gs costs nothing, the `T` at 15 is the one allowed mismatch, and the
    /// `G` at 14 then becomes the new `firstGPos`. So 14 bases go, not 12.
    /// That is fastp's rule, and it is why polyG trimming can take a base or
    /// two of real insert with it.
    #[test]
    fn trims_a_polyg_tail_back_to_the_reachable_g() {
        let seq = b"ACGTACGTACGTACGTGGGGGGGGGGGG";
        let mut r = ReadRec::from_seq(seq);
        let removed = trim_poly_g(&mut r, 10);
        assert_eq!(removed, 14);
        assert_eq!(r.bases(), b"ACGTACGTACGTAC");
    }

    #[test]
    fn trims_a_clean_polyg_tail_exactly() {
        // No G adjacent to the run, so the cut is exactly the run.
        let seq = b"ACGTACGTACGTACAAGGGGGGGGGGGG";
        let mut r = ReadRec::from_seq(seq);
        assert_eq!(trim_poly_g(&mut r, 10), 12);
        assert_eq!(r.bases(), b"ACGTACGTACGTACAA");
    }

    #[test]
    fn leaves_a_clean_read_alone() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGT";
        let mut r = ReadRec::from_seq(seq);
        assert_eq!(trim_poly_g(&mut r, 10), 0);
        assert_eq!(r.bases(), seq);
        let mut r2 = ReadRec::from_seq(seq);
        assert!(trim_poly_x(&mut r2, 10).is_none());
        assert_eq!(r2.bases(), seq);
    }

    #[test]
    fn empty_read_is_safe() {
        let mut r = ReadRec::from_seq(b"");
        assert_eq!(trim_poly_g(&mut r, 10), 0);
        assert!(trim_poly_x(&mut r, 10).is_none());
    }
}
