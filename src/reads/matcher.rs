//! One-gap matching.
//!
//! Ported from fastp 1.3.6 `src/matcher.cpp` (MIT, (c) 2016 OpenGene).
//!
//! Asks: if `ins_data` carries exactly one extra base somewhere, can it be
//! aligned to `normal_data` within `diff_limit` mismatches? The answer comes
//! from a prefix scan and a suffix scan, so a split point can be evaluated in
//! O(1) at every position.
//!
//! fastp's companion `matchWithOneInsertion` is not ported: its only callers
//! were the adapter gap passes, which we dropped. This routine survives as the
//! engine behind the overlap analyser's optional gap path
//! (`--allow_gap_overlap_trimming`), which is off by default — so if that
//! option is never wanted, this whole module can go too.
//!
//! Fidelity note: fastp declares both accumulators as uninitialised C99 VLAs
//! and `break`s out of the fill loops early, so the untouched tail is read as
//! whatever was on the stack. Reading uninitialised memory is not expressible
//! in safe Rust, and the breaks only fire once the running total has already
//! passed `diff_limit`, so we seed the arrays with `diff_limit + 1` — the
//! value fastp's own suffix loop writes when it takes the same exit. Every
//! comparison against those cells is a `> diff_limit` test, so the outcome is
//! identical while the behaviour becomes defined.

/// Reusable accumulators. The DP is re-run at every candidate offset, so
/// allocating these per call dominated the whole trimmer: a read with no
/// adapter tries ~146 offsets in each of the two gap passes, and each call
/// used to heap-allocate two vectors. Keeping them in a scratch and refilling
/// costs the same memset without the malloc.
#[derive(Default)]
pub struct MatchScratch {
    left: Vec<i32>,
    right: Vec<i32>,
}

/// Prefix and suffix mismatch accumulators for a one-insertion alignment.
///
/// `acc_left[i]`  — mismatches of `ins_data[0..=i]` against `normal_data[0..=i]`
/// `acc_right[i]` — mismatches of `ins_data[i+1..=cmplen]` against `normal_data[i..]`
///
/// A split at `i` therefore costs `acc_left[i-1] + acc_right[i]`.
fn accumulate(
    scratch: &mut MatchScratch,
    ins_data: &[u8],
    normal_data: &[u8],
    cmplen: usize,
    diff_limit: i32,
) -> bool {
    if cmplen == 0 || ins_data.len() < cmplen + 1 || normal_data.len() < cmplen {
        return false;
    }
    let over = diff_limit + 1;
    // clear + resize refills every slot with the sentinel and reuses the
    // capacity, so this allocates only on the first (and longest) call.
    scratch.left.clear();
    scratch.left.resize(cmplen, over);
    scratch.right.clear();
    scratch.right.resize(cmplen, over);
    let acc_left = &mut scratch.left;
    let acc_right = &mut scratch.right;

    acc_left[0] = i32::from(ins_data[0] != normal_data[0]);
    acc_right[cmplen - 1] = i32::from(ins_data[cmplen] != normal_data[cmplen - 1]);

    for i in 1..cmplen {
        acc_left[i] = acc_left[i - 1] + i32::from(ins_data[i] != normal_data[i]);
        if acc_left[i] + acc_right[cmplen - 1] > diff_limit {
            break;
        }
    }
    for i in (0..cmplen - 1).rev() {
        acc_right[i] = acc_right[i + 1] + i32::from(ins_data[i + 1] != normal_data[i]);
        if acc_right[i] + acc_left[0] > diff_limit {
            for slot in acc_right.iter_mut().take(i) {
                *slot = over;
            }
            break;
        }
    }
    true
}

/// Smallest mismatch count achievable with one insertion, or `-1` when every
/// split point exceeds `diff_limit`.
pub fn diff_with_one_insertion(
    scratch: &mut MatchScratch,
    ins_data: &[u8],
    normal_data: &[u8],
    cmplen: usize,
    diff_limit: i32,
) -> i32 {
    if !accumulate(scratch, ins_data, normal_data, cmplen, diff_limit) {
        return -1;
    }
    let (acc_left, acc_right) = (&scratch.left, &scratch.right);
    let mut min_diff = 100_000_000;
    for i in 1..cmplen {
        if acc_left[i - 1] + acc_right[cmplen - 1] > diff_limit {
            return -1;
        }
        let diff = acc_left[i - 1] + acc_right[i];
        if diff <= min_diff {
            min_diff = diff;
        }
    }
    min_diff
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The early-out fires as soon as the *ungapped* prefix alone exceeds the
    /// budget, and returns -1 even when a zero-cost split was already scored.
    /// Widen the budget past what the shifted prefix accumulates and the same
    /// input reports 0. Values traced against the C++ accumulators.
    #[test]
    fn early_out_discards_a_split_it_already_found() {
        let ins = b"ACGTXACGT";
        let norm = b"ACGTACGT";
        assert_eq!(diff_with_one_insertion(&mut MatchScratch::default(), ins, norm, norm.len(), 0), -1);
        assert_eq!(diff_with_one_insertion(&mut MatchScratch::default(), ins, norm, norm.len(), 1), -1);
        assert_eq!(diff_with_one_insertion(&mut MatchScratch::default(), ins, norm, norm.len(), 2), -1);
        assert_eq!(diff_with_one_insertion(&mut MatchScratch::default(), ins, norm, norm.len(), 3), 0);
    }

    /// When the ungapped prefix stays inside the budget, a cost comes back.
    #[test]
    fn reports_a_cost_when_the_prefix_stays_within_budget() {
        assert_eq!(diff_with_one_insertion(&mut MatchScratch::default(), b"AAAAAAAAA", b"AAAAAAAA", 8, 0), 0);
    }

    #[test]
    fn rejects_when_over_the_limit() {
        let mut sc = MatchScratch::default();
        assert_eq!(diff_with_one_insertion(&mut sc, b"AAAAAAAAA", b"TTTTTTTT", 8, 1), -1);
    }

    #[test]
    fn short_input_is_rejected_rather_than_panicking() {
        let mut sc = MatchScratch::default();
        assert_eq!(diff_with_one_insertion(&mut sc, b"AC", b"ACGT", 4, 1), -1);
        assert_eq!(diff_with_one_insertion(&mut sc, b"", b"", 0, 1), -1);
    }

    /// The scratch is reused across calls, so a later short call must not read
    /// values left over from an earlier long one.
    #[test]
    fn scratch_reuse_does_not_leak_between_calls() {
        let mut sc = MatchScratch::default();
        assert_eq!(diff_with_one_insertion(&mut sc, b"AAAAAAAAAAAAAAAAA", b"AAAAAAAAAAAAAAAA", 16, 0), 0);
        assert_eq!(diff_with_one_insertion(&mut sc, b"AAAAAAAAA", b"AAAAAAAA", 8, 0), 0);
        assert_eq!(diff_with_one_insertion(&mut sc, b"ACGTXACGT", b"ACGTACGT", 8, 0), -1);
    }
}
