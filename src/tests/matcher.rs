//! Unit tests for [`matcher`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

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
