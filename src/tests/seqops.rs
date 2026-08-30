//! Unit tests for [`seqops`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;

#[test]
fn revcomp_matches_fastp_table() {
    assert_eq!(reverse_complement(b"ATCGN"), b"NCGAT".to_vec());
    // lowercase shares the nibble, unknown bases fall to N
    assert_eq!(reverse_complement(b"atcg"), b"CGAT".to_vec());
    assert_eq!(reverse_complement(b"AXA"), b"TNT".to_vec());
}

#[test]
fn mismatch_counters() {
    assert_eq!(count_mismatches(b"AAAA", b"AATA", 4), 1);
    assert!(count_mismatches_bounded(b"AAAA", b"TTTT", 4, 1) > 1);
    assert!(count_mismatches_bounded(b"AAAA", b"TTTT", 4, 0) > 0);
    assert_eq!(count_mismatches_bounded(b"AAAA", b"AAAA", 4, 0), 0);
}

/// The vector path must agree with the scalar one everywhere the result is
/// load-bearing: exact totals, and exact counts up to the limit. Lengths
/// span the 16-byte block boundary in both directions.
#[test]
fn simd_agrees_with_scalar() {
    let mut seed = 0x12345678u64;
    let mut rand = move || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        (seed >> 33) as usize
    };
    let bases = b"ACGT";
    for len in 0..80usize {
        for _ in 0..20 {
            let a: Vec<u8> = (0..len).map(|_| bases[rand() % 4]).collect();
            let mut b = a.clone();
            for _ in 0..rand() % 6 {
                if len > 0 {
                    let i = rand() % len;
                    b[i] = bases[rand() % 4];
                }
            }
            let want = count_mismatches_scalar(&a, &b);
            assert_eq!(count_mismatches(&a, &b, len), want, "len {len}");
            for limit in 0..6usize {
                let got = count_mismatches_bounded(&a, &b, len, limit);
                let want_b = count_mismatches_bounded_scalar(&a, &b, limit);
                if want_b <= limit {
                    assert_eq!(got, want_b, "exact under limit, len {len} limit {limit}");
                } else {
                    assert!(got > limit, "must still reject, len {len} limit {limit}");
                }
            }
        }
    }
}

#[test]
fn adjacent_diffs_is_the_complexity_measure() {
    assert_eq!(count_adjacent_diffs(b"AAAA"), 0);
    assert_eq!(count_adjacent_diffs(b"ATAT"), 3);
    assert_eq!(count_adjacent_diffs(b""), 0);
    assert_eq!(count_adjacent_diffs(b"A"), 0);
}

#[test]
fn quality_metrics_sum_phred_not_bytes() {
    // 'I' is phred 40, '#' is phred 2
    let m = count_quality_metrics(b"II#I", b"ACNG", i64::from(b'5'));
    assert_eq!(m.total_qual, 40 + 40 + 2 + 40);
    assert_eq!(m.low_qual, 1);
    assert_eq!(m.n_bases, 1);
}
