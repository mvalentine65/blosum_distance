//! Unit tests for [`periodicity`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;

fn rep(unit: &[u8], len: usize) -> Vec<u8> {
    (0..len).map(|i| unit[i % unit.len()]).collect()
}

#[test]
fn catches_a_homopolymer() {
    assert_eq!(repeat_period(&rep(b"A", 128), 0.8), Some(1));
}

#[test]
fn catches_di_and_trinucleotide_repeats() {
    assert_eq!(repeat_period(&rep(b"AG", 128), 0.8), Some(2));
    assert_eq!(repeat_period(&rep(b"CAG", 128), 0.8), Some(3));
}

/// The case entropy misses: four distinct bases keep the diversity up.
#[test]
fn catches_the_satellite_entropy_lets_through() {
    assert_eq!(repeat_period(&rep(b"CATA", 128), 0.8), Some(4));
    assert_eq!(repeat_period(&rep(b"GACAAAT", 128), 0.8), Some(7));
}

/// The point of normalising by composition: AT-biased but aperiodic.
#[test]
fn leaves_at_rich_aperiodic_sequence_alone() {
    // 80% A/T, no periodicity -- roughly weevil mitochondrial composition.
    let mut seed = 0x5EEDu64;
    let mut rand = move || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        (seed >> 33) as usize
    };
    let biased = b"AATTAATTAC"; // draw pool: 80% A/T
    let s: Vec<u8> = (0..200).map(|_| biased[rand() % biased.len()]).collect();
    assert_eq!(repeat_period(&s, 0.8), None, "{}", String::from_utf8_lossy(&s));
}

#[test]
fn leaves_ordinary_sequence_alone() {
    let s = b"CATAACAAAAGGCAGAGACCTGATGAAGGATGATGGTGTGATGAATTCCTGGGACGAAGCATGGCATAAG";
    assert!(repeat_period(s, 0.8).is_none());
}

/// An imperfect repeat still counts: real satellites carry miscalls. Each
/// substitution breaks two lag-2 comparisons, so 1 in 25 scores 0.825.
#[test]
fn tolerates_an_imperfect_repeat() {
    let mut s = rep(b"AG", 128);
    for i in (0..s.len()).step_by(25) {
        s[i] = b'C';
    }
    assert_eq!(repeat_period(&s, 0.8), Some(2));
}

/// The other side of that boundary: 1 in 17 scores 0.772 and is let go.
#[test]
fn a_heavily_degraded_repeat_is_let_through() {
    let mut s = rep(b"AG", 128);
    for i in (0..s.len()).step_by(17) {
        s[i] = b'C';
    }
    assert_eq!(repeat_period(&s, 0.8), None);
}

#[test]
fn reports_the_shortest_period_not_a_multiple() {
    // (AG)n matches at 2, 4, 6, ... -- the unit length is the useful answer.
    assert_eq!(repeat_period(&rep(b"AG", 128), 0.8), Some(2));
}

#[test]
fn short_input_is_safe() {
    assert!(repeat_period(b"", 0.8).is_none());
    assert!(repeat_period(b"AC", 0.8).is_none());
    assert!(repeat_period(b"ACG", 0.8).is_none());
}

#[test]
fn ns_do_not_crash_the_composition_estimate() {
    let s = vec![b'N'; 128];
    assert!(repeat_period(&s, 0.8).is_none());
}

#[test]
fn threshold_moves_the_boundary() {
    let mut s = rep(b"AG", 128);
    // Degrade it until it sits between the two thresholds.
    for i in (0..s.len()).step_by(5) {
        s[i] = b'T';
    }
    let loose = repeat_period(&s, 0.6).is_some();
    let strict = repeat_period(&s, 0.95).is_some();
    assert!(loose || !strict, "a looser threshold must never flag less");
}

/// Every frame and both strands of the same satellite must give one key.
#[test]
fn rotations_and_strands_fold_together() {
    let k = canonical_unit(b"AG", 2);
    for v in [&b"AG"[..], b"GA", b"CT", b"TC"] {
        assert_eq!(canonical_unit(v, 2), k, "{}", String::from_utf8_lossy(v));
    }
    assert_eq!(String::from_utf8_lossy(&k), "AG");
}

/// A read starting mid-unit is still filed under the right satellite.
#[test]
fn a_read_starting_mid_unit_is_not_mislabelled() {
    let mut s = b"AA".to_vec();
    s.extend(rep(b"GA", 130));
    assert_eq!(repeat_period(&s, 0.8).unwrap(), 2);
    assert_eq!(canonical_unit(&s, 2), b"AG".to_vec());
}

/// Caught at a multiple, the tally must still name the real satellite.
#[test]
fn a_unit_caught_at_a_multiple_reduces_to_its_primitive() {
    let s = rep(b"AG", 132);
    for p in [2usize, 4, 6, 12] {
        assert_eq!(canonical_unit(&s, p), b"AG".to_vec(), "period {p}");
    }
    assert_eq!(canonical_unit(&rep(b"A", 132), 12), b"A".to_vec());
}

#[test]
fn distinct_units_stay_distinct() {
    assert_ne!(canonical_unit(b"AAT", 3), canonical_unit(b"ACT", 3));
}

#[test]
fn unit_is_taken_at_the_period_not_the_whole_read() {
    assert_eq!(canonical_unit(b"CATCATCATCAT", 3).len(), 3);
}

#[test]
fn canonical_unit_is_safe_on_edges() {
    assert!(canonical_unit(b"", 3).is_empty());
    assert!(canonical_unit(b"AC", 0).is_empty());
    assert_eq!(canonical_unit(b"AC", 9).len(), 2);
}
