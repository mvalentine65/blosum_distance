use super::*;

fn seq(n: usize, seed: u64) -> Vec<u8> {
    let mut s = seed;
    (0..n).map(|_| { s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); ((s >> 33) % 4) as u8 }).collect()
}

// a substitution at every `every`-th base
fn mutate(s: &[u8], every: usize) -> Vec<u8> {
    s.iter().enumerate().map(|(i, &c)| if i % every == every - 1 { (c + 1) % 4 } else { c }).collect()
}

#[test]
fn test_an_exact_copy_spans_the_unit() {
    let e = seq(150, 1);
    assert_eq!(banded(&e[20..110], &e, 20, 16), (1.0, 0, 90));
}

#[test]
fn test_substitutions_and_a_gap_lower_identity() {
    let e = seq(200, 2);
    let (ident, s, t) = banded(&mutate(&e[..105], 10), &e, 0, 16);
    assert_eq!((s, t), (0, 105));
    assert!((ident - 95.0 / 105.0).abs() < 1e-9);
    let u: Vec<u8> = e[..60].iter().chain(e[63..120].iter()).copied().collect();
    let (ident, s, t) = banded(&u, &e, 0, 16);
    assert_eq!(t - s, 117);
    assert!((ident - 117.0 / 120.0).abs() < 1e-9);
}

#[test]
fn test_each_unit_finds_the_exon_it_copies() {
    let exons = vec![seq(120, 10), seq(150, 11), seq(130, 12)];
    let units = vec![exons[1][10..100].to_vec(), mutate(&exons[2], 30), seq(140, 20)];
    let got = score_group(&exons, &units);
    assert_eq!((got[0].0, got[0].1, got[0].2), (1, 1.0, 1.0));
    assert!((got[0].3 - 90.0 / 150.0).abs() < 1e-9);
    assert_eq!(got[1].0, 2);
    assert!(got[1].1 >= 0.95 && got[1].2 >= 0.99);
    assert_eq!(got[2].0, -1);   // unrelated: too few shared k-mers to be anyone's copy
}

#[test]
fn test_a_diverged_copy_falls_short() {
    let exons = vec![seq(150, 11)];
    let got = score_group(&exons, &[mutate(&exons[0][10..100], 4)]);   // 75% identical
    assert!(got[0].0 == -1 || got[0].1 < MIN_IDENTITY);
}

#[test]
fn test_a_unit_much_longer_than_an_exon_is_no_copy_of_it() {
    let e = seq(100, 3);
    let u: Vec<u8> = e.iter().chain(seq(60, 4).iter()).copied().collect();
    assert_eq!(score_group(&[e], &[u])[0].0, -1);
}
