//! AVX2 Forward/Backward against the generic log-space path, on models whose
//! best paths run long deletions across many lanes.

use super::*;
use crate::genome_splice::align::{align_both, residue_code, Dp, Profile};
use crate::genome_splice::hmm::{read_hmms, Hmm};

const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// A profile of `pep` (its residue at 0.8 per node) with sticky deletions.
fn model(pep: &[u8], tag: &str) -> Hmm {
    let nl = |p: f64| format!("{:.5}", -p.ln());
    let row = |r: Option<u8>| AA.iter().map(|&a| nl(match r { Some(x) if x == a => 0.8, Some(_) => 0.2 / 19.0, None => 0.05 }))
        .collect::<Vec<_>>().join(" ");
    let tr = |last: bool| [nl(0.6), nl(0.05), if last { "*".into() } else { nl(0.35) }, nl(0.5), nl(0.5), nl(0.1), nl(0.9)].join(" ");
    let mut s = format!("HMMER3/f\nNAME  t\nLENG  {}\nALPH  amino\nHMM          A\n            m->m\n  COMPO   {}\n          {}\n          {}\n",
                        pep.len(), row(None), row(None), tr(false));
    for (k, &r) in pep.iter().enumerate() {
        s += &format!("{:7} {} {}\n          {}\n          {}\n", k + 1, row(Some(r)), k + 1, row(None), tr(k + 1 == pep.len()));
    }
    s += "//\n";
    let path = std::env::temp_dir().join(format!("avx_test_{}_{tag}.hmm", std::process::id()));
    std::fs::write(&path, s).unwrap();
    let h = read_hmms(path.to_str().unwrap(), 1).unwrap().remove(0);
    let _ = std::fs::remove_file(&path);
    h
}

#[test]
fn avx_forward_backward_match_generic_across_long_deletions() {
    if !available() { return; }
    let mut x = 7u64;
    let mut rnd = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
    let (mut worst, mut worst_at) = (0.0f32, String::new());
    for m in [12usize, 16, 17, 24, 25, 32, 40, 48] {
        let pep: Vec<u8> = (0..m).map(|_| AA[(rnd() % 20) as usize]).collect();
        let hmm = model(&pep, &m.to_string());
        let p = Profile::new(&hmm, 1, m);
        let om = OProfile::new(&p);
        for del in [7usize, 9, 12] {
            for a in [1usize, 2, 3] {
                if a + del >= m { continue; }
                let pepseq: Vec<u8> = pep[..a].iter().chain(&pep[a + del..]).copied().collect();
                let seq: Vec<u8> = pepseq.iter().map(|&c| residue_code(c)).collect();
                let (g, _) = align_both(&p, &seq, &mut Dp::default());
                let mut dp = AvxDp::default();
                let f = fwd_bits(&om, &seq, &mut dp);
                let _ = align(&p, &om, &seq, &mut dp);
                for (what, v) in [("fwd", f), ("bck", dp.bck_bits)] {
                    let d = (v - g.bits).abs();
                    if d > worst { worst = d; worst_at = format!("M={m} del={del} a={a} {what}: avx {v:.3} generic {:.3}", g.bits); }
                }
            }
        }
    }
    assert!(worst < 0.05, "worst gap {worst:.3} bits at {worst_at}");
}
