//! Unit tests for [`genome_splice::orf`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;
use crate::genome_splice::hmm::read_hmms;
use crate::genome_splice::sites::code;

// straw (Drosophila Laccase2) exon of isoforms A/F/E
const PEP: &[u8] = b"VMNPLDAQCNRQRNDAICVSQLKNALEIDRGILAEKPDVKIFLPFRFFVYRAEDLFQPNTYNRFL";
const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// One codon per residue, so a peptide reads back as itself in frame 0.
fn back(p: &[u8]) -> Vec<u8> {
    p.iter().flat_map(|&a| match a {
        b'A' => *b"GCT", b'C' => *b"TGT", b'D' => *b"GAT", b'E' => *b"GAA", b'F' => *b"TTT", b'G' => *b"GGT",
        b'H' => *b"CAT", b'I' => *b"ATT", b'K' => *b"AAA", b'L' => *b"CTT", b'M' => *b"ATG", b'N' => *b"AAT",
        b'P' => *b"CCT", b'Q' => *b"CAA", b'R' => *b"CGT", b'S' => *b"TCT", b'T' => *b"ACT", b'V' => *b"GTT",
        b'W' => *b"TGG", _ => *b"TAT",
    }).collect()
}

/// A profile of `pep` in HMMER3 text (its residue at 0.8 per node), read back.
fn model(pep: &[u8]) -> Hmm {
    let nl = |p: f64| format!("{:.5}", -p.ln());
    let row = |r: Option<u8>| AA.iter().map(|&a| nl(match r { Some(x) if x == a => 0.8, Some(_) => 0.2 / 19.0, None => 0.05 }))
        .collect::<Vec<_>>().join(" ");
    let tr = |last: bool| [nl(0.9), nl(0.05), if last { "*".into() } else { nl(0.05) }, nl(0.5), nl(0.5), nl(0.5), nl(0.5)].join(" ");
    let mut s = format!("HMMER3/f\nNAME  straw\nLENG  {}\nALPH  amino\nHMM          A\n            m->m\n  COMPO   {}\n          {}\n          {}\n",
                        pep.len(), row(None), row(None), tr(false));
    for (k, &r) in pep.iter().enumerate() {
        s += &format!("{:7} {} {}\n          {}\n          {}\n", k + 1, row(Some(r)), k + 1, row(None), tr(k + 1 == pep.len()));
    }
    s += "//\n";
    let path = std::env::temp_dir().join(format!("orf_cache_test_{}.hmm", std::process::id()));
    std::fs::write(&path, s).unwrap();
    let h = read_hmms(path.to_str().unwrap(), 1).unwrap().remove(0);
    let _ = std::fs::remove_file(&path);
    h
}

/// orf_window over `w` on `strand`: its rows and the exons it keeps.
fn search(al: &mut Aligner, hmm: &Hmm, w: &[u8], strand: u8, cache: Option<&OrfCache>) -> (String, String) {
    let o = OrfOpts { margin: 5, min_aa: 10, thr: 1.0, near: 100, revthr: 3.0, minm: 12, all: true, decoy: false,
                      overlap: 0, anchored: false, rank: false, runs: true };
    let (mut out, mut kept) = (String::new(), Vec::new());
    orf_window(al, 0, hmm, code(1).unwrap(), &o, "w0", "g\ts\t+", 1000, strand, 1, hmm.m as i64, 1, w.len() as i64, w,
               0.5, false, false, &mut out, &mut kept, cache);
    (out, format!("{kept:?}"))
}

#[test]
fn test_orf_window_reads_the_same_from_its_cache() {
    let hmm = model(PEP);
    // two copies of the exon in their own ORFs among stop-broken spacer, so there are several ORFs per frame
    let mut w = b"TAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_vec();
    w.extend(back(PEP));
    w.extend(b"TAAGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTTAG");
    w.extend(back(&PEP[5..]));
    w.extend(b"TGACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    let rc: Vec<u8> = w.iter().rev().map(|&b| match b { b'A' => b'T', b'T' => b'A', b'C' => b'G', _ => b'C' }).collect();
    let mut al = Aligner::default();
    for (seq, strand) in [(&w, b'+'), (&rc, b'-')] {
        let plain = search(&mut al, &hmm, seq, strand, None);
        let cache = OrfCache::default();
        assert_eq!(search(&mut al, &hmm, seq, strand, Some(&cache)), plain);
        // a window may repeat a peptide (poly-C spacer reads the same in two frames), so a first pass can hit too
        let (hits, misses) = (cache.hits(), cache.misses());
        assert!(misses > 0);
        // the same window again reads every ORF back from the cache
        assert_eq!(search(&mut al, &hmm, seq, strand, Some(&cache)), plain);
        assert_eq!((cache.hits(), cache.misses()), (2 * hits + misses, misses));
        assert!(plain.0.lines().count() >= 2, "{}", plain.0);
    }
}
