//! Unit tests for [`genome_splice::run`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;

fn ex(start: i64, end: i64, k1: i64, k2: i64, src: u8) -> ChainExon {
    ChainExon { start, end, k1, k2, src, codon: 0, bits: 0.0, run: false, owes: false }
}

fn row(gene: &str, ex: Vec<ChainExon>) -> Chain {
    Chain { gene: gene.into(), model: "AA1".into(), scaffold: "NT_033778.4".into(), strand: b'+', klo: 1, khi: 525,
            flank_only: false, passive: false, imx: false, alt_of: None, rebased: None, ex }
}

// Drosophila straw: exon 43037, then X1 (A/F's block) or X2 (D's and G's), then the ends
fn straw() -> Vec<Chain> {
    let e37 = ex(5439570, 5439754, 240, 293, SRC_INPUT);
    let x1 = ex(5440610, 5440805, 297, 345, SRC_ORF);
    let x2 = ex(5453281, 5453476, 297, 345, SRC_SIBFILL);
    vec![
        row("AA1#2", vec![e37, x1, ex(5440858, 5441070, 353, 433, SRC_INPUT), ex(5446535, 5446765, 438, 499, SRC_INPUT)]),
        row("AA1#2_1", vec![e37, x2, ex(5453559, 5454151, 354, 521, SRC_INPUT)]),
        row("AA1#2_2", vec![e37, x1, ex(5454421, 5455013, 357, 521, SRC_INPUT)]),
    ]
}

#[test]
fn test_an_end_takes_the_copy_beside_it() {
    let mut cs = straw();
    let log = nearest_copies(&mut cs, 0.75);
    // 43049's row gives up X1 (13.6 kb away) for X2 (945 nt away): G, not E
    assert_eq!((cs[2].ex[1].start, cs[2].ex[1].end, cs[2].ex[1].src), (5453281, 5453476, SRC_SIBFILL));
    assert_eq!(log.lines().count(), 1);
    // A/F keeps X1 (X2 lies outside its gap) and D keeps X2 (nearer 43048 than X1)
    assert_eq!(cs[0].ex[1].start, 5440610);
    assert_eq!(cs[1].ex[1].start, 5453281);
}

#[test]
fn test_a_start_keeps_its_local_copy() {
    // P450 alternative start: start B with its own exon 2, the sibling's exon 2 further on
    let body = ex(22639189, 22639279, 92, 125, SRC_INPUT);
    let mut cs = vec![
        row("clan4#30", vec![ex(22637013, 22637169, 5, 53, SRC_INPUT), ex(22638408, 22638517, 57, 93, SRC_SPLICE), body]),
        row("clan4#30_1", vec![ex(22635270, 22635426, 8, 52, SRC_INPUT), ex(22636271, 22636380, 57, 93, SRC_SIBFILL), body]),
    ];
    assert_eq!(nearest_copies(&mut cs, 0.75), "");
    assert_eq!(cs[1].ex[1].start, 22636271);
}

// A row's recovered exon e (1400-1600) whose own flank is b; the sibling shares a and offers y
fn swap_case(y: ChainExon) -> Vec<Chain> {
    let a = ex(1000, 1100, 1, 30, SRC_INPUT);
    vec![
        row("g#1", vec![a, y, ex(3000, 3100, 91, 120, SRC_INPUT)]),
        row("g#1_1", vec![a, ex(1400, 1600, 31, 90, SRC_ORF), ex(1800, 1900, 91, 120, SRC_INPUT)]),
    ]
}

#[test]
fn test_a_copy_overlapping_the_exon_is_not_swapped_in() {
    let mut cs = swap_case(ex(1500, 1700, 31, 90, SRC_ORF));
    assert_eq!(nearest_copies(&mut cs, 0.75), "");
    assert_eq!(cs[1].ex[1].start, 1400);
}

#[test]
fn test_a_copy_of_unlike_length_is_not_swapped_in() {
    // 101 nt for 201: under 0.75
    let mut cs = swap_case(ex(1650, 1750, 31, 90, SRC_ORF));
    assert_eq!(nearest_copies(&mut cs, 0.75), "");
    // 151 nt for 201: just over 0.75, and nearer b
    let mut cs = swap_case(ex(1610, 1760, 31, 90, SRC_ORF));
    assert_eq!(nearest_copies(&mut cs, 0.75).lines().count(), 1);
    assert_eq!((cs[1].ex[1].start, cs[1].ex[1].src), (1610, SRC_SIBFILL));
}

#[test]
fn test_keep_low_breaks_a_tie_by_coding_order() {
    // 1001-1300 cut at 1151-1150: two 150 nt pieces; the coding-first piece stays on either strand
    assert!(keep_low(b'+', 1001, 1300, 1151, 1150));
    assert!(!keep_low(b'-', 1001, 1300, 1151, 1150));
    // otherwise the longer piece stays, whatever the strand
    for st in [b'+', b'-'] {
        assert!(keep_low(st, 1001, 1300, 1154, 1153));
        assert!(!keep_low(st, 1001, 1300, 1148, 1147));
    }
}

fn jx(gap: i64, dis: Vec<(usize, i64, i64, bool)>) -> JxOut {
    JxOut { sites: JxSites::default(), pseudo: Pseudo::default(), status: "ok", gap, dis, stopfree: f64::NAN }
}

#[test]
fn test_a_stop_with_no_room_for_an_intron_is_read_past() {
    let stop = (DIS_STOP as usize, 18416355, 0, true);
    assert!(stop_only(&jx(9, vec![stop])));
    assert!(!stop_only(&jx(9, vec![(DIS_STOP as usize, 18416355, 0, false)])));   // not a strict stop
    assert!(!stop_only(&jx(9, vec![stop, (DIS_FS as usize, 18416300, 1, true)])));  // a frameshift too
    assert!(!stop_only(&jx(MIN_INTRON as i64, vec![stop])));                      // an intron fits
    assert!(!stop_only(&jx(3, Vec::new())));
}

#[test]
fn test_only_a_short_last_piece_goes_past_a_stop() {
    assert!(past_stop(117, 3054));    // the elephant shark Nu tail, 3.8% of its row
    assert!(past_stop(100, 1000));
    assert!(!past_stop(101, 1000));
}

#[test]
fn test_a_stop_needs_a_matched_tail_to_be_real() {
    assert!(real_stop(3, 3) && real_stop(10, 8));
    assert!(!real_stop(2, 2));     // too few residues past it: a hit's edge
    assert!(!real_stop(10, 7));    // past it the hit leaves the model
}

/// A profile of `pep` (its residue at 0.8 per node), as the ORF tests build it.
fn straw_model(pep: &[u8]) -> Hmm {
    const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
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
    let path = std::env::temp_dir().join(format!("stop_in_hit_test_{}.hmm", std::process::id()));
    std::fs::write(&path, s).unwrap();
    let h = super::super::hmm::read_hmms(path.to_str().unwrap(), 1).unwrap().remove(0);
    let _ = std::fs::remove_file(&path);
    h
}

fn codons(p: &[u8]) -> Vec<u8> {
    p.iter().flat_map(|&a| match a {
        b'A' => *b"GCT", b'C' => *b"TGT", b'D' => *b"GAT", b'E' => *b"GAA", b'F' => *b"TTT", b'G' => *b"GGT",
        b'H' => *b"CAT", b'I' => *b"ATT", b'K' => *b"AAA", b'L' => *b"CTT", b'M' => *b"ATG", b'N' => *b"AAT",
        b'P' => *b"CCT", b'Q' => *b"CAA", b'R' => *b"CGT", b'S' => *b"TCT", b'T' => *b"ACT", b'V' => *b"GTT",
        b'W' => *b"TGG", _ => *b"TAT",
    }).collect()
}

#[test]
fn test_a_stop_inside_a_hit_is_kept_only_with_model_residues_past_it() {
    const PEP: &[u8] = b"VMNPLDAQCNRQRNDAICVSQLKNALEIDRGILAEKPDVKIFLPFRFFVYRAEDLFQPNTYNRFL";
    let hmm = straw_model(PEP);
    // hit a: residues 1-50 with residue 31 a TAA stop; then 300 nt of spacer and exon b
    let mut sc = codons(&PEP[..30]);
    sc.extend(b"TAA");
    sc.extend(codons(&PEP[31..50]));
    sc.extend(std::iter::repeat_n(b'C', 300));
    sc.extend(codons(&PEP[50..]));
    let ex = |start: i64, end: i64, k1: i64, k2: i64, src: u8| ChainExon { start, end, k1, k2, src, codon: start, bits: 0.0, run: false, owes: false };
    let a = ex(1, 150, 1, 50, SRC_INPUT);
    let b = ex(451, 495, 51, 65, SRC_ORF);
    let c = Chain { gene: "straw#0".into(), model: "straw".into(), scaffold: "s".into(), strand: b'+', klo: 1, khi: 65,
                    flank_only: false, passive: false, imx: false, alt_of: None, rebased: None, ex: vec![a, b] };
    let jx_at = |g: i64| Junction { status: JxStatus::Ok, s: Vec::new(), w: sc.len() as i64, lo: 0, strand: b'+', gap: 0, d: 0, axe: 0, bxs: 0,
        res: super::super::splice::SpliceResult { dis: vec![super::super::splice::Disable { kind: DIS_STOP, pos: g - 1, n: 0, tight: true }],
                                                 ..Default::default() } };
    let mut al = Aligner::default();
    // 19 residues past the stop, all on the model: the gene's own stop
    assert!(stop_in_hit(&mut al, 0, &hmm, code(1).unwrap(), &c, &a, &b, &jx_at(91), 1, &sc));
    // the same stop two residues before the hit's edge: a retained piece of intron
    let short = ex(1, 99, 1, 33, SRC_INPUT);
    let c2 = Chain { ex: vec![short, b], ..c.clone() };
    assert!(!stop_in_hit(&mut al, 0, &hmm, code(1).unwrap(), &c2, &short, &b, &jx_at(91), 1, &sc));
    // a stop in an exonfill piece is never kept this way
    let orf = ex(1, 150, 1, 50, SRC_ORF);
    let c3 = Chain { ex: vec![orf, b], ..c.clone() };
    assert!(!stop_in_hit(&mut al, 0, &hmm, code(1).unwrap(), &c3, &orf, &b, &jx_at(91), 1, &sc));
}
