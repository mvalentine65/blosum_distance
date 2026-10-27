//! Unit tests for [`column_cull`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;

fn cull(records: &[(&str, &str)]) -> Vec<String> {
    let recs = records.iter().map(|(h, s)| (h.to_string(), s.to_string())).collect();
    let (out, ..) = cull_columns(recs, ".", 0.33, 10, 12, 25, true, 0.5, 10, 1.0, None, 0, None).unwrap();
    out.into_iter().map(|(h, _)| h).collect()
}

#[test]
fn test_min_length_counts_only_residues_in_reference_columns() {
    let gap = "-------------";
    let refs = ["FGDKESRSPHESAPQYAREI", "FWDKESRTPHESAPQYARQI", "FWDKESRSPHESAPQYARPI", "FWDKESRSPHEDAPQNARKI"];
    let tails = ["WEMAAAVAPVQATIRSVINI", "WHMAAAVAPHQALIRSVINI", "WEMAAAVAPHQCSIRSVWNI", "WEMAAVNAPHQATIRSVINI"];
    let rows: Vec<(String, String)> = (0..4).map(|i| (format!("g|ref{i}|x|."), format!("{}{gap}{}", refs[i], tails[i]))).collect();
    let mut recs: Vec<(&str, &str)> = rows.iter().map(|(h, s)| (h.as_str(), s.as_str())).collect();
    let full = format!("FWDKESRSPHESAPQYARKI{gap}WEMAAAVAPHQATIRSVINI");
    // 14 residues, 13 of them in insertion columns: one reference column only
    let frag = format!("{}IPYPHRYPLPVTRF{}", "-".repeat(19), "-".repeat(20));
    recs.push(("g|A|tax|1|1-2|+|c", &full));
    recs.push(("g|B|tax|2|1-2|+|c", &frag));
    let kept = cull(&recs);
    assert!(kept.iter().any(|h| h.starts_with("g|A|")));
    assert!(!kept.iter().any(|h| h.starts_with("g|B|")));
}

#[test]
fn test_splice_check_counts_residues_masked_left_of_the_block() {
    // One codon per residue; only V (GTT) and Q (CAG) contain a G.
    let codon = |a: char| match a {
        'A' => "ACT", 'L' => "CTT", 'K' => "AAA", 'P' => "CCA", 'T' => "ACA", 'N' => "AAC",
        'I' => "ATT", 'H' => "CAC", 'Y' => "TAC", 'F' => "TTC", 'S' => "TCT",
        'V' => "GTT", 'Q' => "CAG", _ => unreachable!(),
    };
    let refs: Vec<(String, String)> = (0..4)
        .map(|r| (format!("ref{r}."), format!("{}{}{}{}", if r == 0 { "AAAA" } else { "----" }, "L".repeat(30), "-".repeat(15), "S".repeat(30))))
        .collect();
    // 4 residues in edge-masked columns, then a stop-free GT...AG insertion.
    let cand = format!("KTNA{}VKPTNIHYFSLKPTQ{}", "LKPTNIHYFS".repeat(3), "SFYHINTPKL".repeat(3));
    let nt: String = cand.chars().map(codon).collect();
    let mut recs = refs.clone();
    recs.push(("cand".to_string(), cand));
    let nt_seqs: HashMap<String, String> = [("cand".to_string(), nt)].into_iter().collect();
    let (_, _, _, gff_intron, exons) =
        cull_columns(recs, ".", 0.33, 10, 12, 25, true, 0.5, 10, 1.0, Some(nt_seqs), 0, None).unwrap();
    assert_eq!(gff_intron.get("cand"), Some(&vec![(34, 15)]));
    assert_eq!(exons.get("cand"), Some(&vec![(0, 34), (49, 30)]));
}
