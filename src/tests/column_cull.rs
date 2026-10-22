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
