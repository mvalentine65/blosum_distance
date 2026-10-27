//! Unit tests for [`dedupe`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;

/// A fresh scratch directory for one test.
fn dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("bd_dedupe_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}


#[test]
fn a_wrapped_fasta_record_reads_as_one_sequence() {
    let d = dir("wrapped");
    let seq: String = (0..200).map(|i| b"ACGT"[(i * 7 + i / 3) % 4] as char).collect();
    for (name, width, nl) in [("lf.fa", 60, "\n"), ("crlf.fa", 80, "\r\n")] {
        let lines: Vec<&str> = seq.as_bytes().chunks(width).map(|c| std::str::from_utf8(c).unwrap()).collect();
        let path = d.join(name);
        std::fs::write(&path, format!(">r{nl}{}{nl}", lines.join(nl))).unwrap();
        let mut src = RecordSource::open(&path).unwrap();
        assert_eq!(src.next_record().unwrap().unwrap().0, seq.as_bytes(), "{name}");
        assert!(src.next_record().unwrap().is_none());
    }
}
