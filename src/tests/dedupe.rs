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

fn fastq(path: &Path, seqs: &[&str]) {
    let body: String = seqs.iter().enumerate()
        .map(|(i, s)| format!("@r{i}\n{s}\n+\n{}\n", "I".repeat(s.len())))
        .collect();
    std::fs::write(path, body).unwrap();
}

fn items(paths: &[PathBuf]) -> Result<Vec<Item>> {
    let mut stream = ItemStream::new(paths);
    let mut out = Vec::new();
    while let Some(it) = stream.next_item()? {
        out.push(it);
    }
    Ok(out)
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

#[test]
fn mate_files_of_unequal_length_are_an_error() {
    let d = dir("unequal");
    let (r1, r2) = (d.join("lib_1.fq"), d.join("lib_2.fq"));
    fastq(&r1, &["ACGTACGTAC", "CCGTACGTAC", "GCGTACGTAC"]);
    fastq(&r2, &["TTGTACGTAC", "TAGTACGTAC"]);
    let paths = vec![r1.clone(), r2.clone()];
    let err = items(&paths).err().expect("the stream must not stop quietly");
    assert!(err.to_string().contains("different numbers of records"), "{err}");
    // the parallel reader, both when the short file ends inside a batch and on a batch edge
    for batch in [2, 4] {
        let res = read_pair_parallel(Box::new(RecordSource::open(&r1).unwrap()),
                                     Box::new(RecordSource::open(&r2).unwrap()), batch, &mut |_| true);
        assert!(res.is_err(), "batch {batch}");
    }
}
