use std::fs::File;
use std::io::Write;
use pyo3::prelude::*;
use flate2;
use flate2::Compression;
use flate2::write::GzEncoder;

#[pyfunction]
pub fn writeFastaCompressed(path: String, records: Vec<(&str, &str)>) {
    let output_path = match &path[path.len()-3..] {
        ".gz" => path,
        _ => path + ".gz"
    };

    let file = File::create(output_path).unwrap();
    let mut encoder = GzEncoder::new(file, Compression::default());
    let mut line;
    let mut buf;
    for (header, sequence) in records {
        line = format!(">{}\n{}\n", header, sequence);
        buf = line.as_bytes();
        encoder.write_all(buf).unwrap();
    }
    encoder.finish().unwrap();
}

#[pyfunction]
pub fn writeFastaUncompressed(path: String, records: Vec<(&str, &str)>) {
    let mut file = File::create(path).unwrap();
    let mut line;
    let mut buf;
    for (header, sequence) in records {
        line = format!(">{}\n;{}\n", header, sequence);
        buf = line.as_bytes();
        file.write_all(buf).unwrap();
    }
    file.flush().unwrap();
}