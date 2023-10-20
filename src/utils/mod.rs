use std::fs::File;
use std::io::Write;
use pyo3::prelude::*;
// use flate2;
// use flate2::Compression;
// use flate2::write::GzEncoder;
use gzip_header;
// use std::time::{SystemTime, UNIX_EPOCH};
use std::path::Path;
use crc32fast;
use fastx::FastX;
use crc;
use deflate::Compression;
use deflate::write::GzEncoder;
#[pyfunction]
pub fn write_fasta_compressed(path: String, records: Vec<(String, String)>) {
    let output_path = match &path[path.len()-3..] {
        ".gz" =>path,
        _ =>  path + ".gz"
    };
    let mut file = File::create(&output_path).unwrap();
    let mut encoder= GzEncoder::new(file, Compression::default());
    let mut line;
    let mut buf;
    for (header, sequence) in records {
        line = format!(">{}\n{}\n", header, sequence);
        buf = line.as_bytes();
        encoder.write_all(buf).unwrap();
    }
    encoder.flush().unwrap();
    encoder.finish().unwrap();
}

#[pyfunction]
pub fn write_fasta_uncompressed(path: String, records: Vec<(String, String)>) {
    let mut file = File::create(path).unwrap();
    let mut line;
    let mut buf;
    for (header, sequence) in records {
        line = format!(">{}\n{}\n", header, sequence);
        buf = line.as_bytes();
        file.write_all(buf).unwrap();
    }
    file.flush().unwrap();
}

fn parse_fasta(path: &str) -> Vec<(String, String)> {
    let mut records = Vec::<(String, String)>::new();
    let mut reader = FastX::reader_from_path(Path::new(path)).unwrap();
    let mut fastx_record = FastX::from_reader(&mut reader).unwrap();
    while let Ok(_some @ 1..=usize::MAX) = fastx_record.read(&mut reader) {
        records.push((fastx_record.id().to_string(), String::from_utf8(fastx_record.seq()).unwrap()));
    }
    records
}

