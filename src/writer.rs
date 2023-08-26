use std::fs::File;
use std::io::Write;
use pyo3::prelude::*;
use flate2;
use flate2::Compression;
use flate2::write::GzEncoder;
use gzip_header;
use std::time::{SystemTime, UNIX_EPOCH};
use std::path::Path;
use crc32fast;

#[pyfunction]
pub fn writeFastaCompressed(path: String, records: Vec<(&str, &str)>) {
    let output_path = match &path[path.len()-3..] {
        ".gz" => path,
        _ => path + ".gz"
    };
    let mut builder = gzip_header::GzBuilder::new();
    // let current_time = SystemTime::now();
    // let duration = current_time.duration_since(UNIX_EPOCH).expect("Time went backwards");
    // let time_as_u32 = duration.as_secs() as u32;

    builder = builder.filename(Path::new(&output_path).file_name().unwrap().to_str().unwrap().as_bytes().to_vec());
    // builder = builder.mtime(time_as_u32);
    let mut file = File::create(output_path).unwrap();
    // file.write_all(&builder.into_header()).unwrap();

    // let mut crc32 = crc32fast::Hasher::new();
    let mut encoder = GzEncoder::new(&file, Compression::best());
    let mut line;
    let mut buf;
    for (header, sequence) in records {
        line = format!(">{}\n{}\n", header, sequence);
        buf = line.as_bytes();
        // crc32.update(&buf);
        encoder.write_all(buf).unwrap();
    }
    // let checksum = crc32.finalize();
    // let checksum_bytes: [u8; 4] = [
    //     (checksum >> 24) as u8,
    //     (checksum >> 16) as u8,
    //     (checksum >> 8) as u8,
    //     checksum as u8,
    // ];
    encoder.flush().unwrap();
    // encoder.write_all(&checksum_bytes).unwrap();
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