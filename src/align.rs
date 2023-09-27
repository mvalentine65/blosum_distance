use std::collections::{HashMap, HashSet};
use std::{fs, process};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use tempfile::NamedTempFile;
use pyo3::prelude::*;
fn has_multiple_sequences(path: &Path) -> bool {
    let file = fs::File::open(&path).expect(&format!("Failed to open file {}", path.to_str().unwrap()));
    let reader = BufReader::new(file);
    let mut count = 0_u8;
    for _line in reader.lines() {
        let line = _line.expect(&format!("Cannot parse fasta line in file {}", path.to_str().unwrap()));
        if line.starts_with('>') {
            count += 1;
            if count >= 2 {
                return true;
            }
        }
    }
    false
}

#[pyfunction]
pub fn process_cluster_file(p2: &str) {
    let p2_path = Path::new(p2);

    if !p2_path.exists() {
        println!("Cluster file '{}' does not exist. Skipping...", p2);
        return;
    }
    if !has_multiple_sequences(p2_path) {
        return;
    }
    let p2_number = p2.trim_start_matches("unaligned_cluster_");
    let intermediate_file = NamedTempFile::new().expect("cannot create temp file");
    let temp_path = intermediate_file.path().to_str().expect("cannot get path to temp file");
    let _status = std::process::Command::new("clustalo")
        .args(&["-i", p2, "-o", temp_path, "--threads=1", "--full", "--force"])
        .status()
        .expect("Failed to execute clustalo");

    let out_file = format!("output_{}.fasta", p2_number);
    let _status = std::process::Command::new("clustalo")
        .args(&["--p1", "references", "--p2", temp_path, "-o", &out_file, "--threads=1", "--full", "--is-profile", "--force"])
        .status()
        .expect("Failed to execute clustalo");
}

pub fn find_kmers(fasta: &HashMap<String, String>, kmer_length: usize) -> HashMap<String, HashSet<String>> {
    let mut kmers:HashMap<String, HashSet<String>> = HashMap::new();
    let mut kmer: &str;
    for (header, sequence) in fasta.iter() {
        for i in 0..(sequence.len() - &kmer_length) {
            kmer = &sequence[i..i + &kmer_length];
            if !kmer.contains('*') && !kmer.contains('-') {
                kmers.entry(header.to_owned())
                    .or_insert(HashSet::<String>::new())
                    .insert(kmer.to_string());

            }
        }
    }
    kmers
}

#[pyfunction]
pub fn generate_clusters(data: HashMap<String, String>, kmer_length: usize, kmer_percent:f32) -> Vec<Vec<String>> {
    let mut cluster_children: HashMap<String, Vec<String>> = data.keys().map(|header| (header.clone(), vec![header.clone()])).collect();
    let mut kmers = find_kmers(&data, kmer_length);
    let gene_headers: Vec<String> = kmers.keys().cloned().collect();
    let empty_set = HashSet::<String>::new();
    let mut child_sets: HashMap<String, HashSet<String>> = data.keys().map(|header| (header.clone(), HashSet::new())).collect();

    for iteration in 0..2 {
        let mut processed_headers = HashSet::new();
        let mut merge_occurred = true;

        while merge_occurred {
            merge_occurred = false;

            for i in (0..gene_headers.len()).rev() {
                let master_header = &gene_headers[i];
                let master = kmers.get(master_header).unwrap_or(&empty_set);

                if !master.is_empty() {
                    for j in 0..i {
                        let candidate_header = &gene_headers[j];

                        if processed_headers.contains(candidate_header) {
                            continue;
                        }

                        let candidate = kmers.get(candidate_header).unwrap_or(&empty_set);
                        let mut matched = false;

                        let mut headers_to_check: Vec<&str> = vec![master_header];
                        headers_to_check.extend(child_sets[master_header].iter().map(|s| s.as_str()));

                        for header_to_check in headers_to_check {
                            let set_to_check = kmers.get(header_to_check).unwrap_or(&empty_set);
                            let similar: HashSet<_> = set_to_check.intersection(candidate).cloned().collect();

                            if !similar.is_empty() {
                                let similarity = similar.len() as f32 / usize::min(set_to_check.len(), candidate.len()) as f32;

                                if similarity >= kmer_percent {
                                    if iteration == 0 || (iteration == 1 && cluster_children[master_header].len() != 1) {
                                        child_sets.get_mut(master_header).unwrap().insert(candidate_header.clone());
                                        let extend_slice = cluster_children.remove(candidate_header).unwrap();
                                        cluster_children.get_mut(master_header).unwrap().extend(extend_slice);

                                        kmers.remove(candidate_header);
                                        processed_headers.insert(candidate_header.clone());

                                        matched = true;
                                        break;
                                    }
                                }
                            }
                        }

                        if matched {
                            merge_occurred = true;
                        }
                    }
                }
            }
        }
    }

    cluster_children.values().cloned().collect()

    // the following lines get the values without cloning the values
    // if they work, it should be much faster
    // let keys = cluster_children.keys().clone().collect();
    // keys.iter().map(|key| cluster_children.remove(key).unwrap()).collect()
}
#[pyfunction]
pub fn seperate_into_clusters(
    cluster_children: Vec<Vec<String>>,
    data: HashMap<String, String>,
    subcluster_at: usize,
    cluster_every: usize,
) -> Vec<Vec<String>> {
    let mut clusters = Vec::new();
    for this_cluster in cluster_children {
        // guard statement: if not greater than subcluster_at, append and continue
        // else make subclusters
        if this_cluster.len() <= subcluster_at {
            clusters.push(this_cluster);
            continue;
        }
        let clusters_to_create = (this_cluster.len() as f32/ cluster_every as f32).ceil() as usize;
        let mut temp_in = NamedTempFile::new().expect("Cannot create temp file for Sigclust input");
        for header in this_cluster.iter() {
            temp_in.write(format!(">{}", header).as_bytes()).unwrap();
            temp_in.write(data.get(header).unwrap().as_bytes()).unwrap();
        }
        temp_in.flush().unwrap();

        let temp_path = temp_in.path().to_str().unwrap();
        let mut sigclust = process::Command::new("siglclust/Sigclust")
            .args(["-k", "8", "-c", &clusters_to_create.to_string(), temp_path])
            .output()
            .expect("Cannot get output from sigclust run")
            .stdout;
        let output = String::from_utf8(
            sigclust).unwrap();
        let split_position = output.find("Writing output").unwrap() + 14; //size of "Writing output"
        let lines = output[split_position..].split('\n').collect::<Vec<&str>>();
        let mut subcluster = HashMap::new();
        for line in lines {
            let nums = line.split(',').collect::<Vec<&str>>();
            let seq_index = nums[0].parse::<usize>().expect("Bad sigclust data");
            let clust_index = nums[1].parse::<usize>().expect("Bad sigclust data");
            subcluster.entry(clust_index).or_insert_with(Vec::new)
                .push(this_cluster[seq_index].clone())
        }
        for value in subcluster.into_values(){
            clusters.push(value)
        };
    };
    clusters
}
