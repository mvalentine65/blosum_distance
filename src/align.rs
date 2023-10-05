use std::collections::{HashMap, HashSet};
use std::{fs, process};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use fastx::FastX;
use fastx::FastX::FastXFormat;
use fastx::FastX::FastXFormat::FASTA;
use tempfile::NamedTempFile;
use pyo3::prelude::*;


fn count_sequences(path: &str) -> usize {
    let mut count = 0_usize;
    let mut file = File::open(path).unwrap();
    let mut reader = BufReader::new(file);

    for line in reader.lines() {
        if let Ok(line) = line {
            if line.starts_with('>') {
                count += 1
            }
        }
    }
    count
}
fn writeFastaFileUncompressed<T: AsRef<str>>(file: &mut File, records: &Vec<(T, T)>) {
    for (header, sequence) in records {
        write!(file, ">{}\n{}\n", header.as_ref(), sequence.as_ref()).unwrap()
    }
}
fn writeFastaPathUncompressed<T: AsRef<str>>(path: &String, records: &Vec<(T, T)>) {
    // let mut target_path = path.to_string();
    // if !target_path.ends_with(".fa"){
    //     target_path = target_path + ".fa";
    // }
    let mut file = File::create(path).unwrap();
    for (header, sequence) in records {
        // line = format!(">{}\n;{}\n", header.as_ref(), sequence.as_ref());
        // buf = line.as_bytes();
        // file.write_all(buf).unwrap();
        file.write(format!(">{}\n", header.as_ref()).as_bytes()).unwrap();
        file.write(format!("{}\n", sequence.as_ref()).as_bytes()).unwrap();
    }
    file.flush().unwrap();

}

fn default_clust_prep(dict: &mut HashMap<String, HashSet<String>>, key: &str) {
   if !dict.contains_key(key) {
       let new = HashSet::new();
       dict.insert(key.to_string(), new);
   }
}

fn parse_fasta(path: &String) -> Vec<(String, String)> {
    let mut records = Vec::<(String, String)>::new();
    let mut reader = FastX::reader_from_path(Path::new(path)).unwrap();
    let mut fastx_record = FastX::from_reader(&mut reader).unwrap();
    while let Ok(_some @ 1..=usize::MAX) = fastx_record.read(&mut reader) {
        records.push((fastx_record.id().to_string(), String::from_utf8(fastx_record.seq()).unwrap()));
    }
    records


}

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
                default_clust_prep(&mut kmers, master_header);
                // let master = match kmers.get(master_header) {
                //     Some(set) => set,
                //     None => continue,
                // };
                let master = kmers.get(master_header).unwrap();
                if !master.is_empty() {
                    for j in 0..i {
                        let candidate_header = &gene_headers[j];

                        if processed_headers.contains(candidate_header) {
                            continue;
                        }
                        default_clust_prep(&mut kmers, candidate_header);
                        let candidate = kmers.get(candidate_header).unwrap();
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

#[pyfunction]
pub fn make_aligned_ingredients(clusters: Vec<Vec<String>>, data: HashMap<String, String>, gene: String, aligned_files_tmp: String
, raw_files_tmp: String, this_intermediates: String) -> Vec<(String, usize, usize)> {
    let mut aligned_ingredients:Vec<(String, usize, usize)> = Vec::new();
    for (cluster_i, cluster) in clusters.iter().enumerate() {
        let cluster_seqs: Vec<(&str, &str)> = cluster.iter()
            .map(|header| (header.as_str(), data.get(header).unwrap().as_str()))
            .collect();

        let cluster_length = cluster.len();

        let this_clus_align = format!("aligned_cluster_{}", cluster_i);
        let mut aligned_cluster = PathBuf::from(&aligned_files_tmp);
        aligned_cluster = aligned_cluster.join(&this_clus_align);
        let mut raw_cluster = PathBuf::from(&raw_files_tmp);
        raw_cluster = raw_cluster.join(format!("{}_cluster{}", gene, cluster_i));
        if cluster_length == 1 {
            writeFastaPathUncompressed(&aligned_cluster.to_string_lossy().to_string(), &cluster_seqs);
            aligned_ingredients.push((aligned_cluster.to_string_lossy().to_string(), 1, cluster_i));
            continue;
        }
        // if cluster_seqs.len() > 1...
        writeFastaPathUncompressed(&raw_cluster.to_string_lossy().to_string(), &cluster_seqs);
        // make_file(&aligned_cluster);
        // let mut out_file = File::create(&aligned_cluster).unwrap();
        // let mut intermediate_path = PathBuf::from(&this_intermediates);
        // intermediate_path = intermediate_path.join(this_clus_align);
        // writeFastaUncompressed(&intermediate_path.to_string_lossy().to_string(), cluster_seqs);
        // command = "clustalo -i {in_file} -o {out_file} --threads=1 --full"

        let _status = process::Command::new("clustalo")
            .args(&["-i", &raw_cluster.to_string_lossy().to_string(),
                "-o", &aligned_cluster.to_string_lossy().to_string(),
                "--threads=1", "--full", "--force"])
            .output()
            .expect("Failed to execute clustalo");
        // let aligned_sequences = parse_fasta(&aligned_cluster.to_string_lossy().to_string());
        let count = count_sequences(&aligned_cluster.to_string_lossy().to_string());
        // if debug:
        //     printv(command, args.verbose, 3)
        // writeFasta(
        //     path.join(
        //         this_intermediates,
        //         this_clus_align,
        //     ),
        //     aligned_sequences,
        // )

        aligned_ingredients.push((aligned_cluster.to_string_lossy().to_string(), count, cluster_i));
    }
    aligned_ingredients
}

#[pyfunction]
pub fn run_intermediate(cluster_file: String, seq_count: usize, cluster_i: usize,  tmp_align: String,
                     parent_tempdir: String, this_intermediates: String, debug: bool) {
    let outfile = Path::new(&parent_tempdir);
    let outfile = outfile.join(format!("part_{}.fa", cluster_i)).to_string_lossy().to_string();
    let args;
    if seq_count >= 1 {
        args = ["--p1", &tmp_align, "--p2", &cluster_file, "-o", &outfile, "--threads=1", "--full", "--is-profile", "--force"];
    } else {
        args = ["--p1", &tmp_align, "--p2", &cluster_file, "-o", &outfile, "--threads=1", "--full", "", "--force"]
    }


    let output = std::process::Command::new("clustalo")
        .args(args)
        .output()
        .expect("Failed to execute clustalo")
        .stdout;

    // if debug:
    //     printv(
    //         f"clustalo --p1 {tmp_aln.name} --p2 {file} -o {out_file} --threads=1 --full {is_profile} --force",
    // args.verbose,
    // 3,
    // )
    if debug {
        let debug_path = Path::new(&this_intermediates);
        let debug_path = debug_path.join(format!("reference_subalignment_{}.fa", cluster_i));
        writeFastaPathUncompressed(&debug_path.to_string_lossy().to_string(), &parse_fasta(&outfile));
    }
}


struct RecordHolder {
    header: Vec<u8>,
    sequence: Vec<u8>,
    // in_header: bool
}

impl RecordHolder {
    fn to_tuple(&self) -> (String, String) {
        (String::from_utf8_lossy(&self.header).to_string(), String::from_utf8_lossy(&self.sequence).to_string())
    }
}


// pub fn parse_stdout(stdout: Vec<u8>) -> Vec<String> {
//     let mut output = Vec::new();
//     // let mut current:RecordHolder = RecordHolder{
//     //     header: Vec::new(),
//     //     sequence: Vec::new(),
//     // };
//     let mut current = &Vec::new();
//     let mut in_header = true;
//     for byte in stdout {
//         if byte == b'>' {
//             output.push(Vec::new());
//             current = output.last().unwrap();
//             in_header = true;
//             }
//             // output.push(&current);
//         } else if byte == b'\n' && in_header == true {
//                 output.push(Vec::new());
//                 in_header = false;
//                 current = output.last().unwrap()
//         } else {
//             current.push(byte);
//         }
//
//     output.into_iter()
//         .map(|record| {
//             String::from_utf8(record).unwrap()
//         })
//         .collect()
// }
#[pyfunction]
pub fn process_clusters(clusters: Vec<Vec<String>>, data: HashMap<String, String>, this_intermediates: String
, parent_tempdir: String, tmp_align: String) -> Vec<(String, String)> {
    // let p2_path = Path::new(p2);

    // if !p2_path.exists() {
    //     println!("Cluster file '{}' does not exist. Skipping...", p2);
    //     return;
    // // }
    // if !has_multiple_sequences(p2_path) {
    //     return;
    // }
    let mut final_sequences = Vec::new();
    for (cluster_i, cluster) in clusters.iter().enumerate() {
        let cluster_seqs: Vec<(&String, &String)> = cluster.iter()
            .map(|header| (header, data.get(header).unwrap()))
            .collect();
        if cluster_seqs.len() == 1 {
            let (header, sequence) = cluster_seqs[0];
            final_sequences.push((header.to_owned(), sequence.to_owned()));
            continue;
        }
        let mut unaligned_in = PathBuf::from(&this_intermediates);
        unaligned_in = unaligned_in.join(format!("unaligned_{}", cluster_i));
        let mut in_file = NamedTempFile::new_in(&this_intermediates).unwrap();
        for (header, seq) in cluster_seqs {
            write!(&mut in_file, ">{}\n{}\n", header, seq).unwrap()
        }
        let in_path = in_file.path().to_str().expect("cannot get path to temp infile");
        // let p2_number = p2.trim_start_matches("unaligned_cluster_");
        let intermediate_file = NamedTempFile::new_in(&this_intermediates).expect("cannot create temp file");
        // (format!("cluster_{}.fasta", p2_number));
        let temp_path = intermediate_file.path().to_str().expect("cannot get path to temp file");
        let _status = std::process::Command::new("clustalo")
            .args(&["-i", in_path, "-o", temp_path, "--threads=1", "--force", "--full"])
            .status()
            .expect("Failed to execute clustalo");

        // if !status.success() {
        //     eprintln!("Failed to process unaligned cluster file '{}'", p2);
        // }
        // println!("first is fine {}",cluster_i);
        let mut out_file = PathBuf::from(&parent_tempdir);
        let out_tail = format!("part_{}.fa", cluster_i);
        out_file = out_file.join(out_tail);
        let out_path = out_file.as_path().to_str().unwrap();
        let args;
        let seq_count = count_sequences(temp_path);
        if seq_count >= 1 {
            args = ["--p1", &tmp_align, "--p2", temp_path, "-o", out_path, "--threads=1", "--full", "--is-profile", "--force"];
        } else {
            args = ["--p1", &tmp_align, "--p2", temp_path, "-o", out_path, "--threads=1", "--full", "", "--force"]
        }
        let _status = std::process::Command::new("clustalo")
            .args(&["--p1", &tmp_align, "--p2", temp_path, "-o", &out_path, "--threads=1", "--full", "--is-profile", "--force"])
            .status()
            .expect("Failed to execute clustalo");

        // if !status.success() {
        //     eprintln!("Failed to process aligned cluster file '{}'", p2);
        // }
    }
    final_sequences
}