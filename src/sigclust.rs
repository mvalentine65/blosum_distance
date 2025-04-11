use pyo3::pyfunction;
use rand::distributions::{Distribution,Uniform};
use rand::{Rng, SeedableRng};
use std::collections::HashSet;
use rand_chacha::ChaCha8Rng;
use std::collections::hash_map::DefaultHasher;
use std::hash::Hasher;
//f"sigclust/SigClust -k 8 -c {clusters_to_create} {this_tmp.name} > {this_out.name}",

const SIGNATURE_WIDTH: usize = 256;
const SIGNATURE_SIZE: usize = SIGNATURE_WIDTH / 64;
const DENSITY: f32 = 1_f32 / 21_f32;
const KMEANS_ITERATIONS: usize = 4;

#[pyfunction]
pub fn sigclust(records: Vec<(String, String)>, k: usize, c: usize) -> Vec<Vec<String>> {
    let sigs = convert_fasta_to_signatures(&records, k);
    let clusters = cluster_signatures(&sigs, c);
    prepare_header_output(&records, &clusters, c)
}

#[pyfunction]
pub fn sigclust_with_sequence(records: Vec<(String, String)>, k: usize, c: usize) -> Vec<Vec<(String, String)>> {
    let sigs = convert_fasta_to_signatures(&records, k);
    let clusters = cluster_signatures(&sigs, c);
    prepare_tuple_output(&records, &clusters, c)
}

fn generate_signature(output: &mut [u64], sequence: &str, kmer_length: usize) {
    let mut seq_vec = sequence.as_bytes().to_vec();
    while seq_vec.len() < kmer_length {
        seq_vec.push(b'X');
    }
    // upper bound in *inclusive* in c++, and *exclusive* in rust
    // cpp code subtracts one from upper bound, dont do that here

    let mut unflattened_signature = vec![0_i32; SIGNATURE_WIDTH];
    let mut uniform_distribution = Uniform::new(
        -64_i32 * SIGNATURE_SIZE as i32,
        64_i32 * SIGNATURE_SIZE as i32,
    );
    let set_bits = (DENSITY * SIGNATURE_WIDTH as f32) as i32;
    for i in 0..(seq_vec.len() - kmer_length + 1) {

        let mut hash = DefaultHasher::new();
        hash.write(&seq_vec[i..i + kmer_length]);
        let seed = hash.finish();
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        for _ in 0..set_bits {
            let bit_pos = uniform_distribution.sample(&mut rng);
            if bit_pos >= 0 {
                unflattened_signature[bit_pos as usize] += 1;
            } else {
                unflattened_signature[(bit_pos + SIGNATURE_WIDTH as i32) as usize] -= 1
            }
        }
    }
    for i in 0..SIGNATURE_SIZE {
        output[i] = 0;
    }
    for i in 0..SIGNATURE_WIDTH {
        if unflattened_signature[i] > 0 {
            output[i / 64] |= 1_u64 << (i % 64);
        }
    }
}

fn convert_fasta_to_signatures(records: &Vec<(String, String)>, kmer_length: usize) -> Vec<u64> {
    let mut output = vec![0_u64; records.len() * SIGNATURE_SIZE];
    for i in 0..records.len() {
        generate_signature(
            &mut output[SIGNATURE_SIZE * i..SIGNATURE_SIZE * i + SIGNATURE_SIZE],
            &records[i].1,
            kmer_length,
        );
    }
    output
}

fn create_cluster_lists(clusters: &Vec<usize>, cluster_count: usize) -> Vec<Vec<usize>> {
    let mut cluster_lists: Vec<Vec<usize>> = vec![Vec::new(); cluster_count];
    for i in 0..clusters.len() {
        cluster_lists[clusters[i]].push(i);
    }
    cluster_lists
}

fn create_cluster_sigs(
    cluster_lists: &Vec<Vec<usize>>,
    sigs: &Vec<u64>,
    cluster_count: usize,
) -> Vec<u64> {
    let mut cluster_sigs = vec![0_u64; SIGNATURE_SIZE * cluster_count];
    let mut unflattened_signature;

    for cluster_i in 0..cluster_lists.len() {
        unflattened_signature = vec![0_i32; SIGNATURE_WIDTH];

        for signature in &cluster_lists[cluster_i] {
            let signature_data =
                &sigs[SIGNATURE_SIZE * signature..SIGNATURE_SIZE * (signature+1)];
            for sig_i in 0..SIGNATURE_WIDTH {
                let signature_mask = 1_u64 << (sig_i % 64);
                if (signature_mask & signature_data[sig_i / 64]) != 0 {
                    unflattened_signature[sig_i] += 1;
                } else {
                    unflattened_signature[sig_i] -= 1
                }
            }
        }
        let mut flattened_signature = &mut cluster_sigs[cluster_i * SIGNATURE_SIZE..];
        for i in 0..SIGNATURE_WIDTH {
            if unflattened_signature[i] > 0 {
                flattened_signature[i / 64] |= 1 << (i % 64);
            }
        }
    }
    cluster_sigs
}

fn recluster_signatures(
    clusters: &mut Vec<usize>,
    mean_sigs: &Vec<u64>,
    sigs: &Vec<u64>,
    cluster_count: usize,
) {
    for sig_i in 0..clusters.len() {
        let source_signature = &sigs[sig_i * SIGNATURE_SIZE..];
        let mut min_hd_cluster = 0_usize;
        let mut min_hd = usize::MAX;
        for cluster_i in 0..cluster_count {
            let cluster_signature = &mean_sigs[cluster_i * SIGNATURE_SIZE..];
            let mut hd: usize = 0;
            for i in 0..SIGNATURE_SIZE {
                hd += (source_signature[i] ^ cluster_signature[i]).count_ones() as usize;
            }
            if hd < min_hd {
                min_hd = hd;
                min_hd_cluster = cluster_i
            }
        }
        clusters[sig_i] = min_hd_cluster;
    }
}

fn create_random_sigs<R: Rng>(rng: &mut R, sigs: &Vec<u64>, cluster_count: usize) -> Vec<u64> {
    let mut cluster_sigs = vec![0_u64; SIGNATURE_SIZE * cluster_count];
    let signature_count = sigs.len() / SIGNATURE_SIZE;
    let uniform_int_distribution = Uniform::new(0_usize, signature_count);
    let mut finished = false;

    let mut unique_sigs_set = HashSet::new(); // line 179 in cpp
    let mut unique_sigs_vec = Vec::new();
    for _sig in 0..signature_count {
        let sig = rng.sample(uniform_int_distribution);
        let sig_data = sigs[sig * SIGNATURE_SIZE..(sig+1) * SIGNATURE_SIZE].to_vec();
        if !unique_sigs_set.contains(&sig_data) {
            unique_sigs_set.insert(sig_data.clone());
            unique_sigs_vec.push(sig_data);
        }

        if unique_sigs_set.len() >= cluster_count {
            finished = true;
            break;
        }
    }
    if !finished {
        panic!("Fewer unique signatures than clusters, please lower cluster count.")
    }

    let mut i = 0_usize;
    for sig in unique_sigs_vec.iter() {
        cluster_sigs[i * SIGNATURE_SIZE..(i + 1) * SIGNATURE_SIZE]
            .copy_from_slice(sig.as_slice());
        i += 1;
    }

    if i != unique_sigs_vec.len() {
        eprintln!(
            "Mismatch {} sigs added and {} unique sigs",
            i,
            unique_sigs_vec.len()
        );
        panic!();
    }
    cluster_sigs
}

fn cluster_signatures(sigs: &Vec<u64>, cluster_count: usize) -> Vec<usize> {
    let mut rng = ChaCha8Rng::seed_from_u64(19780503);
    let mut clusters = vec![0_usize; sigs.len() / SIGNATURE_SIZE];
    let mut clusters_list = Vec::new();
    let mut mean_sigs = create_random_sigs(&mut rng, sigs, cluster_count);
    for _kmeans_iteration in 0..KMEANS_ITERATIONS {
        recluster_signatures(&mut clusters, &mean_sigs, sigs, cluster_count);
        clusters_list = create_cluster_lists(&clusters, cluster_count);
        mean_sigs = create_cluster_sigs(&clusters_list, sigs, cluster_count);
    }

    clusters
}

fn prepare_header_output(
    records: &Vec<(String, String)>,
    clusters: &Vec<usize>,
    cluster_count: usize,
) -> Vec<Vec<String>> {
    let mut output = vec![Vec::new(); cluster_count];
    for i in 0..clusters.len() {
        output[clusters[i]].push(records[i].0.clone());
    }
    output
}

fn prepare_tuple_output(
    records: &Vec<(String, String)>,
    clusters: &Vec<usize>,
    cluster_count: usize,
) -> Vec<Vec<(String, String)>> {
    let mut output = vec![Vec::new(); cluster_count];
    for i in 0..clusters.len() {
        output[clusters[i]].push(records[i].clone());
    }
    output
}
