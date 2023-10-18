use rand::distributions::Uniform;
use rand::rngs::SmallRng;
use rand::{thread_rng, Rng, SeedableRng};

//f"sigclust/SigClust -k 8 -c {clusters_to_create} {this_tmp.name} > {this_out.name}",



const SIGNATURE_WIDTH: usize = 256;
const SIGNATURE_SIZE: usize = SIGNATURE_WIDTH / 64;
const DENSITY: f32 = 1_f32 / 21_f32;
const KMEANS_ITERATIONS: usize = 4;
pub fn sigclust(records: Vec<(String, String)>, k: usize, c: usize) -> Vec<Vec<(String, String)>> {
    let output = Vec::new();

    output
}

fn generate_signature(output: &mut [u64], sequence: &str, kmer_length: usize) {
    let mut seq_vec = sequence.as_bytes().to_vec();
    while seq_vec.len() < kmer_length {
        seq_vec.push(b'X');
    }
    let mut rng;
    // upper bound in *inclusive* in c++, and *exclusive* in rust
    // cpp code subtracts one from upper bound, dont do that here
    let uniform_distribution = Uniform::new(
        -64_i32 * SIGNATURE_WIDTH as i32,
        64_i32 * SIGNATURE_WIDTH as i32,
    );
    let mut unflattened_signature = Vec::<i32>::new();
    let mut seed = [0; 32];
    let set_bits = (DENSITY * SIGNATURE_SIZE as f32 * 64.0) as i32;
    let mut kmer;
    for i in 0..(seq_vec.len() - kmer_length + 1) {
        // let seed = &seq_vec[i..i+kmer_length]
        seed[..kmer_length].copy_from_slice(&seq_vec[i..i + kmer_length]);
        rng = SmallRng::from_seed(seed);
        kmer = &seq_vec[i..i + kmer_length];
        for j in 0..set_bits {
            let bit_pos = rng.sample(uniform_distribution);
            if bit_pos >= 0 {
                unflattened_signature[bit_pos as usize] += 1;
            } else {
                unflattened_signature[(bit_pos + 64 * SIGNATURE_SIZE as i32) as usize] -= 1
            }
        }
    }
    for i in 0..SIGNATURE_SIZE {
        output[i] = 0;
    }
    for i in 0..SIGNATURE_SIZE * 64 {
        if unflattened_signature[i] > 0 {
            output[i / 64] |= 1 << (i % 64);
        }
    }
}

fn convert_fasta_to_signatures(records: &Vec<(String, String)>, kmer_length: usize) -> Vec<u64> {
    let mut output = vec![0_u64;records.len() * 64];
    for i in 0..records.len() {
        generate_signature(&mut output[SIGNATURE_SIZE * i..SIGNATURE_SIZE * i + SIGNATURE_SIZE],
                            &records[i].1, kmer_length);
    }
    output
}

fn create_cluster_lists(clusters: &Vec<usize>, cluster_count: usize) -> Vec<Vec<usize>> {
    let mut cluster_lists:Vec<Vec<usize>> = vec![Vec::new();cluster_count];
    for i in 0..clusters.len() {
        cluster_lists[clusters[i]].push(i);
    }
    cluster_lists
}