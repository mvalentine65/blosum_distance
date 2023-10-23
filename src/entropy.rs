use pyo3::pyfunction;

// Accepts fastas as the standard tuple[str, str] and an entropy threshold.
// Returns two tuple[str, str],
// one for sequences above the entropy threshold, and one for sequences below the threshold
#[pyfunction]
pub fn entropy_filter(records: Vec<String>, entropy_threshold: f64) -> Vec<(String, String)> {
    let mut greater = Vec::with_capacity(records.len());
    for string in records.into_iter() {
        let fields = string.split("\n").collect::<Vec<&str>>();
        let header = fields[0].to_string();
        let seq = fields[1].to_string();
        // Calculate the frequency distribution of nucleotides in the sequence
        let a_count = seq.bytes().filter(|&c| c == b'A').count() as u32;
        let c_count = seq.bytes().filter(|&c| c == b'C').count() as u32;
        let g_count = seq.bytes().filter(|&c| c == b'G').count() as u32;
        let t_count = seq.bytes().filter(|&c| c == b'T').count() as u32;
        let counts = [
            a_count as f64,
            c_count as f64,
            g_count as f64,
            t_count as f64,
        ];

        // Calculate the entropy of the sequence
        let freqs = counts.map(|count| count as f64 / seq.len() as f64);
        let ent = _entropy(&freqs);

        // Filter sequences with low entropy and write the remaining ones to the output FASTA file
        if ent >= entropy_threshold {
            greater.push((header, seq))
        }
    }
    greater
}

// Helper function to calculate entropy of a frequency distribution
fn _entropy(freqs: &[f64]) -> f64 {
    -freqs.iter().map(|&p| p * p.log2()).sum::<f64>()
}
#[pyfunction]
pub fn entropy(freqs: Vec<f64>) -> f64 {
    -freqs.iter().map(|&p| p * p.log2()).sum::<f64>()
}
