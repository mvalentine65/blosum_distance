use bio::pattern_matching::pssm::{Motif, ProtMotif};
use pyo3::prelude::*;
#[pyclass]
struct ScoredPosition {
    #[pyo3(get, set)]
    pub location: usize,
    #[pyo3(get, set)]
    pub sum: f32,
    #[pyo3(get, set)]
    pub scores: Vec<f32>
}

#[pymethods]
impl ScoredPosition {
    #[new]
    pub fn __init__(location: usize, sum: f32, scores: Vec<f32>) -> ScoredPosition {
        ScoredPosition {
            location,
            sum,
            scores,
        }
    }
}


#[pyclass]
struct ProteinMotif {
    pub matrix: ProtMotif,
}

#[pymethods]
impl ProteinMotif {
    #[new]
    pub fn __init__(sequences: Vec<String>) -> ProteinMotif {
        let seqs: [Vec<u8>] = &sequences.iter().map(|&s| s.as_bytes().to_vec()).collect();
        ProteinMotif {
            matrix: ProtMotif::from_seqs(&seqs, None).unwrap()
        }
    }
    
    pub fn get_scores(&self) -> Vec<f32> {
        self.matrix.get_scores().to_owned().into_raw_vec()
    }
    
     pub fn min_score(&self) -> f32 {
        self.matrix.get_min_score()
    }
    
    pub fn max_score(&self) -> f32 {
        self.matrix.get_max_score()
    }
    
    pub fn degenerate_consensus(&self) -> Vec<u8> {
        self.matrix.degenerate_consensus()
    }
    
    pub fn lookup(&self, mono: u8) -> Option<usize> {
        match self.matrix.lookup(mono) {
            Ok(m) => m,
            Err(_) => None,
        }
    }
    
    pub fn len(&self) -> usize {
        self.matrix.len()
    }
    
    pub fn rev_lk(&self, index: usize) -> u8 {
        self.matrix.rev_lk(index)
    }
    
    pub fn is_empty(&self) -> bool {
        self.matrix.is_empty()
    }
    
    pub fn raw_scores(&self, query_sequence: String) -> Option<(usize, f32, Vec<f32>)> {
        match self.matrix.raw_score(query_sequence.iter()) {
            Ok(tuple) => tuple,
            Err(_) => None,
        }
    }
    
    pub fn score(&self, query_sequence: String) -> Option<ScoredPosition> {
        match self.matrix.score(query_sequence.iter()) {
            Ok(scored) => ScoredPosition{
                                        location: scored.loc,
                                        sum: scored.sum,
                                        scores: scored.scores,
                                        },
            Err(_) => None,
        }
    }
    
    pub fn info_content(&self) -> f32 {
        self.matrix.info_content()
    }
    
    pub fn __eq__(&self, other: ProteinMotif) -> bool {
        self.matrix.eq(&other.matrix)
    }
    
    pub fn __ne__(&self, other: ProteinMotif) -> bool {
        self.matrix.ne(&other.matrix)
    }
}