// mod utils;
// use crate::bio_revcomp;

use bio;
use bio::alphabets::dna::revcomp;
use gxhash;
use gxhash::{GxHashMap, GxHashSet};
// use itertools::enumerate;
use pyo3::prelude::*;
use std::collections::{HashSet};
// use std::io::{BufRead, BufReader};
// use std::path::Path;
//use bio::io::fasta::{Record, Records};
// use fastx::R
// use utils::FastaParser;
// fn has_lowercase(sequence: &String) -> bool {
//     sequence.as_bytes().any(|c| (c >= b'a') && (c <= b'z'))
// }

#[pyclass]
pub struct PyGxSet {
    gx_set: HashSet<String, gxhash::GxBuildHasher>,
}

#[pymethods]
impl PyGxSet {
    #[new]
    fn new() -> Self {
        PyGxSet {
            gx_set: GxHashSet::default(),
        }
    }

    fn __contains__(&self, item: String) -> bool {
        self.gx_set.contains(&item)
    }

    fn contains_with_revcomp(&self, item: String) -> bool {
        self.gx_set.contains(&item)
            || self.gx_set.contains(
                &String::from_utf8(revcomp(item.into_bytes())).unwrap(),
            )
    }

    fn add(&mut self, item: String) {
        self.gx_set.insert(item);
    }

    fn add_with_revcomp(&mut self, item: String) {
        self.gx_set.insert(item.clone());
        self.gx_set
            .insert(String::from_utf8(revcomp(item.into_bytes())).unwrap());
    }
}

#[pyclass]
pub struct PyGxDict {
    hashmap: GxHashMap<String, String>
}

#[pymethods]
impl PyGxDict {
    #[new]
    fn new() -> Self {
        PyGxDict {
            hashmap: GxHashMap::default()
        }
    }

    fn __contains__(&self, key: String) -> bool {
        self.hashmap.contains_key(&key)
    }

    fn add(&mut self, key: String, value:String) {
        self.hashmap.insert(key, value);
    }

    fn get(&self, key:String) -> Option<String> {
        match self.get(key) {
            Some(value) => Some(value),
            None => None,
        }
    }
}

#[pyclass]
pub struct PyGxCounter {
    counter: GxHashMap<String, i32>
}

#[pymethods]
impl PyGxCounter {
   #[new]
   fn new() -> Self {
       PyGxCounter {
           counter: GxHashMap::default()
       }
   }
    fn __contains__(&self, key: String) -> bool {
        self.counter.contains_key(&key)
    }
    
    fn add(&mut self, key: String) {
        *self.counter.entry(key).or_insert(0) += 1
    }
    
    fn get(&self, key: String) -> i32 {
        match self.counter.get(&key) {
            Some(value) => *value,
            None => 0
        }
        
    }
}
//
// #[pyclass]
// struct SeqDeduper {
//     original_positions: GxHashMap<K, V>,
//     original_inputs: Vec<String>,
//     lines: Vec<String>,
//     file_index: usize,
//     min_seq_len: usize,
//     verbose: bool,
//     is_assembly: bool,
//     is_genome: bool,
//     overlap_length: usize,
//     rename: bool,
//     transcripts_mapped_to: GxHashMap<String, String>,
//     hash_set: GxHashSet<String>,
//     dupes: usize,
//     rev_comp_save: GxHashMap<String, usize>,
//     duplicates: GxHashMap<String, usize>,
// }
//
// #[pymethods]
// impl SeqDeduper {
//     #[new]
//     fn new(minimum_seq_len: usize, verbosity: bool, overlap_length: usize, rename: bool) -> Self {
//         SeqDeduper {
//             original_positions: GxHashMap::default(),
//             original_inputs: Vec::new(),
//             lines: Vec::new(),
//             file_index: 0,
//             min_seq_len: minimum_seq_len,
//             verbose: verbosity,
//             is_assembly: false,
//             is_genome: false,
//             overlap_length,
//             rename,
//             transcripts_mapped_to: GxHashMap::default(),
//             hash_set: GxHashSet::default(),
//             dupes: 0,
//             rev_comp_save: GxHashMap::default(),
//             duplicates: GxHashMap::default(),
//         }
//     }
//
//     fn __call__(
//         &mut self,
//         fa_file_path: String,
//         duplicates: HashMap<String, usize>,
//         rev_comp_save: HashMap<String, usize>,
//         this_index: usize,
//     ) {
//         const chomp_len: usize = 750;
//         const chomp_cutoff: usize = 200000;
//         const assembly_len: usize = 750;
//
//         self.original_inputs.push(fa_file_path.clone());
//         // TODO: progress bar here
//         let reader = fastx::FastX::reader_from_path(Path::new(&fa_file_path)).unwrap();
//         let mut requires = false;
//         for (line_index, record) in enumerate(reader.records()) {
//             let record = record.unwrap();
//             let header = record.id();
//             let sequence = record.sequence();
//             if line_index == 0 {
//                 requires = has_lowercase(&sequence);
//             }
//             if len
//         }
//     }
// }
//
// #[pyclass]
// struct GxDeduper {
//     #[pyo3(get, set)]
//     chomp_len: usize,
//     #[pyo3(get, set)]
//     chomp_cutoff: usize,
//     #[pyo3(get, set)]
//     assembly_len: usize,
//     #[pyo3(get, set)]
//     rename: bool,
//     #[pyo3(get, set)]
//     min_seq_len: usize,
//     #[pyo3(get, set)]
//     this_index: usize,
//     // header_template: String,
//     append_index_template: String,
//     sequence_template: String,
//     duplicates: GxHashMap<String, usize>,
//     transcripts_mapped_to: GxHashMap<String, String>,
//     dupe_count: usize,
//     rev_comp_map: GxHashMap<String, i64>,
//     #[pyo3(get, set)]
//     is_assembly: bool,
//     is_genome: bool,
// }
//
// #[pymethods]
// impl GxDeduper {
//     #[new]
//     fn new(rename: bool, min_seq_len: usize, is_assembly: bool, is_genome: bool) -> Self {
//         GxDeduper {
//             chomp_len: 750,
//             chomp_cutoff: 200000,
//             assembly_len: 750,
//             rename,
//             min_seq_len,
//             this_index: 0,
//             // header_template: "NODE_{}".to_string(),
//             append_index_template: "{}_{}".to_string(),
//             sequence_template: ">{}\n{}\n".to_string(),
//             duplicates: GxHashMap::<String, usize>::default(),
//             transcripts_mapped_to: GxHashMap::<String, String>::default(),
//             dupe_count: 0,
//             rev_comp_map: GxHashMap::<i32, i32>::default(),
//             is_assembly,
//             is_genome,
//         }
//     }
//
//     fn n_trim(&self, parent_seq: String) -> Vec<String> {
//         parent_seq
//             .split("N")
//             .map(|str| str.to_string())
//             .filter(|chunk| chunk.len() > self.min_seq_len)
//             .collect()
//     }
//     fn apply(&mut self, _header: String, parent_seq: String, line_index: usize) {
//         if parent_seq.len() < self.min_seq_len {
//             return;
//         }
//         let parent_seq = parent_seq.to_uppercase();
//         let n_sequences = self.n_trim(parent_seq);
//         let header: String;
//         match self.rename {
//             true => header = format!("NODE_{}", self.this_index),
//             false => header = header.split(" ").take(1).collect(),
//         }
//         for seq in n_sequences.into_iter() {
//             // if self.transcripts_mapped_to.contains_key(&seq) {
//             //     self.duplicates[self.transcripts_mapped_to]
//             // }
//             let mapped_header = self.transcripts_mapped_to.get(&seq);
//             if mapped_header.is_some() {
//                 *self
//                     .duplicates
//                     .entry(mapped_header.unwrap().to_string())
//                     .or_insert(0) += 1;
//                 self.dupe_count += 1;
//                 continue;
//             }
//             self.transcripts_mapped_to
//                 .insert(seq.to_string(), header.clone());
//
//             let rev_comp = bio_revcomp(seq.clone());
//             let rev_header = self.transcripts_mapped_to.get(&rev_comp);
//             if rev_header.is_some() {
//                 *self
//                     .duplicates
//                     .entry(rev_header.unwrap().to_string())
//                     .or_insert(0) += 1;
//                 self.dupe_count += 1;
//                 continue;
//             }
//             self.transcripts_mapped_to
//                 .insert(rev_comp.clone(), header.clone());
//             if (!self.is_assembly && !self.is_genome) && seq.len() >= self.assembly_len {
//                 self.is_assembly = true;
//             }
//             if seq.len() > self.chomp_cutoff {
//                 if !self.is_genome {
//                     self.is_genome = true;
//                     self.is_assembly = false;
//                 }
//             }
//         }
//     }
// }
