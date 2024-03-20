
use bio;
use gxhash;
use pyo3::prelude::*;
use std::collections::HashSet;

#[pyclass]
pub struct GxHasher {
    gx_set: HashSet<String, gxhash::GxBuildHasher>
}


#[pymethods]
impl GxHasher {

    #[new]
    fn new() -> Self {
        GxHasher {
            gx_set: gxhash::GxHashSet::default()
        }
    }


    fn __contains__(&self, item: String) -> bool {
        self.gx_set.contains(&item)
    }

    fn contains_with_revcomp(&self, item: String) -> bool {
        self.gx_set.contains(&item) ||
            self.gx_set.contains(&String::from_utf8(bio::alphabets::dna::revcomp(
                item.into_bytes()
            )).unwrap())
    }

    fn add(&mut self, item: String) {
        self.gx_set.insert(item);
    }

    fn add_with_revcomp(&mut self, item: String) {
        self.gx_set.insert(item.clone());
        self.gx_set.insert(String::from_utf8(bio::alphabets::dna::revcomp(
            item.into_bytes()
        )).unwrap());
    }

}



