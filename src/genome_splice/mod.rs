//! exonfill in Rust: anchored exon fill, splice-site refinement and pseudogene
//! calls, ported from the C version built inside BATH.

pub mod align;
#[cfg(target_arch = "x86_64")]
pub mod avx;
pub mod chain;
pub mod gate;
pub mod gff;
pub mod hmm;
pub mod junction;
pub mod module;
pub mod orf;
pub mod pseudo;
pub mod refine;
pub mod run;
pub mod score;
pub mod sites;
pub mod splice;
