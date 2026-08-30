//! The record the trimmers operate on.
//!
//! fastp mutates a `Read` holding `std::string` members, so every trim stage
//! calls `substr`/`resize` and reallocates. Trimming only ever removes bases
//! from one end or the other, so we carry a window over the borrowed record
//! instead and never copy: the ingest hands the dedup table
//! `rec.bases()` once, at the end of the chain.

/// A window over one read's sequence and quality.
///
/// `qual` is empty for FASTA input; every quality-dependent stage checks
/// [`has_qual`](ReadRec::has_qual) and no-ops without it.
#[derive(Clone, Copy, Debug)]
pub struct ReadRec<'a> {
    seq: &'a [u8],
    qual: &'a [u8],
    start: usize,
    end: usize,
}

impl<'a> ReadRec<'a> {
    pub fn new(seq: &'a [u8], qual: &'a [u8]) -> Self {
        debug_assert!(qual.is_empty() || qual.len() == seq.len());
        ReadRec { seq, qual, start: 0, end: seq.len() }
    }

    /// Sequence with no quality, as FASTA input arrives. Test convenience.
    #[allow(dead_code)]
    pub fn from_seq(seq: &'a [u8]) -> Self {
        ReadRec::new(seq, b"")
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[inline]
    pub fn has_qual(&self) -> bool {
        !self.qual.is_empty()
    }

    #[inline]
    pub fn bases(&self) -> &'a [u8] {
        &self.seq[self.start..self.end]
    }

    #[inline]
    pub fn quals(&self) -> &'a [u8] {
        if self.qual.is_empty() {
            b""
        } else {
            &self.qual[self.start..self.end]
        }
    }

    /// Keep the first `n` bases, fastp's `Read::resize`.
    #[inline]
    pub fn resize(&mut self, n: usize) {
        self.end = self.start + n.min(self.len());
    }

    /// Drop `n` bases from the front, fastp's `erase(0, n)`.
    #[inline]
    pub fn trim_front(&mut self, n: usize) {
        self.start += n.min(self.len());
    }

    /// Offset of this window's first base within the original record. The
    /// overlap trimmer needs it as fastp's `frontTrimmed`.
    #[inline]
    pub fn front_trimmed(&self) -> usize {
        self.start
    }
}

#[cfg(test)]
#[path = "../tests/read.rs"]
mod tests;
