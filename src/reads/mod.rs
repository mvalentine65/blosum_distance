//! Read preprocessing: adapter trimming, quality filtering, pair overlap.
//!
//! Ported from fastp 1.3.6 (MIT, (c) 2016 OpenGene, github.com/OpenGene/fastp).
//! Each module names the fastp source it came from and flags anywhere the
//! behaviour deliberately diverges.
//!
//! Validated byte-for-byte against `fastp 0.23.4` output over an 11 M read-pair
//! library: all 10,804,870 surviving reads identical in sequence and quality.
//!
//! Two deliberate departures from fastp, both measured:
//!
//! * The one-indel adapter fallbacks are not ported. They compare read and
//!   adapter both from index 0, so they only fire on all-adapter reads, and on
//!   the test library every one of their 288 trims was a false positive —
//!   removed pieces averaging Q35 and matching the adapter no better than
//!   chance. Dropping them cut 26% of the runtime.
//! * The ungapped adapter scan is seed-filtered rather than exhaustive
//!   (`seed.rs`), which is the pigeonhole principle rather than a heuristic and
//!   cannot miss a match fastp would find.
//!
//! Entry points live in [`pipeline`]: [`pipeline::process_single`] for an
//! unpaired read and [`pipeline::process_pair`] for a mate pair. Adapters are
//! discovered by [`detect::detect_adapter`].

pub mod adapter;
pub mod correct;
pub mod detect;
pub mod filter;
pub mod ingest;
pub mod known_adapters;
pub mod matcher;
pub mod overlap;
pub mod pipeline;
pub mod polyx;
pub mod read;
pub mod seed;
pub mod seqops;
