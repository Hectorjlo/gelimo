//! `gelimo`: Engine for Generalized Phylogenetic Linear Models in Rust.
//!
//! Compatible with input and output specifications of R's `caper` package.

pub mod data;
pub mod math;
pub mod pipeline;
pub mod stats;
pub mod tree;

pub use pipeline::{run_pgls_pipeline, PglsConfig, PglsSummary};
