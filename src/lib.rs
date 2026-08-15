//! `gelimo`: Engine for Generalized Phylogenetic Linear Models in Rust.
//!
//! Compatible with input and output specifications of R's `caper` package.

pub mod data;
pub mod genomics;
pub mod math;
pub mod opt;
pub mod pipeline;
pub mod stats;
pub mod tree;

pub use genomics::{
    align_tripartite, filter_orthogroups, read_orthogroup_gene_counts, FilteredOrthoData,
    OrthoFilterOptions, OrthoFilterReport, OrthoGeneCountTable, TripartiteComparativeData,
    TripartiteReport,
};
pub use opt::{grid_search_lambda_parallel, optimize_pagel_lambda, LambdaOptResult};
pub use pipeline::{run_pgls_pipeline, PglsConfig, PglsSummary};
