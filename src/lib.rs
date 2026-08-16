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
    adjust_p_values, align_tripartite, export_ortho_results_tsv, filter_orthogroups,
    read_orthogroup_gene_counts, run_ortho_pgls_parallel, run_ortho_pipeline,
    AnnotatedOrthoResult, FilteredOrthoData, OrthoFilterOptions, OrthoFilterReport,
    OrthoGeneCountTable, OrthoModelMode, OrthoPipelineConfig, OrthoPipelineSummary, OrthoPglsConfig,
    OrthoPglsResult, TripartiteComparativeData, TripartiteReport,
};
pub use opt::{grid_search_lambda_parallel, optimize_pagel_lambda, LambdaOptResult};
pub use pipeline::{run_pgls_pipeline, PglsConfig, PglsSummary};
