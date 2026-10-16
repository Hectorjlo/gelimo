//! `gelimo`: Engine for Generalized Phylogenetic Linear Models in Rust.
//!
//! Compatible with input and output specifications of R's `caper` package.

pub mod cli;
pub mod data;
pub mod genomics;
pub mod math;
pub mod opt;
pub mod pipeline;
#[cfg(feature = "r-bindings")]
pub mod r_bindings;
pub mod stats;
pub mod tree;

pub use data::{ComparativeData, ComparativeDataOptions, PruneReport, TraitTable};
pub use genomics::{
    adjust_p_values, align_tripartite, export_ortho_results_tsv, filter_orthogroups,
    read_orthogroup_gene_counts, run_ortho_pgls_parallel, run_ortho_pipeline,
    AnnotatedOrthoResult, FilteredOrthoData, OrthoFilterOptions, OrthoFilterReport,
    OrthoGeneCountTable, OrthoModelMode, OrthoPipelineConfig, OrthoPipelineSummary, OrthoPglsConfig,
    OrthoPglsResult, TripartiteComparativeData, TripartiteReport,
};
pub use math::{apply_pagel_lambda, apply_pagel_lambda_tree, compute_vcv, solve_gls, GlsResult};
pub use opt::{optimize_pagel_lambda, LambdaOptResult};
pub use pipeline::{run_pgls_pipeline, PglsConfig, PglsSummary};
pub use tree::{normalize_species_name, parse_newick, Node, PhyloTree};
