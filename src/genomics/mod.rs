//! Comparative genomics and high-throughput orthogroup analysis module.

pub mod engine;
pub mod filter;
pub mod multiple_testing;
pub mod ortho_reader;
pub mod pipeline;

pub use engine::{run_ortho_pgls_parallel, OrthoPglsConfig, OrthoPglsResult};
pub use filter::{
    filter_orthogroups, FilteredOrthoData, OrthoFilterOptions, OrthoFilterReport,
};
pub use multiple_testing::{
    adjust_p_values, export_ortho_results_tsv, AnnotatedOrthoResult,
};
pub use ortho_reader::{
    align_tripartite, read_orthogroup_gene_counts, OrthoGeneCountTable,
    TripartiteComparativeData, TripartiteReport,
};
pub use pipeline::{run_ortho_pipeline, OrthoPipelineConfig, OrthoPipelineSummary};
