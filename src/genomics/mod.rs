//! Comparative genomics and high-throughput orthogroup analysis module.

pub mod engine;
pub mod filter;
pub mod ortho_reader;

pub use engine::{run_ortho_pgls_parallel, OrthoPglsConfig, OrthoPglsResult};
pub use filter::{
    filter_orthogroups, FilteredOrthoData, OrthoFilterOptions, OrthoFilterReport,
};
pub use ortho_reader::{
    align_tripartite, read_orthogroup_gene_counts, OrthoGeneCountTable,
    TripartiteComparativeData, TripartiteReport,
};
