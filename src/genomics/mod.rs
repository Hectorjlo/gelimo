//! Comparative genomics and high-throughput orthogroup analysis module.

pub mod filter;
pub mod ortho_reader;

pub use filter::{
    filter_orthogroups, FilteredOrthoData, OrthoFilterOptions, OrthoFilterReport,
};
pub use ortho_reader::{
    align_tripartite, read_orthogroup_gene_counts, OrthoGeneCountTable,
    TripartiteComparativeData, TripartiteReport,
};
