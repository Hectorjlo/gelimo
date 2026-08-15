//! Comparative genomics and high-throughput orthogroup analysis module.

pub mod ortho_reader;

pub use ortho_reader::{
    read_orthogroup_gene_counts, align_tripartite, OrthoGeneCountTable,
    TripartiteComparativeData, TripartiteReport,
};
