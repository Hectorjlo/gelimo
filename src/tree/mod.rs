//! Tree module for phylogenetic node structures and Newick parsing.

pub mod node;
pub mod parser;

pub use node::{Node, PhyloTree};
pub use parser::parse_newick;
