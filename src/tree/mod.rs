//! Tree module for phylogenetic node structures and Newick parsing.

pub mod node;
pub mod parser;

pub use crate::math::pagel::apply_pagel_lambda_tree;
pub use node::{normalize_species_name, Node, PhyloTree};
pub use parser::parse_newick;
