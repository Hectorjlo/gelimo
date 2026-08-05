//! Dynamic phylogenetic node and tree representation with species name normalization.

/// Represents a single node in a phylogenetic tree.
#[derive(Debug, Clone, PartialEq)]
pub struct Node {
    pub id: usize,
    pub name: Option<String>,
    pub branch_length: Option<f64>,
    pub children: Vec<usize>,
    pub parent: Option<usize>,
    pub is_tip: bool,
}

/// Represents a complete phylogenetic tree.
#[derive(Debug, Clone, PartialEq)]
pub struct PhyloTree {
    pub nodes: Vec<Node>,
    pub root: usize,
    pub tips: Vec<usize>,
}

impl PhyloTree {
    /// Returns raw tip names in the tree.
    pub fn tip_names(&self) -> Vec<String> {
        self.tips
            .iter()
            .map(|&idx| {
                self.nodes[idx]
                    .name
                    .clone()
                    .unwrap_or_else(|| format!("node_{}", idx))
            })
            .collect()
    }

    /// Returns normalized tip names (spaces replaced by underscores, common biological suffixes removed).
    pub fn normalized_tip_names(&self) -> Vec<String> {
        self.tip_names()
            .iter()
            .map(|name| normalize_species_name(name))
            .collect()
    }
}

/// Normalizes a species name by replacing spaces with underscores and stripping common suffixes.
pub fn normalize_species_name(raw: &str) -> String {
    let mut name = raw.trim().replace(' ', "_");
    
    // List of biological suffixes to strip
    let suffixes = [".faa", ".fasta", "_protein", "_pep", "_cds", "_genome", ".fa"];
    for suffix in suffixes {
        if name.ends_with(suffix) {
            name = name[..name.len() - suffix.len()].to_string();
        }
    }
    name
}
