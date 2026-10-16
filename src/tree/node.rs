//! Dynamic phylogenetic node and tree representation with species name normalization.

use std::collections::HashSet;

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

    /// Prunes the tree, keeping only tips whose normalized name is in `species_to_keep`.
    /// Single-child internal nodes resulting from pruning are collapsed.
    pub fn prune(&self, species_to_keep: &HashSet<String>) -> Result<PhyloTree, String> {
        let keep_node_ids: HashSet<usize> = self
            .tips
            .iter()
            .copied()
            .filter(|&tip_id| {
                let tip_name = self.nodes[tip_id].name.as_deref().unwrap_or("");
                species_to_keep.contains(&normalize_species_name(tip_name))
            })
            .collect();

        if keep_node_ids.is_empty() {
            return Err("Pruning resulted in empty tree (0 tips remaining)".to_string());
        }

        fn prune_rec(
            tree: &PhyloTree,
            curr_idx: usize,
            keep_node_ids: &HashSet<usize>,
            nodes: &mut Vec<Node>,
            tips: &mut Vec<usize>,
        ) -> Option<usize> {
            let node = &tree.nodes[curr_idx];
            if node.is_tip {
                if !keep_node_ids.contains(&curr_idx) {
                    return None;
                }
                let id = nodes.len();
                nodes.push(Node {
                    id,
                    name: node.name.clone(),
                    branch_length: node.branch_length,
                    children: Vec::new(),
                    parent: None,
                    is_tip: true,
                });
                tips.push(id);
                return Some(id);
            }

            let kept_children: Vec<usize> = node
                .children
                .iter()
                .filter_map(|&c| prune_rec(tree, c, keep_node_ids, nodes, tips))
                .collect();

            match kept_children.len() {
                0 => None,
                1 => {
                    let child_id = kept_children[0];
                    let merged_len = node.branch_length.unwrap_or(0.0)
                        + nodes[child_id].branch_length.unwrap_or(0.0);
                    nodes[child_id].branch_length = Some(merged_len);
                    Some(child_id)
                }
                _ => {
                    let id = nodes.len();
                    for &c_id in &kept_children {
                        nodes[c_id].parent = Some(id);
                    }
                    nodes.push(Node {
                        id,
                        name: node.name.clone(),
                        branch_length: node.branch_length,
                        children: kept_children,
                        parent: None,
                        is_tip: false,
                    });
                    Some(id)
                }
            }
        }

        let mut nodes = Vec::new();
        let mut tips = Vec::new();
        let root = prune_rec(self, self.root, &keep_node_ids, &mut nodes, &mut tips)
            .ok_or_else(|| "Failed to prune tree root".to_string())?;

        Ok(PhyloTree { nodes, root, tips })
    }

    /// Serializes the phylogenetic tree back into standard Newick format string.
    pub fn to_newick(&self) -> String {
        let mut buf = String::new();
        self.format_node_newick(self.root, &mut buf);
        buf.push(';');
        buf
    }

    fn format_node_newick(&self, node_idx: usize, buf: &mut String) {
        let node = &self.nodes[node_idx];
        if !node.children.is_empty() {
            buf.push('(');
            for (i, &child_idx) in node.children.iter().enumerate() {
                if i > 0 {
                    buf.push(',');
                }
                self.format_node_newick(child_idx, buf);
            }
            buf.push(')');
        }
        if let Some(ref name) = node.name {
            buf.push_str(name);
        }
        if let Some(bl) = node.branch_length {
            use std::fmt::Write;
            let _ = write!(buf, ":{}", bl);
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::parser::parse_newick;

    #[test]
    fn test_tree_pruning_and_collapsing() {
        let input = "((SpA:1.0,SpB:1.0):1.0,SpC:2.0);";
        let tree = parse_newick(input).unwrap();

        let mut keep = HashSet::new();
        keep.insert("SpB".to_string());
        keep.insert("SpC".to_string());

        let pruned = tree.prune(&keep).unwrap();

        assert_eq!(pruned.tips.len(), 2);
        assert_eq!(pruned.normalized_tip_names(), vec!["SpB", "SpC"]);

        // SpB branch length should be collapsed: 1.0 (from SpB) + 1.0 (from internal node) = 2.0
        let spb_node = &pruned.nodes[pruned.tips[0]];
        assert_eq!(spb_node.branch_length, Some(2.0));
    }

    #[test]
    fn test_tree_to_newick_roundtrip() {
        let input = "((SpA:1.5,SpB:2.5):0.5,SpC:3);";
        let tree = parse_newick(input).unwrap();
        let newick_out = tree.to_newick();
        
        let roundtrip_tree = parse_newick(&newick_out).unwrap();
        assert_eq!(roundtrip_tree.tips.len(), 3);
        assert_eq!(roundtrip_tree.normalized_tip_names(), vec!["SpA", "SpB", "SpC"]);
    }
}
