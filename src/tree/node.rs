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
    //  Clean recursive node pruning and 1-child node collapsing
    pub fn prune(&self, species_to_keep: &HashSet<String>) -> Result<PhyloTree, String> {
        let mut keep_node_ids = HashSet::new();
        for &tip_id in &self.tips {
            let tip_name = self.nodes[tip_id].name.as_deref().unwrap_or("");
            let norm = normalize_species_name(tip_name);
            if species_to_keep.contains(&norm) {
                keep_node_ids.insert(tip_id);
            }
        }

        if keep_node_ids.is_empty() {
            return Err("Pruning resulted in empty tree (0 tips remaining)".to_string());
        }

        struct SimpleNode {
            name: Option<String>,
            branch_length: Option<f64>,
            children: Vec<SimpleNode>,
        }

        fn build_simple(tree: &PhyloTree, curr_idx: usize, keep_node_ids: &HashSet<usize>) -> Option<SimpleNode> {
            let node = &tree.nodes[curr_idx];
            if node.is_tip {
                if keep_node_ids.contains(&curr_idx) {
                    Some(SimpleNode {
                        name: node.name.clone(),
                        branch_length: node.branch_length,
                        children: Vec::new(),
                    })
                } else {
                    None
                }
            } else {
                let mut kept_children = Vec::new();
                for &child_idx in &node.children {
                    if let Some(child_simple) = build_simple(tree, child_idx, keep_node_ids) {
                        kept_children.push(child_simple);
                    }
                }

                if kept_children.is_empty() {
                    None
                } else if kept_children.len() == 1 {
                    let mut child = kept_children.pop().unwrap();
                    let my_len = node.branch_length.unwrap_or(0.0);
                    let child_len = child.branch_length.unwrap_or(0.0);
                    child.branch_length = Some(my_len + child_len);
                    Some(child)
                } else {
                    Some(SimpleNode {
                        name: node.name.clone(),
                        branch_length: node.branch_length,
                        children: kept_children,
                    })
                }
            }
        }

        let root_simple = build_simple(self, self.root, &keep_node_ids)
            .ok_or_else(|| "Failed to prune tree root".to_string())?;

        let mut nodes = Vec::new();
        let mut tips = Vec::new();

        fn flatten(
            snode: SimpleNode,
            parent: Option<usize>,
            nodes: &mut Vec<Node>,
            tips: &mut Vec<usize>,
        ) -> usize {
            let id = nodes.len();
            let is_tip = snode.children.is_empty();

            nodes.push(Node {
                id,
                name: snode.name,
                branch_length: snode.branch_length,
                children: Vec::new(),
                parent,
                is_tip,
            });

            if is_tip {
                tips.push(id);
            } else {
                let child_ids: Vec<usize> = snode
                    .children
                    .into_iter()
                    .map(|child| flatten(child, Some(id), nodes, tips))
                    .collect();
                nodes[id].children = child_ids;
            }

            id
        }

        let root_idx = flatten(root_simple, None, &mut nodes, &mut tips);

        Ok(PhyloTree {
            nodes,
            root: root_idx,
            tips,
        })
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
}
