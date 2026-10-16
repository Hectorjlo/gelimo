//! Variance-Covariance Matrix (VCV) generation from phylogenetic trees.

use faer::Mat;
use std::collections::HashSet;
use crate::tree::PhyloTree;

/// Computes root-to-node path distance for all nodes in the tree.
fn compute_root_distances(tree: &PhyloTree) -> Vec<f64> {
    let mut dists = vec![0.0; tree.nodes.len()];
    for (i, d) in dists.iter_mut().enumerate() {
        let mut curr = i;
        let mut sum = 0.0;
        while let Some(parent_idx) = tree.nodes[curr].parent {
            sum += tree.nodes[curr].branch_length.unwrap_or(0.0);
            curr = parent_idx;
        }
        *d = sum;
    }
    dists
}

/// Finds the Most Recent Common Ancestor (MRCA) node index between two node indices.
fn find_mrca(tree: &PhyloTree, mut u: usize, mut v: usize) -> usize {
    let mut ancestors_u = HashSet::new();
    ancestors_u.insert(u);

    while let Some(p) = tree.nodes[u].parent {
        ancestors_u.insert(p);
        u = p;
    }

    while !ancestors_u.contains(&v) {
        if let Some(p) = tree.nodes[v].parent {
            v = p;
        } else {
            break;
        }
    }
    v
}

/// Computes the VCV matrix $V$ from a given `PhyloTree` for its tips.
pub fn compute_vcv(tree: &PhyloTree) -> Mat<f64> {
    let n = tree.tips.len();
    let root_dists = compute_root_distances(tree);
    let mut vcv = Mat::zeros(n, n);

    for i in 0..n {
        let tip_i_node = tree.tips[i];
        vcv[(i, i)] = root_dists[tip_i_node];

        for j in 0..i {
            let tip_j_node = tree.tips[j];
            let mrca_node = find_mrca(tree, tip_i_node, tip_j_node);
            let cov = root_dists[mrca_node];

            vcv[(i, j)] = cov;
            vcv[(j, i)] = cov;
        }
    }

    vcv
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::parser::parse_newick;
    use approx::assert_relative_eq;

    #[test]
    fn test_compute_vcv_3sp() {
        let newick = "((SpA:1.0,SpB:1.0):1.0,SpC:2.0);";
        let tree = parse_newick(newick).expect("Failed to parse Newick");
        let vcv = compute_vcv(&tree);

        assert_eq!(vcv.nrows(), 3);
        assert_eq!(vcv.ncols(), 3);

        // V = [[2.0, 1.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
        assert_relative_eq!(vcv[(0, 0)], 2.0, epsilon = 1e-12);
        assert_relative_eq!(vcv[(0, 1)], 1.0, epsilon = 1e-12);
        assert_relative_eq!(vcv[(0, 2)], 0.0, epsilon = 1e-12);

        assert_relative_eq!(vcv[(1, 0)], 1.0, epsilon = 1e-12);
        assert_relative_eq!(vcv[(1, 1)], 2.0, epsilon = 1e-12);
        assert_relative_eq!(vcv[(1, 2)], 0.0, epsilon = 1e-12);

        assert_relative_eq!(vcv[(2, 0)], 0.0, epsilon = 1e-12);
        assert_relative_eq!(vcv[(2, 1)], 0.0, epsilon = 1e-12);
        assert_relative_eq!(vcv[(2, 2)], 2.0, epsilon = 1e-12);
    }
}
