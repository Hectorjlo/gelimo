//! Pagel's evolutionary parameter transformations (lambda, kappa, delta).

use faer::Mat;
use crate::tree::PhyloTree;

/// Applies Pagel's $\lambda$ transformation to a VCV matrix $V$.
///
/// $V(\lambda)_{i,j} = V_{i,i}$ if $i = j$, and $\lambda \cdot V_{i,j}$ if $i \neq j$.
// ponytail: elementwise transformation preserving original diagonal
pub fn apply_pagel_lambda(vcv: &Mat<f64>, lambda: f64) -> Mat<f64> {
    let n = vcv.nrows();
    let mut v_lambda = vcv.clone();

    for i in 0..n {
        for j in 0..n {
            if i != j {
                v_lambda[(i, j)] *= lambda;
            }
        }
    }

    v_lambda
}

/// Applies Pagel's $\kappa$ transformation to a `PhyloTree`, scaling branch lengths by $b^\kappa$.
// ponytail: pure functional branch length exponentiation
pub fn apply_pagel_kappa(tree: &PhyloTree, kappa: f64) -> PhyloTree {
    let mut new_tree = tree.clone();
    for node in &mut new_tree.nodes {
        if let Some(b) = node.branch_length {
            node.branch_length = Some(b.powf(kappa));
        }
    }
    new_tree
}

/// Applies Pagel's $\delta$ transformation to a `PhyloTree`, scaling node heights $h$ by $h^\delta$.
// ponytail: node height scaling and branch length delta recalculation
pub fn apply_pagel_delta(tree: &PhyloTree, delta: f64) -> PhyloTree {
    let n = tree.nodes.len();
    let mut heights = vec![0.0; n];

    // Compute root-to-node heights
    for i in 0..n {
        let mut curr = i;
        let mut sum = 0.0;
        while let Some(parent_idx) = tree.nodes[curr].parent {
            sum += tree.nodes[curr].branch_length.unwrap_or(0.0);
            curr = parent_idx;
        }
        heights[i] = sum;
    }

    let mut new_tree = tree.clone();
    for i in 0..n {
        if let Some(p) = tree.nodes[i].parent {
            let h_child = heights[i].powf(delta);
            let h_parent = heights[p].powf(delta);
            new_tree.nodes[i].branch_length = Some((h_child - h_parent).max(0.0));
        }
    }

    new_tree
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::parser::parse_newick;
    use approx::assert_relative_eq;

    #[test]
    fn test_apply_pagel_lambda() {
        let mut v = Mat::zeros(2, 2);
        v[(0, 0)] = 2.0;
        v[(0, 1)] = 1.0;
        v[(1, 0)] = 1.0;
        v[(1, 1)] = 2.0;

        let v_lam = apply_pagel_lambda(&v, 0.5);

        assert_relative_eq!(v_lam[(0, 0)], 2.0, epsilon = 1e-12);
        assert_relative_eq!(v_lam[(0, 1)], 0.5, epsilon = 1e-12);
        assert_relative_eq!(v_lam[(1, 0)], 0.5, epsilon = 1e-12);
        assert_relative_eq!(v_lam[(1, 1)], 2.0, epsilon = 1e-12);
    }

    #[test]
fn test_apply_pagel_lambda_asymmetric_3x3() {
    // Matriz de covarianza de 3 especies
    let mut v = Mat::zeros(3, 3);
    v[(0, 0)] = 3.0; v[(0, 1)] = 1.5; v[(0, 2)] = 0.5;
    v[(1, 0)] = 1.5; v[(1, 1)] = 2.5; v[(1, 2)] = 0.5;
    v[(2, 0)] = 0.5; v[(2, 1)] = 0.5; v[(2, 2)] = 4.0;

    let lambda = 0.8;
    let v_lam = apply_pagel_lambda(&v, lambda);

    // 1. Comprobar que la diagonal NO cambió
    assert_relative_eq!(v_lam[(0, 0)], 3.0, epsilon = 1e-12);
    assert_relative_eq!(v_lam[(1, 1)], 2.5, epsilon = 1e-12);
    assert_relative_eq!(v_lam[(2, 2)], 4.0, epsilon = 1e-12);

    // 2. Comprobar fuera de la diagonal (escalados por lambda)
    assert_relative_eq!(v_lam[(0, 1)], 1.2, epsilon = 1e-12);
    assert_relative_eq!(v_lam[(1, 0)], 1.2, epsilon = 1e-12);

    assert_relative_eq!(v_lam[(0, 2)], 0.4, epsilon = 1e-12);
    assert_relative_eq!(v_lam[(2, 0)], 0.4, epsilon = 1e-12);

    assert_relative_eq!(v_lam[(1, 2)], 0.4, epsilon = 1e-12);
    assert_relative_eq!(v_lam[(2, 1)], 0.4, epsilon = 1e-12);

    // 3. Verificar que la matriz original sigue intacta (inmutabilidad)
    assert_relative_eq!(v[(0, 1)], 1.5, epsilon = 1e-12);
    }

    #[test]
    fn test_apply_pagel_kappa() {
        let newick = "((SpA:4.0,SpB:4.0):1.0,SpC:9.0);";
        let tree = parse_newick(newick).unwrap();
        // kappa = 0.5 -> branch lengths become sqrt(b): 4.0 -> 2.0, 1.0 -> 1.0, 9.0 -> 3.0
        let kappa_tree = apply_pagel_kappa(&tree, 0.5);

        let spa_node = &kappa_tree.nodes[kappa_tree.tips[0]];
        assert_relative_eq!(spa_node.branch_length.unwrap(), 2.0, epsilon = 1e-12);
    }
}
