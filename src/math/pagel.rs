//! Pagel's evolutionary parameter transformations (lambda).

use faer::Mat;
use crate::tree::PhyloTree;

/// Applies Pagel's $\lambda$ transformation to a VCV matrix $V$.
///
/// $V(\lambda)_{i,j} = V_{i,i}$ if $i = j$, and $\lambda \cdot V_{i,j}$ if $i \neq j$.
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

/// Applies Pagel's $\lambda$ transformation to a `PhyloTree` by scaling internal branch lengths.
pub fn apply_pagel_lambda_tree(tree: &PhyloTree, lambda: f64) -> PhyloTree {
    let mut new_tree = tree.clone();
    for node in &mut new_tree.nodes {
        if !node.is_tip && node.id != tree.root && node.branch_length.is_some() {
            node.branch_length = node.branch_length.map(|b| b * lambda);
        }
    }
    new_tree
}

#[cfg(test)]
mod tests {
    use super::*;
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
        let mut v = Mat::zeros(3, 3);
        v[(0, 0)] = 3.0; v[(0, 1)] = 1.5; v[(0, 2)] = 0.5;
        v[(1, 0)] = 1.5; v[(1, 1)] = 2.5; v[(1, 2)] = 0.5;
        v[(2, 0)] = 0.5; v[(2, 1)] = 0.5; v[(2, 2)] = 4.0;

        let lambda = 0.8;
        let v_lam = apply_pagel_lambda(&v, lambda);

        assert_relative_eq!(v_lam[(0, 0)], 3.0, epsilon = 1e-12);
        assert_relative_eq!(v_lam[(1, 1)], 2.5, epsilon = 1e-12);
        assert_relative_eq!(v_lam[(2, 2)], 4.0, epsilon = 1e-12);

        assert_relative_eq!(v_lam[(0, 1)], 1.2, epsilon = 1e-12);
        assert_relative_eq!(v_lam[(1, 0)], 1.2, epsilon = 1e-12);

        assert_relative_eq!(v_lam[(0, 2)], 0.4, epsilon = 1e-12);
        assert_relative_eq!(v_lam[(2, 0)], 0.4, epsilon = 1e-12);

        assert_relative_eq!(v_lam[(1, 2)], 0.4, epsilon = 1e-12);
        assert_relative_eq!(v_lam[(2, 1)], 0.4, epsilon = 1e-12);

        assert_relative_eq!(v[(0, 1)], 1.5, epsilon = 1e-12);
    }
}
