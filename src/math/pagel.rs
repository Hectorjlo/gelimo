//! Pagel's evolutionary parameter transformations (lambda, kappa).

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

/// Applies Pagel's $\lambda$ transformation to a `PhyloTree` by scaling internal branch lengths.
// ponytail: tree branch length lambda transformation
pub fn apply_pagel_lambda_tree(tree: &PhyloTree, lambda: f64) -> PhyloTree {
    let mut new_tree = tree.clone();
    for node in &mut new_tree.nodes {
        if !node.is_tip && node.id != tree.root && node.branch_length.is_some() {
            node.branch_length = node.branch_length.map(|b| b * lambda);
        }
    }
    new_tree
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
        let kappa_tree = apply_pagel_kappa(&tree, 0.5);

        let spa_node = &kappa_tree.nodes[kappa_tree.tips[0]];
        assert_relative_eq!(spa_node.branch_length.unwrap(), 2.0, epsilon = 1e-12);
    }
}
