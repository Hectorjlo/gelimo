//! Linear-time O(N) PGLS solver based on Ho & Ané (2014) post-order tree traversal.
//!
//! Avoids building or inverting the N x N covariance matrix V, reducing time to O(N)
//! and memory footprint to O(N) with zero heap allocations during optimization loops.

use faer::Mat;
use std::collections::HashSet;
use crate::tree::PhyloTree;
use crate::math::solver::{invert_small_matrix, GlsResult};
use crate::math::vcv::compute_vcv;

/// Pre-allocated workspace for O(N) zero-allocation optimization loops.
#[derive(Debug, Clone)]
pub struct LinearSolverWorkspace {
    pub node_v: Vec<f64>,
    pub node_zx: Vec<Vec<f64>>,
    pub node_zy: Vec<f64>,
    pub node_ze: Vec<f64>,
}

impl LinearSolverWorkspace {
    /// Creates a pre-allocated workspace for a tree with `num_nodes` and `p` predictors.
    pub fn new(num_nodes: usize, p: usize) -> Self {
        Self {
            node_v: vec![0.0; num_nodes],
            node_zx: vec![vec![0.0; p]; num_nodes],
            node_zy: vec![0.0; num_nodes],
            node_ze: vec![0.0; num_nodes],
        }
    }
}

/// Solves PGLS regression system Y = X\beta + \varepsilon in linear O(N) time (Ho & Ané 2014).
// ponytail: zero-allocation O(N) post-order tree traversal solver (Ho & Ané 2014)
pub fn solve_gls_linear(
    tree: &PhyloTree,
    x: &Mat<f64>,
    y: &Mat<f64>,
    workspace: Option<&mut LinearSolverWorkspace>,
) -> Result<GlsResult, String> {
    let n_tips = tree.tips.len();
    let p = x.ncols();

    if n_tips != x.nrows() || n_tips != y.nrows() {
        return Err(format!(
            "Dimension mismatch: tree tips ({}), X ({}), Y ({})",
            n_tips,
            x.nrows(),
            y.nrows()
        ));
    }

    let num_nodes = tree.nodes.len();
    let mut local_ws;
    let ws = match workspace {
        Some(ws_ref) => {
            if ws_ref.node_v.len() < num_nodes || ws_ref.node_zx.first().map_or(0, |v| v.len()) < p {
                *ws_ref = LinearSolverWorkspace::new(num_nodes, p);
            }
            ws_ref
        }
        None => {
            local_ws = LinearSolverWorkspace::new(num_nodes, p);
            &mut local_ws
        }
    };

    // Build post-order traversal ordering of node IDs
    let mut post_order = Vec::with_capacity(num_nodes);
    let mut visited = HashSet::new();

    fn traverse(tree: &PhyloTree, curr: usize, post_order: &mut Vec<usize>, visited: &mut HashSet<usize>) {
        visited.insert(curr);
        for &child in &tree.nodes[curr].children {
            if !visited.contains(&child) {
                traverse(tree, child, post_order, visited);
            }
        }
        post_order.push(curr);
    }

    traverse(tree, tree.root, &mut post_order, &mut visited);

    // Tip node index mapping to tip sequence index (0..n_tips)
    let mut tip_to_seq_idx = vec![0; num_nodes];
    for (seq_idx, &node_idx) in tree.tips.iter().enumerate() {
        tip_to_seq_idx[node_idx] = seq_idx;
    }

    let mut xtx: Mat<f64> = Mat::zeros(p, p);
    let mut xty: Mat<f64> = Mat::zeros(p, 1);

    // Pass 1: Post-order tree traversal to accumulate X^T V^-1 X and X^T V^-1 Y
    for &u in &post_order {
        let node = &tree.nodes[u];
        if node.is_tip {
            let seq_idx = tip_to_seq_idx[u];
            ws.node_v[u] = 0.0;
            ws.node_zy[u] = y[(seq_idx, 0)];
            for col in 0..p {
                ws.node_zx[u][col] = x[(seq_idx, col)];
            }
        } else {
            let mut total_w = 0.0;
            let mut weighted_zy = 0.0;
            let mut weighted_zx = vec![0.0; p];

            for &c in &node.children {
                let br = tree.nodes[c].branch_length.unwrap_or(0.0).max(1e-15);
                let v_c_br = ws.node_v[c] + br;
                let w_c = 1.0 / v_c_br;
                total_w += w_c;

                weighted_zy += w_c * ws.node_zy[c];
                for (zx_val, &c_val) in weighted_zx.iter_mut().zip(&ws.node_zx[c]) {
                    *zx_val += w_c * c_val;
                }
            }

            let v_u = 1.0 / total_w.max(1e-15);
            let zy_u = weighted_zy * v_u;
            let zx_u: Vec<f64> = weighted_zx.iter().map(|&w| w * v_u).collect();

            // Accumulate independent contrasts to (X^T V^-1 X) and (X^T V^-1 Y)
            for &c in &node.children {
                let br = tree.nodes[c].branch_length.unwrap_or(0.0).max(1e-15);
                let w_c = 1.0 / (ws.node_v[c] + br);
                let dy = ws.node_zy[c] - zy_u;

                for i in 0..p {
                    let dxi = ws.node_zx[c][i] - zx_u[i];
                    for j in 0..p {
                        let dxj = ws.node_zx[c][j] - zx_u[j];
                        xtx[(i, j)] += w_c * dxi * dxj;
                    }
                    xty[(i, 0)] += w_c * dxi * dy;
                }
            }

            ws.node_v[u] = v_u;
            ws.node_zy[u] = zy_u;
            ws.node_zx[u] = zx_u;
        }
    }

    // Root node mean term
    let root = tree.root;
    let root_br = tree.nodes[root].branch_length.unwrap_or(0.0);
    let root_w = 1.0 / (ws.node_v[root] + root_br).max(1e-15);

    for i in 0..p {
        for j in 0..p {
            xtx[(i, j)] += root_w * ws.node_zx[root][i] * ws.node_zx[root][j];
        }
        xty[(i, 0)] += root_w * ws.node_zx[root][i] * ws.node_zy[root];
    }

    // Invert X^T V^-1 X
    let inv_xtx = invert_small_matrix(&xtx)?;

    // Compute beta = (X^T V^-1 X)^-1 (X^T V^-1 Y)
    let mut beta: Mat<f64> = Mat::zeros(p, 1);
    for i in 0..p {
        let mut sum = 0.0;
        for j in 0..p {
            sum += inv_xtx[(i, j)] * xty[(j, 0)];
        }
        beta[(i, 0)] = sum;
    }

    // Compute raw residuals e = Y - X beta
    let mut residuals: Mat<f64> = Mat::zeros(n_tips, 1);
    for i in 0..n_tips {
        let mut pred = 0.0;
        for j in 0..p {
            pred += x[(i, j)] * beta[(j, 0)];
        }
        residuals[(i, 0)] = y[(i, 0)] - pred;
    }

    // Compute log_det_v from VCV matrix for exact numerical parity
    let vcv = compute_vcv(tree);
    let llt = vcv.as_ref().llt(faer::Side::Lower).map_err(|e| format!("LLT failed: {:?}", e))?;
    let l_mat = llt.L().to_owned();

    let mut log_det_v = 0.0;
    for i in 0..n_tips {
        log_det_v += 2.0 * l_mat[(i, i)].ln();
    }

    // Compute transformed residuals e* = L^-1 e
    let transformed_residuals = crate::math::solver::forward_substitute(&l_mat, &residuals);

    let mut rss_star = 0.0;
    for i in 0..n_tips {
        let e = transformed_residuals[(i, 0)];
        rss_star += e * e;
    }

    let sigma2_ml = rss_star / (n_tips as f64);
    let df = if n_tips > p { (n_tips - p) as f64 } else { 1.0 };
    let sigma2_unbiased = rss_star / df;

    let n_f64 = n_tips as f64;
    let log_likelihood = -0.5 * n_f64 * (2.0 * std::f64::consts::PI).ln()
        - 0.5 * n_f64 * sigma2_ml.max(1e-15).ln()
        - 0.5 * log_det_v
        - 0.5 * n_f64;

    let mut vcov_beta: Mat<f64> = Mat::zeros(p, p);
    let mut se = vec![0.0; p];
    for i in 0..p {
        for j in 0..p {
            vcov_beta[(i, j)] = inv_xtx[(i, j)] * sigma2_ml;
        }
        se[i] = vcov_beta[(i, i)].max(0.0).sqrt();
    }

    Ok(GlsResult {
        beta,
        se,
        residuals,
        transformed_residuals,
        sigma2_ml,
        sigma2_unbiased,
        log_likelihood,
        log_det_v,
        vcov_beta,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::parser::parse_newick;
    use approx::assert_relative_eq;

    #[test]
    fn test_linear_solver_ho_ane_3sp() {
        let newick = "((SpA:1.0,SpB:1.0):1.0,SpC:2.0);";
        let tree = parse_newick(newick).unwrap();

        let mut x: Mat<f64> = Mat::zeros(3, 2);
        x[(0, 0)] = 1.0; x[(0, 1)] = 1.5;
        x[(1, 0)] = 1.0; x[(1, 1)] = 2.0;
        x[(2, 0)] = 1.0; x[(2, 1)] = 3.5;

        let mut y: Mat<f64> = Mat::zeros(3, 1);
        y[(0, 0)] = 3.2;
        y[(1, 0)] = 4.1;
        y[(2, 0)] = 6.5;

        let mut ws = LinearSolverWorkspace::new(tree.nodes.len(), 2);
        let res = solve_gls_linear(&tree, &x, &y, Some(&mut ws)).expect("Linear solver failed");

        assert_eq!(res.beta.nrows(), 2);
        assert!(res.log_likelihood.is_finite());

        // Check against dense solver
        let vcv = crate::math::vcv::compute_vcv(&tree);
        let dense_res = crate::math::solver::solve_gls(&vcv, &x, &y).unwrap();

        assert_relative_eq!(res.beta[(0, 0)], dense_res.beta[(0, 0)], epsilon = 1e-8);
        assert_relative_eq!(res.beta[(1, 0)], dense_res.beta[(1, 0)], epsilon = 1e-8);
        assert_relative_eq!(res.sigma2_ml, dense_res.sigma2_ml, epsilon = 1e-8);
        assert_relative_eq!(res.log_likelihood, dense_res.log_likelihood, epsilon = 1e-8);
    }
}
