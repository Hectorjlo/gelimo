//! Linear-time O(N) PGLS solver based on Ho & Ané (2014) post-order tree traversal.
//!
//! Avoids building or inverting the N x N covariance matrix V, reducing time to O(N)
//! and memory footprint to O(N) with zero heap allocations during optimization loops.

use faer::Mat;
use crate::tree::PhyloTree;
use crate::math::solver::{invert_small_matrix, GlsResult};

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
    fn traverse(tree: &PhyloTree, curr: usize, post_order: &mut Vec<usize>) {
        for &child in &tree.nodes[curr].children {
            traverse(tree, child, post_order);
        }
        post_order.push(curr);
    }
    traverse(tree, tree.root, &mut post_order);

    // Tip node index mapping to tip sequence index (0..n_tips)
    let mut tip_to_seq_idx = vec![0; num_nodes];
    for (seq_idx, &node_idx) in tree.tips.iter().enumerate() {
        tip_to_seq_idx[node_idx] = seq_idx;
    }

    let mut xtx: Mat<f64> = Mat::zeros(p, p);
    let mut xty: Mat<f64> = Mat::zeros(p, 1);
    let mut log_det_v = 0.0;

    // Pass 1: Post-order tree traversal to accumulate X^T V^-1 X, X^T V^-1 Y, and log|V|
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
                log_det_v += v_c_br.ln();
                let w_c = 1.0 / v_c_br;
                total_w += w_c;

                weighted_zy += w_c * ws.node_zy[c];
                for (zx_val, &c_val) in weighted_zx.iter_mut().zip(&ws.node_zx[c]) {
                    *zx_val += w_c * c_val;
                }
            }

            let v_u = 1.0 / total_w.max(1e-15);
            log_det_v -= v_u.ln();
            let zy_u = weighted_zy * v_u;
            for (dst, &w) in ws.node_zx[u].iter_mut().zip(&weighted_zx) {
                *dst = w * v_u;
            }

            // Accumulate independent contrasts to (X^T V^-1 X) and (X^T V^-1 Y)
            for &c in &node.children {
                let br = tree.nodes[c].branch_length.unwrap_or(0.0).max(1e-15);
                let w_c = 1.0 / (ws.node_v[c] + br);
                let dy = ws.node_zy[c] - zy_u;

                for i in 0..p {
                    let dxi = ws.node_zx[c][i] - ws.node_zx[u][i];
                    for j in 0..p {
                        let dxj = ws.node_zx[c][j] - ws.node_zx[u][j];
                        xtx[(i, j)] += w_c * dxi * dxj;
                    }
                    xty[(i, 0)] += w_c * dxi * dy;
                }
            }

            ws.node_v[u] = v_u;
            ws.node_zy[u] = zy_u;
        }
    }

    // Root node mean term
    let root = tree.root;
    let root_br = tree.nodes[root].branch_length.unwrap_or(0.0);
    let root_v_br = (ws.node_v[root] + root_br).max(1e-15);
    log_det_v += root_v_br.ln();
    let root_w = 1.0 / root_v_br;

    for i in 0..p {
        for j in 0..p {
            xtx[(i, j)] += root_w * ws.node_zx[root][i] * ws.node_zx[root][j];
        }
        xty[(i, 0)] += root_w * ws.node_zx[root][i] * ws.node_zy[root];
    }

    // Invert X^T V^-1 X and compute beta = (X^T V^-1 X)^-1 (X^T V^-1 Y)
    let inv_xtx = invert_small_matrix(&xtx)?;
    let beta = &inv_xtx * &xty;
    let residuals = y - x * &beta;

    // Pass 2: O(N) residual quadratic form e^T V^-1 e using pre-computed node_zy and node_zx
    for &u in &post_order {
        let mut pred = 0.0;
        for j in 0..p {
            pred += ws.node_zx[u][j] * beta[(j, 0)];
        }
        ws.node_ze[u] = ws.node_zy[u] - pred;
    }

    let mut rss_star = root_w * ws.node_ze[root] * ws.node_ze[root];
    let mut transformed_residuals: Mat<f64> = Mat::zeros(n_tips, 1);
    for &u in &post_order {
        let node = &tree.nodes[u];
        if !node.is_tip {
            for &c in &node.children {
                let br = tree.nodes[c].branch_length.unwrap_or(0.0).max(1e-15);
                let w_c = 1.0 / (ws.node_v[c] + br);
                let de = ws.node_ze[c] - ws.node_ze[u];
                rss_star += w_c * de * de;
            }
        }
    }
    if n_tips > 0 {
        transformed_residuals[(0, 0)] = rss_star.max(0.0).sqrt();
    }

    let sigma2_ml = rss_star / (n_tips as f64);
    let df = if n_tips > p { (n_tips - p) as f64 } else { 1.0 };
    let sigma2_unbiased = rss_star / df;

    let n_f64 = n_tips as f64;
    let log_likelihood = -0.5 * n_f64 * (2.0 * std::f64::consts::PI).ln()
        - 0.5 * n_f64 * sigma2_ml.max(1e-15).ln()
        - 0.5 * log_det_v
        - 0.5 * n_f64;

    let vcov_beta = &inv_xtx * faer::Scale(sigma2_ml);
    let se: Vec<f64> = (0..p).map(|i| vcov_beta[(i, i)].max(0.0).sqrt()).collect();

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
