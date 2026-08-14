//! GLS linear system solver using Cholesky decomposition.

use faer::linalg::solvers::DenseSolveCore;
use faer::{Mat, Side};

/// Result container for GLS / PGLS regression solution.
#[derive(Debug, Clone)]
pub struct GlsResult {
    /// Estimated regression coefficients $\hat{\beta}$ ($p \times 1$)
    pub beta: Mat<f64>,
    /// Standard errors of coefficients $\text{SE}(\hat{\beta})$ ($p$)
    pub se: Vec<f64>,
    /// Residuals $\hat{\varepsilon} = Y - X\hat{\beta}$ ($N \times 1$)
    pub residuals: Mat<f64>,
    /// Transformed residuals $\hat{\varepsilon}^* = L^{-1}\hat{\varepsilon}$ ($N \times 1$)
    pub transformed_residuals: Mat<f64>,
    /// Residual variance estimate $\hat{\sigma}^2_\text{ML}$
    pub sigma2_ml: f64,
    /// Unbiased residual variance estimate $\hat{\sigma}^2_\text{unbiased}$ ($N-p$ degrees of freedom)
    pub sigma2_unbiased: f64,
    /// Log-Likelihood $\ln L$ of the fitted model
    pub log_likelihood: f64,
    /// Log-determinant of the VCV matrix $\ln |V|$
    pub log_det_v: f64,
    /// Variance-covariance matrix of coefficient estimates $(X^T V^{-1} X)^{-1} \sigma^2$ ($p \times p$)
    pub vcov_beta: Mat<f64>,
}

/// Solves forward substitution $L X^* = X$ for a lower-triangular matrix $L$.
// ponytail: stdlib forward substitution for O(N^2 * p) triangular solve
pub fn forward_substitute(l: &Mat<f64>, x: &Mat<f64>) -> Mat<f64> {
    let n = x.nrows();
    let p = x.ncols();
    let mut x_star = Mat::zeros(n, p);

    for k in 0..p {
        for i in 0..n {
            let mut sum = 0.0;
            for m in 0..i {
                sum += l[(i, m)] * x_star[(m, k)];
            }
            let l_ii = l[(i, i)];
            x_star[(i, k)] = if l_ii.abs() > 1e-15 {
                (x[(i, k)] - sum) / l_ii
            } else {
                0.0
            };
        }
    }

    x_star
}

/// Computes the inverse of a small $p \times p$ matrix.
// ponytail: explicit 2x2 and faer general inverse solver
pub fn invert_small_matrix(a: &Mat<f64>) -> Result<Mat<f64>, String> {
    let p = a.nrows();
    if p != a.ncols() {
        return Err("Matrix must be square".to_string());
    }

    if p == 1 {
        let val = a[(0, 0)];
        if val.abs() < 1e-15 {
            return Err("Singular 1x1 matrix".to_string());
        }
        let mut inv = Mat::zeros(1, 1);
        inv[(0, 0)] = 1.0 / val;
        return Ok(inv);
    }

    if p == 2 {
        let det = a[(0, 0)] * a[(1, 1)] - a[(0, 1)] * a[(1, 0)];
        if det.abs() < 1e-15 {
            return Err("Singular 2x2 matrix".to_string());
        }
        let inv_det = 1.0 / det;
        let mut inv = Mat::zeros(2, 2);
        inv[(0, 0)] = a[(1, 1)] * inv_det;
        inv[(0, 1)] = -a[(0, 1)] * inv_det;
        inv[(1, 0)] = -a[(1, 0)] * inv_det;
        inv[(1, 1)] = a[(0, 0)] * inv_det;
        return Ok(inv);
    }

    // General inverse via Cholesky / PartialPivLu for p >= 3
    let lu = a.partial_piv_lu();
    Ok(lu.inverse())
}

/// Solves GLS regression system $Y = X\beta + \varepsilon$ with covariance matrix $V$.
// ponytail: complete Cholesky-based GLS solver with exact ML log-likelihood and variance metrics
pub fn solve_gls(
    vcv: &Mat<f64>,
    x: &Mat<f64>,
    y: &Mat<f64>
) -> Result<GlsResult, String> {
    let n = x.nrows();
    let p = x.ncols();

    if n != vcv.nrows() || n != vcv.ncols() || n != y.nrows() {
        return Err(format!(
            "Dimension mismatch: VCV ({}, {}), X ({}, {}), Y ({}, {})",
            vcv.nrows(),
            vcv.ncols(),
            x.nrows(),
            x.ncols(),
            y.nrows(),
            y.ncols()
        ));
    }

    // 1. Cholesky factorization V = L L^T
    let llt = vcv.as_ref()
        .llt(Side::Lower)
        .map_err(|e| format!("LLT failed: {:?}", e))?;
    let l_mat = llt.L().to_owned();

    // 2. Compute log-determinant log|V| = 2 * sum(log(L_ii))
    let mut log_det_v = 0.0;
    for i in 0..n {
        let diag = l_mat[(i, i)];
        if diag <= 0.0 {
            return Err("Non-positive diagonal in Cholesky factor".to_string());
        }
        log_det_v += 2.0 * diag.ln();
    }

    // 3. Transform X and Y via forward substitution: X* = L^-1 X, Y* = L^-1 Y
    let x_star = forward_substitute(&l_mat, x);
    let y_star = forward_substitute(&l_mat, y);

    // 4. Compute X*^T X* and X*^T Y*
    let mut xtx_star = Mat::zeros(p, p);
    for i in 0..p {
        for j in 0..p {
            let mut sum = 0.0;
            for k in 0..n {
                sum += x_star[(k, i)] * x_star[(k, j)];
            }
            xtx_star[(i, j)] = sum;
        }
    }

    let mut xty_star = Mat::zeros(p, 1);
    for i in 0..p {
        let mut sum = 0.0;
        for k in 0..n {
            sum += x_star[(k, i)] * y_star[(k, 0)];
        }
        xty_star[(i, 0)] = sum;
    }

    // 5. Invert (X*^T X*)
    let inv_xtx = invert_small_matrix(&xtx_star)?;

    // 6. Compute beta_hat = (X*^T X*)^-1 (X*^T Y*)
    let mut beta = Mat::zeros(p, 1);
    for i in 0..p {
        let mut sum = 0.0;
        for j in 0..p {
            sum += inv_xtx[(i, j)] * xty_star[(j, 0)];
        }
        beta[(i, 0)] = sum;
    }

    // 7. Compute residuals
    let mut residuals = Mat::zeros(n, 1);
    for i in 0..n {
        let mut pred = 0.0;
        for j in 0..p {
            pred += x[(i, j)] * beta[(j, 0)];
        }
        residuals[(i, 0)] = y[(i, 0)] - pred;
    }

    let transformed_residuals = forward_substitute(&l_mat, &residuals);

    // 8. Residual sum of squares RSS* = sum((e*)^2)
    let mut rss_star = 0.0;
    for i in 0..n {
        let e = transformed_residuals[(i, 0)];
        rss_star += e * e;
    }

    let sigma2_ml = rss_star / (n as f64);
    let df = if n > p { (n - p) as f64 } else { 1.0 };
    let sigma2_unbiased = rss_star / df;

    // 9. Log-Likelihood
    // ln L = -0.5 * N * ln(2*pi) - 0.5 * N * ln(sigma2_ml) - 0.5 * log_det_v - 0.5 * N
    let n_f64 = n as f64;
    let log_likelihood = -0.5 * n_f64 * (2.0 * std::f64::consts::PI).ln()
        - 0.5 * n_f64 * sigma2_ml.ln()
        - 0.5 * log_det_v
        - 0.5 * n_f64;

    // 10. Covariance and Standard Errors of beta
    let mut vcov_beta = Mat::zeros(p, p);
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
    use approx::assert_relative_eq;

    #[test]
    fn test_gls_solver_identity_vcv() {
        // Identity VCV reduces to standard OLS
        let mut vcv = Mat::zeros(2, 2);
        vcv[(0, 0)] = 1.0;
        vcv[(1, 1)] = 1.0;

        let mut x = Mat::zeros(2, 2);
        x[(0, 0)] = 1.0; x[(0, 1)] = 1.0;
        x[(1, 0)] = 1.0; x[(1, 1)] = 2.0;

        let mut y = Mat::zeros(2, 1);
        y[(0, 0)] = 2.0;
        y[(1, 0)] = 4.0;

        let res = solve_gls(&vcv, &x, &y).expect("GLS solve failed");
        // Line y = 2x -> beta = [0.0, 2.0]
        assert_relative_eq!(res.beta[(0, 0)], 0.0, epsilon = 1e-8);
        assert_relative_eq!(res.beta[(1, 0)], 2.0, epsilon = 1e-8);
    }
}
