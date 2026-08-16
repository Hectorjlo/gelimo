//! Maximum Likelihood estimation of Pagel's lambda using bounded 1D optimization and Rayon parallelism.

use crate::math::pagel::apply_pagel_lambda;
use crate::math::solver::{solve_gls, GlsResult};
use faer::Mat;
use rayon::prelude::*;

/// Output summary for Pagel's $\lambda$ Maximum Likelihood optimization.
#[derive(Debug, Clone)]
pub struct LambdaOptResult {
    /// Maximum Likelihood estimate $\hat{\lambda}_{\text{ML}}$.
    pub optimal_lambda: f64,
    /// Log-Likelihood at $\hat{\lambda}_{\text{ML}}$.
    pub log_likelihood_ml: f64,
    /// Log-Likelihood under $\lambda = 0.0$ (no phylogenetic signal / OLS).
    pub log_likelihood_lambda0: f64,
    /// Log-Likelihood under $\lambda = 1.0$ (strict Brownian Motion).
    pub log_likelihood_lambda1: f64,
    /// Likelihood Ratio Test $p$-value comparing $\hat{\lambda}_{\text{ML}}$ vs $\lambda = 0.0$.
    pub lrt_p_value_vs_0: f64,
    /// Likelihood Ratio Test $p$-value comparing $\hat{\lambda}_{\text{ML}}$ vs $\lambda = 1.0$.
    pub lrt_p_value_vs_1: f64,
    /// Model solution at optimal $\lambda$.
    pub gls_result: GlsResult,
}

/// Evaluates Log-Likelihood across a grid of $\lambda$ values in parallel using Rayon (zero-copy).
// ponytail: zero-copy data-parallel grid search using rayon
pub fn grid_search_lambda_parallel(
    vcv: &Mat<f64>,
    x: &Mat<f64>,
    y: &Mat<f64>,
    grid_points: usize,
    lower: f64,
    upper: f64,
) -> Vec<(f64, f64)> {
    let steps = grid_points.max(2);
    let step_size = (upper - lower) / ((steps - 1) as f64);

    let lambdas: Vec<f64> = (0..steps)
        .map(|i| lower + (i as f64) * step_size)
        .collect();

    lambdas
        .into_par_iter()
        .map(|lam| {
            let v_lam = apply_pagel_lambda(vcv, lam);
            let log_lik = match solve_gls(&v_lam, x, y) {
                Ok(res) => res.log_likelihood,
                Err(_) => f64::NEG_INFINITY,
            };
            (lam, log_lik)
        })
        .collect()
}

/// Optimizes Pagel's $\lambda$ by Maximum Likelihood within `[lower, upper]` bounds using Golden-Section search.
// ponytail: zero-allocation 1D bounded optimization for Pagel lambda
pub fn optimize_pagel_lambda(
    vcv: &Mat<f64>,
    x: &Mat<f64>,
    y: &Mat<f64>,
    lower: f64,
    upper: f64,
    tolerance: f64,
) -> Result<LambdaOptResult, String> {
    let tol = tolerance.max(1e-8);
    let phi = (5.0_f64.sqrt() - 1.0) / 2.0; // ~0.6180339887

    let mut a = lower.max(0.0);
    let mut b = upper.max(a + 1e-4);

    let mut c = b - phi * (b - a);
    let mut d = a + phi * (b - a);

    let eval_lik = |lam: f64| -> f64 {
        let v_lam = apply_pagel_lambda(vcv, lam);
        match solve_gls(&v_lam, x, y) {
            Ok(res) => res.log_likelihood,
            Err(_) => f64::NEG_INFINITY,
        }
    };

    let mut fc = eval_lik(c);
    let mut fd = eval_lik(d);

    let max_iter = 100;
    for _ in 0..max_iter {
        if (b - a).abs() < tol {
            break;
        }

        if fc > fd {
            b = d;
            d = c;
            fd = fc;
            c = b - phi * (b - a);
            fc = eval_lik(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + phi * (b - a);
            fd = eval_lik(d);
        }
    }

    let optimal_lambda = (a + b) / 2.0;
    let v_opt = apply_pagel_lambda(vcv, optimal_lambda);
    let gls_result = solve_gls(&v_opt, x, y)?;
    let log_likelihood_ml = gls_result.log_likelihood;

    // Evaluate baseline models lambda = 0.0 and lambda = 1.0 for LRT tests
    let v_0 = apply_pagel_lambda(vcv, 0.0);
    let log_likelihood_lambda0 = solve_gls(&v_0, x, y)?.log_likelihood;

    let v_1 = apply_pagel_lambda(vcv, 1.0);
    let log_likelihood_lambda1 = solve_gls(&v_1, x, y)?.log_likelihood;

    // Likelihood ratio tests: LRT = 2 * (logL_ML - logL_null) ~ ChiSquare(df=1)
    let lrt_0 = (2.0 * (log_likelihood_ml - log_likelihood_lambda0)).max(0.0);
    let lrt_p_value_vs_0 = compute_chisq1_p_value(lrt_0);

    let lrt_1 = (2.0 * (log_likelihood_ml - log_likelihood_lambda1)).max(0.0);
    let lrt_p_value_vs_1 = compute_chisq1_p_value(lrt_1);

    Ok(LambdaOptResult {
        optimal_lambda,
        log_likelihood_ml,
        log_likelihood_lambda0,
        log_likelihood_lambda1,
        lrt_p_value_vs_0,
        lrt_p_value_vs_1,
        gls_result,
    })
}

/// Computes upper-tail p-value for 1 degree of freedom Chi-Square distribution.
fn compute_chisq1_p_value(stat: f64) -> f64 {
    if stat <= 0.0 {
        return 1.0;
    }
    // For df = 1: P(ChiSq >= x) = 2 * (1 - Phi(sqrt(x)))
    let z = stat.sqrt();
    let p = 2.0 * (1.0 - crate::stats::metrics::approx_normal_cdf(z));
    p.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::ComparativeData;
    use crate::math::compute_vcv;
    use approx::assert_relative_eq;

    #[test]
    fn test_optimize_pagel_lambda_3sp() {
        let comp_data = ComparativeData::from_files("data/benchmark_3sp.tree", "data/benchmark_3sp.tsv").unwrap();
        let y = comp_data.to_response_vector("Maximum longevity (yrs)").unwrap();
        let x = comp_data.to_design_matrix(&["Adult weight (g)"]).unwrap();
        let vcv = compute_vcv(&comp_data.tree);

        let opt_res = optimize_pagel_lambda(&vcv, &x, &y, 0.0, 1.0, 1e-6).unwrap();

        assert!(opt_res.optimal_lambda >= 0.0 && opt_res.optimal_lambda <= 1.0);
        assert!(opt_res.log_likelihood_ml >= opt_res.log_likelihood_lambda0 - 1e-4);
        assert!(opt_res.log_likelihood_ml >= opt_res.log_likelihood_lambda1 - 1e-4);
    }

    #[test]
    fn test_grid_search_lambda_parallel_3sp() {
        let comp_data = ComparativeData::from_files("data/benchmark_3sp.tree", "data/benchmark_3sp.tsv").unwrap();
        let y = comp_data.to_response_vector("Maximum longevity (yrs)").unwrap();
        let x = comp_data.to_design_matrix(&["Adult weight (g)"]).unwrap();
        let vcv = compute_vcv(&comp_data.tree);

        let grid = grid_search_lambda_parallel(&vcv, &x, &y, 11, 0.0, 1.0);

        assert_eq!(grid.len(), 11);
        assert_relative_eq!(grid[0].0, 0.0, epsilon = 1e-12);
        assert_relative_eq!(grid[10].0, 1.0, epsilon = 1e-12);
        for (_, loglik) in &grid {
            assert!(loglik.is_finite());
        }
    }
}
