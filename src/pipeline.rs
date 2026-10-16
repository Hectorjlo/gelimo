//! High-level execution pipeline for PGLS models.

use crate::data::{ComparativeData, ComparativeDataOptions};
use crate::math::{
    apply_pagel_lambda, apply_pagel_lambda_tree, compute_vcv, solve_gls_with_mode, GlsResult,
    GlsSolverMode,
};
use crate::stats::{compute_aic, compute_r_squared};

use crate::opt::{optimize_pagel_lambda, LambdaOptResult};

/// Configuration options for running a PGLS model.
#[derive(Debug, Clone)]
pub struct PglsConfig {
    pub tree_path: String,
    pub data_path: String,
    pub response_col: String,
    pub predictor_cols: Vec<String>,
    pub lambda: f64,
    pub optimize_lambda: bool,
    pub prune: bool,
    pub solver_mode: GlsSolverMode,
}

/// Output summary containing aligned data, model fit results, and statistics.
#[derive(Debug, Clone)]
pub struct PglsSummary {
    pub comp_data: ComparativeData,
    pub gls_result: GlsResult,
    pub r_squared: f64,
    pub aic: f64,
    pub p_values: Vec<f64>,
    pub lambda_used: f64,
    pub lambda_opt: Option<LambdaOptResult>,
}

/// Executes the full PGLS pipeline from file ingestion to model fitting and statistical summary.
pub fn run_pgls_pipeline(config: &PglsConfig) -> Result<PglsSummary, String> {
    let options = ComparativeDataOptions {
        prune: config.prune,
        fuzzy_taxonomy: true,
    };

    let comp_data = ComparativeData::from_files_with_options(
        &config.tree_path,
        &config.data_path,
        &options,
    )
    .map_err(|e| format!("ComparativeData error: {}", e))?;

    let y = comp_data
        .to_response_vector(&config.response_col)
        .map_err(|e| format!("Failed to extract response vector Y: {}", e))?;

    let predictor_refs: Vec<&str> = config.predictor_cols.iter().map(|s| s.as_str()).collect();
    let x = comp_data
        .to_design_matrix(&predictor_refs)
        .map_err(|e| format!("Failed to extract design matrix X: {}", e))?;

    let base_vcv = compute_vcv(&comp_data.tree);

    let (lambda_used, gls_result, vcv_lambda, lambda_opt) = if config.optimize_lambda {
        let opt_res = optimize_pagel_lambda(&base_vcv, &x, &y, 0.0, 1.0, 1e-6)?;
        let lam = opt_res.optimal_lambda;
        let v_opt = apply_pagel_lambda(&base_vcv, lam);
        let gls_res = opt_res.gls_result.clone();
        (lam, gls_res, v_opt, Some(opt_res))
    } else {
        let lam = config.lambda;
        let scaled_tree = if (lam - 1.0).abs() > 1e-12 {
            apply_pagel_lambda_tree(&comp_data.tree, lam)
        } else {
            comp_data.tree.clone()
        };
        let v_lam = apply_pagel_lambda(&base_vcv, lam);
        let gls_res = solve_gls_with_mode(
            config.solver_mode,
            Some(&scaled_tree),
            Some(&v_lam),
            &x,
            &y,
        )?;
        (lam, gls_res, v_lam, None)
    };

    // Compute R^2 with exact projection on transformed intercept
    let l_mat = vcv_lambda
        .as_ref()
        .llt(faer::Side::Lower)
        .map_err(|e| format!("LLT failed: {:?}", e))?
        .L()
        .to_owned();
    let y_star = crate::math::solver::forward_substitute(&l_mat, &y);
    let x_star = crate::math::solver::forward_substitute(&l_mat, &x);
    let r_squared = compute_r_squared(&gls_result, &x_star, &y_star);

    let num_params = config.predictor_cols.len() + 2 + if config.optimize_lambda { 1 } else { 0 };
    let aic = compute_aic(gls_result.log_likelihood, num_params);

    // Compute two-tailed p-values for coefficients
    let mut p_values = Vec::with_capacity(gls_result.se.len());
    for (i, &se) in gls_result.se.iter().enumerate() {
        let b = gls_result.beta[(i, 0)];
        let t_stat = if se > 0.0 { (b / se).abs() } else { 0.0 };
        let p_val = 2.0 * (1.0 - crate::stats::metrics::approx_normal_cdf(t_stat));
        p_values.push(p_val.clamp(0.0, 1.0));
    }

    Ok(PglsSummary {
        comp_data,
        gls_result,
        r_squared,
        aic,
        p_values,
        lambda_used,
        lambda_opt,
    })
}
