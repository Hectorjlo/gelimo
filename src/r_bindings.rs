//! R bindings for `gelimo` using `extendr-api`.

use extendr_api::prelude::*;
use crate::genomics::{run_ortho_pipeline, OrthoPipelineConfig};
use crate::math::solver::GlsSolverMode;
use crate::pipeline::{run_pgls_pipeline, PglsConfig};

/// Runs high-throughput parallel PGLS on OrthoFinder orthogroups from R.
///
/// Returns an R data.frame with coefficients, p-values, FDR q-values, and Partial r effect sizes.
#[extendr]
#[allow(clippy::too_many_arguments)]
fn r_pgls_orthogroups(
    tree: String,
    data: String,
    orthogroups: String,
    response: String,
    covariates: Nullable<Vec<String>>,
    output_tsv: Nullable<String>,
    log_count: Nullable<bool>,
    gfs_as_response: Nullable<bool>,
    lambda: Nullable<f64>,
    opt_lambda: Nullable<bool>,
    min_species: Nullable<i32>,
    min_prevalence: Nullable<f64>,
    min_variance: Nullable<f64>,
    min_count_spread: Nullable<f64>,
    min_count_range: Nullable<f64>,
    top_ogs: Nullable<i32>,
    save_tree: Nullable<String>,
    prune: Nullable<bool>,
) -> extendr_api::Result<Robj> {
    let covariate_cols = covariates.into_option().unwrap_or_default();
    let output_path = output_tsv.into_option().unwrap_or_default();
    let log_count_mode = log_count.into_option().unwrap_or(false);
    let gfs_resp = gfs_as_response.into_option().unwrap_or(false);
    let lam = lambda.into_option().unwrap_or(1.0);
    let optimize_lambda = opt_lambda.into_option().unwrap_or(false);
    let min_species_count = min_species.into_option().unwrap_or(5).max(0) as usize;
    let min_prev = min_prevalence.into_option();
    let min_var = min_variance.into_option().unwrap_or(1e-6);
    let min_spread = min_count_spread.into_option();
    let min_range = min_count_range.into_option();
    let top_n = top_ogs.into_option().map(|n| n.max(1) as usize);
    let save_tree_path = save_tree.into_option();
    let prune_data = prune.into_option().unwrap_or(true);

    let config = OrthoPipelineConfig {
        tree_path: tree,
        data_path: data,
        ortho_path: orthogroups,
        response_col: response,
        covariate_cols,
        output_path,
        lambda: lam,
        optimize_lambda,
        min_species_count,
        min_prevalence: min_prev,
        min_variance: min_var,
        min_count_spread: min_spread,
        min_count_range: min_range,
        top_ogs: top_n,
        save_tree_path,
        eval_ranges: None,
        binary_mode: false,
        log_count_mode,
        gfs_as_response: gfs_resp,
        prune: prune_data,
        fuzzy_taxonomy: true,
    };

    let summary = run_ortho_pipeline(&config).map_err(extendr_api::Error::Other)?;

    let mut og_ids = Vec::with_capacity(summary.all_results.len());
    let mut betas = Vec::with_capacity(summary.all_results.len());
    let mut ses = Vec::with_capacity(summary.all_results.len());
    let mut t_stats = Vec::with_capacity(summary.all_results.len());
    let mut p_vals = Vec::with_capacity(summary.all_results.len());
    let mut q_vals = Vec::with_capacity(summary.all_results.len());
    let mut p_bonfs = Vec::with_capacity(summary.all_results.len());
    let mut partial_rs = Vec::with_capacity(summary.all_results.len());
    let mut r2s = Vec::with_capacity(summary.all_results.len());
    let mut sigmas = Vec::with_capacity(summary.all_results.len());
    let mut mean_copies = Vec::with_capacity(summary.all_results.len());
    let mut species_counts = Vec::with_capacity(summary.all_results.len());

    for hit in &summary.all_results {
        og_ids.push(hit.orthogroup_id.clone());
        betas.push(hit.beta_og);
        ses.push(hit.se_og);
        t_stats.push(hit.t_stat);
        p_vals.push(hit.p_value);
        q_vals.push(hit.q_value_fdr);
        p_bonfs.push(hit.p_bonferroni);
        partial_rs.push(hit.partial_r);
        r2s.push(hit.r_squared);
        sigmas.push(hit.sigma2);
        mean_copies.push(hit.mean_copies);
        species_counts.push(hit.species_with_gene as i32);
    }

    let df = data_frame!(
        Orthogroup = og_ids,
        Beta_OG = betas,
        StdError = ses,
        t_stat = t_stats,
        p_value = p_vals,
        q_value_fdr = q_vals,
        p_bonferroni = p_bonfs,
        Partial_r = partial_rs,
        R_squared = r2s,
        Sigma2 = sigmas,
        Mean_Copies = mean_copies,
        Species_With_Gene = species_counts
    );

    Ok(df)
}

/// Fits a single PGLS model on comparative species traits from R.
///
/// Returns an R list containing coefficients, R-squared, AIC, log-likelihood, and Pagel's lambda.
#[extendr]
#[allow(clippy::too_many_arguments)]
fn r_pgls_fit(
    tree: String,
    data: String,
    response: String,
    predictors: Nullable<Vec<String>>,
    lambda: Nullable<f64>,
    opt_lambda: Nullable<bool>,
    prune: Nullable<bool>,
    solver: Nullable<String>,
) -> extendr_api::Result<Robj> {
    let predictor_cols = predictors.into_option().unwrap_or_default();
    let lam = lambda.into_option().unwrap_or(1.0);
    let optimize_lambda = opt_lambda.into_option().unwrap_or(false);
    let prune_data = prune.into_option().unwrap_or(true);
    let solver_str = solver.into_option().unwrap_or_else(|| "dense".into());
    let solver_mode = match solver_str.to_lowercase().as_str() {
        "linear" | "ho_ane" => GlsSolverMode::LinearHoAne,
        _ => GlsSolverMode::DenseCholesky,
    };

    let config = PglsConfig {
        tree_path: tree,
        data_path: data,
        response_col: response,
        predictor_cols: predictor_cols.clone(),
        lambda: lam,
        optimize_lambda,
        prune: prune_data,
        solver_mode,
    };

    let summary = run_pgls_pipeline(&config).map_err(extendr_api::Error::Other)?;

    let mut param_names = vec!["(Intercept)".to_string()];
    param_names.extend(predictor_cols);

    let k = summary.gls_result.beta.nrows();
    let mut names = Vec::with_capacity(k);
    let mut estimates = Vec::with_capacity(k);
    let mut std_errors = Vec::with_capacity(k);
    let mut t_values = Vec::with_capacity(k);
    let mut p_values = Vec::with_capacity(k);

    for (i, p_name) in param_names.into_iter().enumerate() {
        if i < k {
            let b = summary.gls_result.beta[(i, 0)];
            let se = summary.gls_result.se[i];
            let t_val = if se > 0.0 { b / se } else { 0.0 };
            let p_val = summary.p_values[i];

            names.push(p_name);
            estimates.push(b);
            std_errors.push(se);
            t_values.push(t_val);
            p_values.push(p_val);
        }
    }

    let coef_df = data_frame!(
        Parameter = names,
        Estimate = estimates,
        StdError = std_errors,
        t_value = t_values,
        p_value = p_values
    );

    let opt_lambda_val = summary
        .lambda_opt
        .as_ref()
        .map(|o| o.optimal_lambda)
        .unwrap_or(summary.lambda_used);

    let res_list = list!(
        coefficients = coef_df,
        r_squared = summary.r_squared,
        aic = summary.aic,
        log_likelihood = summary.gls_result.log_likelihood,
        sigma2 = summary.gls_result.sigma2_ml,
        lambda = summary.lambda_used,
        optimal_lambda = opt_lambda_val,
        matched_species = summary.comp_data.species_names.len() as i32
    );

    Ok(res_list.into())
}

extendr_module! {
    mod gelimo;
    fn r_pgls_orthogroups;
    fn r_pgls_fit;
}
