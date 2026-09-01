//! End-to-end execution pipeline for high-throughput OrthoFinder PGLS screening.

use crate::data::comparative::{ComparativeData, ComparativeDataOptions};
use crate::genomics::engine::{run_ortho_pgls_parallel, OrthoPglsConfig};
use crate::genomics::filter::{filter_orthogroups, OrthoFilterOptions, OrthoFilterReport};
use crate::genomics::multiple_testing::{
    adjust_p_values, export_ortho_results_tsv, AnnotatedOrthoResult,
};
use crate::genomics::ortho_reader::{align_tripartite, TripartiteReport};
use crate::math::vcv::compute_vcv;
use crate::opt::optimize_pagel_lambda;

/// Configuration options for the genomic orthogroup PGLS pipeline.
#[derive(Debug, Clone)]
pub struct OrthoPipelineConfig {
    pub tree_path: String,
    pub data_path: String,
    pub ortho_path: String,
    pub response_col: String,
    pub covariate_cols: Vec<String>,
    pub output_path: String,
    pub lambda: f64,
    pub optimize_lambda: bool,
    pub min_species_count: usize,
    pub min_prevalence: Option<f64>,
    pub min_variance: f64,
    pub min_count_spread: Option<f64>,
    pub min_count_range: Option<f64>,
    pub top_ogs: Option<usize>,
    pub save_tree_path: Option<String>,
    pub eval_ranges: Option<Vec<f64>>,
    pub binary_mode: bool,
    pub log_count_mode: bool,
    pub gfs_as_response: bool,
    pub prune: bool,
    pub fuzzy_taxonomy: bool,
}

/// Summary report of the genomic PGLS pipeline execution.
#[derive(Debug, Clone)]
pub struct OrthoPipelineSummary {
    pub tri_report: TripartiteReport,
    pub filter_report: OrthoFilterReport,
    pub total_analyzed: usize,
    pub significant_fdr_05: usize,
    pub significant_fdr_01: usize,
    pub significant_bonferroni: usize,
    pub lambda_used: f64,
    pub top_hits: Vec<AnnotatedOrthoResult>,
    pub all_results: Vec<AnnotatedOrthoResult>,
}

/// Executes the full end-to-end genomic PGLS pipeline across all orthogroups.
// ponytail: clean orchestration of alignment, QC, ML lambda, parallel PGLS, FDR and TSV export
pub fn run_ortho_pipeline(config: &OrthoPipelineConfig) -> Result<OrthoPipelineSummary, String> {
    // 1. Tripartite alignment (Tree + Phenotype + OrthoFinder)
    let pheno_options = ComparativeDataOptions {
        prune: config.prune,
        fuzzy_taxonomy: config.fuzzy_taxonomy,
    };

    let tri_data = align_tripartite(
        &config.tree_path,
        &config.data_path,
        &pheno_options,
        &config.ortho_path,
    )?;

    // 2. QC Filtering
    let filter_options = OrthoFilterOptions {
        min_species_count: config.min_species_count,
        min_prevalence: config.min_prevalence,
        min_variance: config.min_variance,
        min_count_spread: config.min_count_spread,
        min_count_range: config.min_count_range,
        top_ogs: config.top_ogs,
        binary_mode: config.binary_mode,
        log_count_mode: config.log_count_mode,
        species_total_genes: None,
    };
    let filtered_ortho = filter_orthogroups(&tri_data.orthogroups, &filter_options);

    // 3. Extract Phenotype Response Y and Covariates Z
    let comp_helper = ComparativeData {
        tree: tri_data.tree.clone(),
        species_names: tri_data.species_names.clone(),
        data: tri_data.pheno_data.clone(),
        report: Default::default(),
    };

    let y = comp_helper
        .to_response_vector(&config.response_col)
        .map_err(|e| format!("Failed to extract response variable '{}': {}", config.response_col, e))?;

    let cov = if !config.covariate_cols.is_empty() {
        let cov_refs: Vec<&str> = config.covariate_cols.iter().map(|s| s.as_str()).collect();
        let z = comp_helper
            .to_design_matrix(&cov_refs)
            .map_err(|e| format!("Failed to extract covariates: {}", e))?;
        // Remove intercept column from design matrix to get pure covariates
        let n = z.nrows();
        let num_covs = config.covariate_cols.len();
        let mut pure_cov = faer::Mat::zeros(n, num_covs);
        for i in 0..n {
            for j in 0..num_covs {
                pure_cov[(i, j)] = z[(i, 1 + j)];
            }
        }
        Some(pure_cov)
    } else {
        None
    };

    // 4. Determine or Optimize Pagel's Lambda
    let lambda_used = if config.optimize_lambda {
        let base_vcv = compute_vcv(&tri_data.tree);
        let cov_refs: Vec<&str> = config.covariate_cols.iter().map(|s| s.as_str()).collect();
        let x_base = if cov_refs.is_empty() {
            let mut ones = faer::Mat::zeros(tri_data.tree.tips.len(), 1);
            for i in 0..ones.nrows() {
                ones[(i, 0)] = 1.0;
            }
            ones
        } else {
            comp_helper
                .to_design_matrix(&cov_refs)
                .map_err(|e| e.to_string())?
        };

        let opt_res = optimize_pagel_lambda(&base_vcv, &x_base, &y, 0.0, 1.0, 1e-6)?;
        opt_res.optimal_lambda
    } else {
        config.lambda
    };

    // 5. Run Parallel PGLS Engine
    let pgls_config = OrthoPglsConfig {
        lambda: lambda_used,
        model_mode: if config.gfs_as_response {
            crate::genomics::engine::OrthoModelMode::GfsResponse
        } else {
            crate::genomics::engine::OrthoModelMode::TraitResponse
        },
    };
    let raw_results = run_ortho_pgls_parallel(
        &tri_data.tree,
        &y,
        cov.as_ref(),
        &filtered_ortho,
        &pgls_config,
    )?;

    // 6. Multiple Testing FDR & Bonferroni Corrections
    let annotated_results = adjust_p_values(raw_results);

    // 7. Export TSV Results (if output path provided)
    if !config.output_path.trim().is_empty() {
        export_ortho_results_tsv(&annotated_results, &config.output_path)?;
    }

    // 8. Statistics Summary
    let total_analyzed = annotated_results.len();
    let significant_fdr_05 = annotated_results.iter().filter(|r| r.q_value_fdr < 0.05).count();
    let significant_fdr_01 = annotated_results.iter().filter(|r| r.q_value_fdr < 0.01).count();
    let significant_bonferroni = annotated_results.iter().filter(|r| r.p_bonferroni < 0.05).count();
    let top_hits = annotated_results.iter().take(10).cloned().collect();

    Ok(OrthoPipelineSummary {
        tri_report: tri_data.report,
        filter_report: filtered_ortho.filter_report,
        total_analyzed,
        significant_fdr_05,
        significant_fdr_01,
        significant_bonferroni,
        lambda_used,
        top_hits,
        all_results: annotated_results,
    })
}
