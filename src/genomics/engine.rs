//! High-throughput parallel PGLS screening engine for tens of thousands of orthogroups using Rayon.

use faer::Mat;
use rayon::prelude::*;

use crate::genomics::filter::FilteredOrthoData;
use crate::math::pagel::apply_pagel_lambda;
use crate::math::solver::{forward_substitute, invert_small_matrix};
use crate::math::vcv::compute_vcv;
use crate::tree::PhyloTree;

/// Single orthogroup PGLS regression result.
#[derive(Debug, Clone)]
pub struct OrthoPglsResult {
    /// Orthogroup identifier (e.g. `OG0000001`).
    pub orthogroup_id: String,
    /// Estimated regression slope for the orthogroup ($\hat{\beta}_{\text{OG}}$).
    pub beta_og: f64,
    /// Standard error $\text{SE}(\hat{\beta}_{\text{OG}})$.
    pub se_og: f64,
    /// $t$-statistic.
    pub t_stat: f64,
    /// Two-tailed $p$-value.
    pub p_value: f64,
    /// Partial correlation coefficient / effect size: $r = t / \sqrt{t^2 + df}$ (Dr. Benjamin).
    pub partial_r: f64,
    /// Phylogenetic $R^2$.
    pub r_squared: f64,
    /// Unbiased residual variance $\hat{\sigma}^2$.
    pub sigma2: f64,
    /// Average gene copy number across species.
    pub mean_copies: f64,
    /// Number of species possessing at least 1 copy of the gene.
    pub species_with_gene: usize,
}

/// Model formulation mode for genomic PGLS.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OrthoModelMode {
    /// Standard GWAS mode: Trait ~ OG + Covariates (Trait is response Y).
    #[default]
    TraitResponse,
    /// Dr. Benjamin's inverse mode: GFS ~ Phenotypes (Gene family size is response Y).
    GfsResponse,
}

/// Configuration options for the high-throughput genomic PGLS engine.
#[derive(Debug, Clone)]
pub struct OrthoPglsConfig {
    /// Pagel's lambda parameter for tree branch scaling (default: 1.0 or ML estimated).
    pub lambda: f64,
    /// Model formulation mode (Trait ~ GFS vs GFS ~ Traits).
    pub model_mode: OrthoModelMode,
}

impl Default for OrthoPglsConfig {
    fn default() -> Self {
        Self {
            lambda: 1.0,
            model_mode: OrthoModelMode::TraitResponse,
        }
    }
}

/// Executes parallel PGLS screening across all filtered orthogroups in `FilteredOrthoData`.
pub fn run_ortho_pgls_parallel(
    tree: &PhyloTree,
    y: &Mat<f64>,
    covariates: Option<&Mat<f64>>,
    ortho_data: &FilteredOrthoData,
    config: &OrthoPglsConfig,
) -> Result<Vec<OrthoPglsResult>, String> {
    let n = y.nrows();
    if n != ortho_data.num_species() {
        return Err(format!(
            "Dimension mismatch: Y ({}) vs Ortho species ({})",
            n,
            ortho_data.num_species()
        ));
    }

    let base_vcv = compute_vcv(tree);
    let vcv_lambda = apply_pagel_lambda(&base_vcv, config.lambda);
    let l_mat = vcv_lambda
        .as_ref()
        .llt(faer::Side::Lower)
        .map_err(|e| format!("Cholesky factorization failed: {:?}", e))?
        .L()
        .to_owned();

    let ones = Mat::from_fn(n, 1, |_, _| 1.0);
    let x0_star = forward_substitute(&l_mat, &ones);
    let dot_x0_x0: f64 = (0..n).map(|i| x0_star[(i, 0)] * x0_star[(i, 0)]).sum();

    match config.model_mode {
        OrthoModelMode::TraitResponse => {
            run_trait_response_pgls(&l_mat, &x0_star, dot_x0_x0, y, covariates, ortho_data)
        }
        OrthoModelMode::GfsResponse => {
            run_gfs_response_pgls(&l_mat, &x0_star, dot_x0_x0, y, covariates, ortho_data)
        }
    }
}

fn compute_null_tss(x0_star: &Mat<f64>, dot_x0_x0: f64, y_star: &Mat<f64>) -> f64 {
    let n = y_star.nrows();
    let dot_x0_y: f64 = (0..n).map(|i| x0_star[(i, 0)] * y_star[(i, 0)]).sum();
    let mu_null = if dot_x0_x0 > 1e-15 { dot_x0_y / dot_x0_x0 } else { 0.0 };
    (0..n)
        .map(|i| {
            let diff = y_star[(i, 0)] - x0_star[(i, 0)] * mu_null;
            diff * diff
        })
        .sum()
}

fn extract_og_vector(og_slice: &[f64]) -> (Mat<f64>, f64, usize) {
    let n = og_slice.len();
    let mut species_with_gene = 0;
    let mut sum_copies = 0.0;
    let mut g_vec = Mat::zeros(n, 1);
    for (i, &count) in og_slice.iter().enumerate() {
        g_vec[(i, 0)] = count;
        if count > 0.0 {
            species_with_gene += 1;
        }
        sum_copies += count;
    }
    (g_vec, sum_copies / (n as f64), species_with_gene)
}

#[allow(clippy::too_many_arguments)]
fn build_ortho_result(
    orthogroup_id: String,
    beta_og: f64,
    inv_xtx_11: f64,
    rss_star: f64,
    tss_star: f64,
    df: f64,
    mean_copies: f64,
    species_with_gene: usize,
) -> OrthoPglsResult {
    let sigma2 = rss_star / df;
    let se_og = (inv_xtx_11 * sigma2).max(0.0).sqrt();
    let t_stat = if se_og > 0.0 { beta_og / se_og } else { 0.0 };
    let p_value = (2.0 * (1.0 - crate::stats::metrics::approx_normal_cdf(t_stat.abs()))).clamp(0.0, 1.0);
    let partial_r = if t_stat.is_finite() {
        t_stat / (t_stat * t_stat + df).sqrt()
    } else {
        0.0
    };
    let r_squared = if tss_star > 1e-15 {
        (1.0 - rss_star / tss_star).clamp(0.0, 1.0)
    } else {
        0.0
    };

    OrthoPglsResult {
        orthogroup_id,
        beta_og,
        se_og,
        t_stat,
        p_value,
        partial_r,
        r_squared,
        sigma2,
        mean_copies,
        species_with_gene,
    }
}

/// Mode A: Trait ~ GFS + Covariates (Trait as Response)
fn run_trait_response_pgls(
    l_mat: &Mat<f64>,
    x0_star: &Mat<f64>,
    dot_x0_x0: f64,
    y: &Mat<f64>,
    covariates: Option<&Mat<f64>>,
    ortho_data: &FilteredOrthoData,
) -> Result<Vec<OrthoPglsResult>, String> {
    let n = y.nrows();
    let y_star = forward_substitute(l_mat, y);
    let cov_star = covariates.map(|cov| forward_substitute(l_mat, cov));
    let num_covs = covariates.map_or(0, |c| c.ncols());

    let p = 2 + num_covs;
    let df = if n > p { (n - p) as f64 } else { 1.0 };
    let tss_star = compute_null_tss(x0_star, dot_x0_x0, &y_star);
    let num_ogs = ortho_data.num_orthogroups();

    let results: Vec<OrthoPglsResult> = (0..num_ogs)
        .into_par_iter()
        .map(|og_idx| {
            let og_id = &ortho_data.orthogroup_ids[og_idx];
            let (g_vec, mean_copies, species_with_gene) =
                extract_og_vector(ortho_data.get_og_counts(og_idx));
            let g_star = forward_substitute(l_mat, &g_vec);

            let mut x_star = Mat::zeros(n, p);
            for i in 0..n {
                x_star[(i, 0)] = x0_star[(i, 0)];
                x_star[(i, 1)] = g_star[(i, 0)];
                if let Some(ref c_star) = cov_star {
                    for j in 0..num_covs {
                        x_star[(i, 2 + j)] = c_star[(i, j)];
                    }
                }
            }

            let xtx = x_star.transpose() * &x_star;
            let xty = x_star.transpose() * &y_star;

            let Ok(inv_xtx) = invert_small_matrix(&xtx) else {
                return OrthoPglsResult {
                    orthogroup_id: og_id.clone(),
                    beta_og: 0.0,
                    se_og: f64::NAN,
                    t_stat: 0.0,
                    p_value: 1.0,
                    partial_r: 0.0,
                    r_squared: 0.0,
                    sigma2: f64::NAN,
                    mean_copies,
                    species_with_gene,
                };
            };

            let beta = &inv_xtx * &xty;
            let err = &y_star - &x_star * &beta;
            let rss_star: f64 = (0..n).map(|i| err[(i, 0)] * err[(i, 0)]).sum();

            build_ortho_result(
                og_id.clone(),
                beta[(1, 0)],
                inv_xtx[(1, 1)],
                rss_star,
                tss_star,
                df,
                mean_copies,
                species_with_gene,
            )
        })
        .collect();

    Ok(results)
}

/// Mode B: Dr. Benjamin's Inverse Mode (GFS ~ Phenotypes)
fn run_gfs_response_pgls(
    l_mat: &Mat<f64>,
    x0_star: &Mat<f64>,
    dot_x0_x0: f64,
    pheno_y: &Mat<f64>,
    covariates: Option<&Mat<f64>>,
    ortho_data: &FilteredOrthoData,
) -> Result<Vec<OrthoPglsResult>, String> {
    let n = pheno_y.nrows();
    let num_covs = covariates.map_or(0, |c| c.ncols());
    let p = 2 + num_covs;
    let df = if n > p { (n - p) as f64 } else { 1.0 };

    let mut x_raw = Mat::zeros(n, p);
    for i in 0..n {
        x_raw[(i, 0)] = 1.0;
        x_raw[(i, 1)] = pheno_y[(i, 0)];
        if let Some(cov) = covariates {
            for j in 0..num_covs {
                x_raw[(i, 2 + j)] = cov[(i, j)];
            }
        }
    }

    let x_star = forward_substitute(l_mat, &x_raw);
    let xtx = x_star.transpose() * &x_star;
    let inv_xtx = invert_small_matrix(&xtx)
        .map_err(|e| format!("Singular phenotype design matrix: {}", e))?;

    // Precompute projection matrix P = inv(X*^T X*) * X*^T of size (p x N) ONCE
    let p_mat = &inv_xtx * x_star.transpose();
    let num_ogs = ortho_data.num_orthogroups();

    let results: Vec<OrthoPglsResult> = (0..num_ogs)
        .into_par_iter()
        .map(|og_idx| {
            let og_id = &ortho_data.orthogroup_ids[og_idx];
            let (gfs_raw, mean_copies, species_with_gene) =
                extract_og_vector(ortho_data.get_og_counts(og_idx));
            let gfs_star = forward_substitute(l_mat, &gfs_raw);

            let beta = &p_mat * &gfs_star;
            let err = &gfs_star - &x_star * &beta;
            let rss_star: f64 = (0..n).map(|i| err[(i, 0)] * err[(i, 0)]).sum();
            let tss_star = compute_null_tss(x0_star, dot_x0_x0, &gfs_star);

            build_ortho_result(
                og_id.clone(),
                beta[(1, 0)],
                inv_xtx[(1, 1)],
                rss_star,
                tss_star,
                df,
                mean_copies,
                species_with_gene,
            )
        })
        .collect();

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::comparative::ComparativeData;
    use crate::genomics::filter::{filter_orthogroups, OrthoFilterOptions};
    use crate::genomics::ortho_reader::OrthoGeneCountTable;

    #[test]
    fn test_ortho_pgls_parallel_engine() {
        let comp_data = ComparativeData::from_files("data/benchmark_3sp.tree", "data/benchmark_3sp.tsv").unwrap();
        let y = comp_data.to_response_vector("Maximum longevity (yrs)").unwrap();

        let table = OrthoGeneCountTable {
            orthogroup_ids: vec!["OG_test1".into(), "OG_test2".into()],
            species_names: comp_data.species_names.clone(),
            counts: vec![
                1.0, 2.0, 4.0,
                3.0, 1.0, 0.0,
            ],
        };

        let filter_opts = OrthoFilterOptions {
            min_species_count: 1,
            min_variance: 1e-6,
            ..Default::default()
        };
        let filtered = filter_orthogroups(&table, &filter_opts);
        let config = OrthoPglsConfig::default();

        let results = run_ortho_pgls_parallel(&comp_data.tree, &y, None, &filtered, &config).unwrap();

        assert_eq!(results.len(), 2);
        assert_eq!(results[0].orthogroup_id, "OG_test1");
        assert!(results[0].beta_og.is_finite());
        assert!(results[0].partial_r.is_finite());
        assert!(results[0].p_value >= 0.0 && results[0].p_value <= 1.0);
    }
}
