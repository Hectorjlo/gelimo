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
// ponytail: one-time Cholesky factorization + Rayon parallel solve sweep
pub fn run_ortho_pgls_parallel(
    tree: &PhyloTree,
    y: &Mat<f64>,
    covariates: Option<&Mat<f64>>,
    ortho_data: &FilteredOrthoData,
    config: &OrthoPglsConfig,
) -> Result<Vec<OrthoPglsResult>, String> {
    match config.model_mode {
        OrthoModelMode::TraitResponse => {
            run_trait_response_pgls(tree, y, covariates, ortho_data, config)
        }
        OrthoModelMode::GfsResponse => {
            run_gfs_response_pgls(tree, y, covariates, ortho_data, config)
        }
    }
}

/// Mode A: Trait ~ GFS + Covariates (Trait as Response)
fn run_trait_response_pgls(
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

    let llt = vcv_lambda
        .as_ref()
        .llt(faer::Side::Lower)
        .map_err(|e| format!("Cholesky factorization failed: {:?}", e))?;
    let l_mat = llt.L().to_owned();

    let y_star = forward_substitute(&l_mat, y);

    let mut ones = Mat::zeros(n, 1);
    for i in 0..n {
        ones[(i, 0)] = 1.0;
    }
    let x0_star = forward_substitute(&l_mat, &ones);

    let cov_star = covariates.map(|cov| forward_substitute(&l_mat, cov));
    let num_covs = covariates.map_or(0, |c| c.ncols());

    let p = 2 + num_covs;
    let df = if n > p { (n - p) as f64 } else { 1.0 };

    let mut dot_x0_y = 0.0;
    let mut dot_x0_x0 = 0.0;
    for i in 0..n {
        let x0 = x0_star[(i, 0)];
        let y_val = y_star[(i, 0)];
        dot_x0_y += x0 * y_val;
        dot_x0_x0 += x0 * x0;
    }
    let mu_null = if dot_x0_x0 > 1e-15 {
        dot_x0_y / dot_x0_x0
    } else {
        0.0
    };

    let mut tss_star = 0.0;
    for i in 0..n {
        let diff = y_star[(i, 0)] - x0_star[(i, 0)] * mu_null;
        tss_star += diff * diff;
    }

    let num_ogs = ortho_data.num_orthogroups();

    let results: Vec<OrthoPglsResult> = (0..num_ogs)
        .into_par_iter()
        .map(|og_idx| {
            let og_slice = ortho_data.get_og_counts(og_idx);
            let og_id = &ortho_data.orthogroup_ids[og_idx];

            let mut species_with_gene = 0;
            let mut sum_copies = 0.0;
            let mut g_vec = Mat::zeros(n, 1);
            for i in 0..n {
                let count = og_slice[i];
                g_vec[(i, 0)] = count;
                if count > 0.0 {
                    species_with_gene += 1;
                }
                sum_copies += count;
            }
            let mean_copies = sum_copies / (n as f64);

            let g_star = forward_substitute(&l_mat, &g_vec);

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

            let mut xtx = Mat::zeros(p, p);
            let mut xty = Mat::zeros(p, 1);
            for i in 0..p {
                for j in 0..p {
                    let mut sum = 0.0;
                    for k in 0..n {
                        sum += x_star[(k, i)] * x_star[(k, j)];
                    }
                    xtx[(i, j)] = sum;
                }
                let mut sum_y = 0.0;
                for k in 0..n {
                    sum_y += x_star[(k, i)] * y_star[(k, 0)];
                }
                xty[(i, 0)] = sum_y;
            }

            let inv_xtx = match invert_small_matrix(&xtx) {
                Ok(inv) => inv,
                Err(_) => {
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
                }
            };

            let mut beta = Mat::zeros(p, 1);
            for i in 0..p {
                let mut sum = 0.0;
                for j in 0..p {
                    sum += inv_xtx[(i, j)] * xty[(j, 0)];
                }
                beta[(i, 0)] = sum;
            }

            let mut rss_star = 0.0;
            for i in 0..n {
                let mut pred = 0.0;
                for j in 0..p {
                    pred += x_star[(i, j)] * beta[(j, 0)];
                }
                let e = y_star[(i, 0)] - pred;
                rss_star += e * e;
            }

            let sigma2 = rss_star / df;
            let se_og = (inv_xtx[(1, 1)] * sigma2).max(0.0).sqrt();
            let beta_og = beta[(1, 0)];

            let t_stat = if se_og > 0.0 { beta_og / se_og } else { 0.0 };
            let p_value = 2.0 * (1.0 - approx_normal_cdf(t_stat.abs())).clamp(0.0, 1.0);
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
                orthogroup_id: og_id.clone(),
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
        })
        .collect();

    Ok(results)
}

/// Mode B: Dr. Benjamin's Inverse Mode (GFS ~ Phenotypes)
// ponytail: precomputed projection matrix P = (X*^T X*)^-1 X*^T for 10ms genome sweep
fn run_gfs_response_pgls(
    tree: &PhyloTree,
    pheno_y: &Mat<f64>,
    covariates: Option<&Mat<f64>>,
    ortho_data: &FilteredOrthoData,
    config: &OrthoPglsConfig,
) -> Result<Vec<OrthoPglsResult>, String> {
    let n = pheno_y.nrows();
    if n != ortho_data.num_species() {
        return Err(format!(
            "Dimension mismatch: Phenotype ({}) vs Ortho species ({})",
            n,
            ortho_data.num_species()
        ));
    }

    let base_vcv = compute_vcv(tree);
    let vcv_lambda = apply_pagel_lambda(&base_vcv, config.lambda);

    let llt = vcv_lambda
        .as_ref()
        .llt(faer::Side::Lower)
        .map_err(|e| format!("Cholesky factorization failed: {:?}", e))?;
    let l_mat = llt.L().to_owned();

    let num_covs = covariates.map_or(0, |c| c.ncols());
    let p = 2 + num_covs; // Intercept + Main Phenotype + Covariates
    let df = if n > p { (n - p) as f64 } else { 1.0 };

    // Build design matrix X: [1, Phenotype_Main, Covariates...] (N x p)
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

    let x_star = forward_substitute(&l_mat, &x_raw);

    // Compute (X*^T X*) and its inverse (p x p) ONCE for all orthogroups!
    let mut xtx = Mat::zeros(p, p);
    for i in 0..p {
        for j in 0..p {
            let mut sum = 0.0;
            for k in 0..n {
                sum += x_star[(k, i)] * x_star[(k, j)];
            }
            xtx[(i, j)] = sum;
        }
    }

    let inv_xtx = invert_small_matrix(&xtx)
        .map_err(|e| format!("Singular phenotype design matrix: {}", e))?;

    // Precompute projection matrix P = inv(X*^T X*) * X*^T of size (p x N) ONCE!
    let mut p_mat = Mat::zeros(p, n);
    for i in 0..p {
        for k in 0..n {
            let mut sum = 0.0;
            for j in 0..p {
                sum += inv_xtx[(i, j)] * x_star[(k, j)];
            }
            p_mat[(i, k)] = sum;
        }
    }

    // Precompute x0_star for null model TSS calculation
    let mut x0_raw = Mat::zeros(n, 1);
    for i in 0..n {
        x0_raw[(i, 0)] = 1.0;
    }
    let x0_star = forward_substitute(&l_mat, &x0_raw);
    let dot_x0_x0 = (0..n).map(|i| x0_star[(i, 0)] * x0_star[(i, 0)]).sum::<f64>();

    let num_ogs = ortho_data.num_orthogroups();

    // Fast parallel sweep across M orthogroups
    let results: Vec<OrthoPglsResult> = (0..num_ogs)
        .into_par_iter()
        .map(|og_idx| {
            let og_slice = ortho_data.get_og_counts(og_idx);
            let og_id = &ortho_data.orthogroup_ids[og_idx];

            let mut species_with_gene = 0;
            let mut sum_copies = 0.0;
            let mut gfs_raw = Mat::zeros(n, 1);
            for i in 0..n {
                let count = og_slice[i];
                gfs_raw[(i, 0)] = count;
                if count > 0.0 {
                    species_with_gene += 1;
                }
                sum_copies += count;
            }
            let mean_copies = sum_copies / (n as f64);

            // Forward substitute GFS response vector Y_k* = L^-1 Y_k
            let gfs_star = forward_substitute(&l_mat, &gfs_raw);

            // Compute beta = P * gfs_star (size p x 1)
            let mut beta = Mat::zeros(p, 1);
            for i in 0..p {
                let mut sum = 0.0;
                for k in 0..n {
                    sum += p_mat[(i, k)] * gfs_star[(k, 0)];
                }
                beta[(i, 0)] = sum;
            }

            // Compute RSS* = ||gfs_star - X_star * beta||^2
            let mut rss_star = 0.0;
            for i in 0..n {
                let mut pred = 0.0;
                for j in 0..p {
                    pred += x_star[(i, j)] * beta[(j, 0)];
                }
                let e = gfs_star[(i, 0)] - pred;
                rss_star += e * e;
            }

            // Compute null model TSS*
            let dot_x0_gfs = (0..n).map(|i| x0_star[(i, 0)] * gfs_star[(i, 0)]).sum::<f64>();
            let mu_null = if dot_x0_x0 > 1e-15 { dot_x0_gfs / dot_x0_x0 } else { 0.0 };
            let mut tss_star = 0.0;
            for i in 0..n {
                let diff = gfs_star[(i, 0)] - x0_star[(i, 0)] * mu_null;
                tss_star += diff * diff;
            }

            let sigma2 = rss_star / df;
            // Coef 1 is the main phenotype of interest
            let se_og = (inv_xtx[(1, 1)] * sigma2).max(0.0).sqrt();
            let beta_og = beta[(1, 0)];

            let t_stat = if se_og > 0.0 { beta_og / se_og } else { 0.0 };
            let p_value = 2.0 * (1.0 - approx_normal_cdf(t_stat.abs())).clamp(0.0, 1.0);
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
                orthogroup_id: og_id.clone(),
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
        })
        .collect();

    Ok(results)
}

fn approx_normal_cdf(x: f64) -> f64 {
    0.5 * (1.0 + erff(x / std::f64::consts::SQRT_2))
}

fn erff(x: f64) -> f64 {
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;
    let p = 0.3275911;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x_abs = x.abs();
    let t = 1.0 / (1.0 + p * x_abs);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x_abs * x_abs).exp();

    sign * y
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
            binary_mode: false,
            log_count_mode: false,
            species_total_genes: None,
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
