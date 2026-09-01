//! Quality Control (QC) filtering and transformations for OrthoFinder orthogroups.

use crate::genomics::ortho_reader::OrthoGeneCountTable;

/// Configuration options for filtering and transforming orthogroup gene counts.
#[derive(Debug, Clone)]
pub struct OrthoFilterOptions {
    /// Minimum number of species that must possess at least 1 gene copy.
    pub min_species_count: usize,
    /// Minimum prevalence as a proportion (0.0 to 1.0) or percentage (0 to 100).
    pub min_prevalence: Option<f64>,
    /// Minimum variance of gene counts across species.
    pub min_variance: f64,
    /// Minimum interdecile count spread (P95 - P05) across species.
    pub min_count_spread: Option<f64>,
    /// Minimum absolute count range (max - min) across species.
    pub min_count_range: Option<f64>,
    /// Maximum number of top most variable orthogroups to retain.
    pub top_ogs: Option<usize>,
    /// Convert gene counts into binary presence/absence (0 vs 1).
    pub binary_mode: bool,
    /// Apply log10(count + 1) transformation to gene counts.
    pub log_count_mode: bool,
    /// Optional total gene counts per species for genome-size normalization (per 10k genes).
    pub species_total_genes: Option<Vec<f64>>,
}

impl Default for OrthoFilterOptions {
    fn default() -> Self {
        Self {
            min_species_count: 0,
            min_prevalence: None,
            min_variance: 1e-7,
            min_count_spread: None,
            min_count_range: None,
            top_ogs: None,
            binary_mode: false,
            log_count_mode: false,
            species_total_genes: None,
        }
    }
}

/// QC filtering summary report.
#[derive(Debug, Clone, Default)]
pub struct OrthoFilterReport {
    pub initial_orthogroups: usize,
    pub retained_orthogroups: usize,
    pub discarded_invariant: usize,
    pub discarded_low_prevalence: usize,
    pub discarded_low_variance: usize,
    pub discarded_low_range: usize,
    pub discarded_low_spread: usize,
    pub discarded_top_ogs: usize,
}

/// Filtered and transformed orthogroup matrix ready for high-throughput PGLS.
#[derive(Debug, Clone)]
pub struct FilteredOrthoData {
    pub orthogroup_ids: Vec<String>,
    pub original_indices: Vec<usize>,
    pub species_names: Vec<String>,
    /// Filtered flat matrix of shape (M_retained x N_species).
    pub counts: Vec<f64>,
    pub filter_report: OrthoFilterReport,
}

impl FilteredOrthoData {
    /// Number of retained orthogroups (rows).
    pub fn num_orthogroups(&self) -> usize {
        self.orthogroup_ids.len()
    }

    /// Number of species (columns).
    pub fn num_species(&self) -> usize {
        self.species_names.len()
    }

    /// Returns a slice of gene counts across species for a retained orthogroup index.
    // ponytail: zero-copy row slice access
    pub fn get_og_counts(&self, og_idx: usize) -> &[f64] {
        let n = self.species_names.len();
        let start = og_idx * n;
        &self.counts[start..start + n]
    }
}

/// Applies QC filtering and optional transformations to an `OrthoGeneCountTable`.
// ponytail: single-pass streaming QC filter with zero allocation per discarded OG
pub fn filter_orthogroups(
    raw_table: &OrthoGeneCountTable,
    options: &OrthoFilterOptions,
) -> FilteredOrthoData {
    let num_ogs = raw_table.num_orthogroups();
    let num_species = raw_table.num_species();
    let n_f64 = num_species as f64;

    // Calculate effective minimum species threshold
    let effective_min_species = match (options.min_prevalence, options.min_species_count) {
        (Some(prev), count) => {
            let prop = if prev > 1.0 { prev / 100.0 } else { prev };
            let from_prev = (prop * n_f64).ceil() as usize;
            count.max(from_prev).min(num_species)
        }
        (None, count) => count.min(num_species),
    };

    let mut candidate_ids = Vec::with_capacity(num_ogs / 2);
    let mut candidate_indices = Vec::with_capacity(num_ogs / 2);
    let mut candidate_counts = Vec::with_capacity((num_ogs / 2) * num_species);
    let mut candidate_variances = if options.top_ogs.is_some() {
        Vec::with_capacity(num_ogs / 2)
    } else {
        Vec::new()
    };

    let mut discarded_invariant = 0;
    let mut discarded_low_prevalence = 0;
    let mut discarded_low_variance = 0;
    let mut discarded_low_range = 0;
    let mut discarded_low_spread = 0;

    let mut temp_buf = vec![0.0; num_species];
    let mut sort_buf = vec![0.0; num_species];

    for m in 0..num_ogs {
        let raw_slice = raw_table.get_og_counts(m);

        let mut species_with_gene = 0;
        let mut sum = 0.0;
        let mut min_val = f64::INFINITY;
        let mut max_val = f64::NEG_INFINITY;

        for (i, &raw_val) in raw_slice.iter().enumerate() {
            let val = if options.binary_mode {
                if raw_val > 0.0 { 1.0 } else { 0.0 }
            } else if options.log_count_mode {
                (raw_val + 1.0).log10()
            } else if let Some(ref totals) = options.species_total_genes {
                let total = totals.get(i).copied().unwrap_or(1.0).max(1.0);
                (raw_val / total) * 10_000.0
            } else {
                raw_val
            };

            temp_buf[i] = val;
            if val > 0.0 {
                species_with_gene += 1;
            }
            if val < min_val {
                min_val = val;
            }
            if val > max_val {
                max_val = val;
            }
            sum += val;
        }

        // Check 1: Invariant check (all 0s or all identical)
        if species_with_gene == 0 {
            discarded_invariant += 1;
            continue;
        }

        // Check 2: Minimum prevalence
        if species_with_gene < effective_min_species {
            discarded_low_prevalence += 1;
            continue;
        }

        // Check 3: Minimum absolute count range (max - min)
        if let Some(min_range) = options.min_count_range {
            if (max_val - min_val) < min_range {
                discarded_low_range += 1;
                continue;
            }
        }

        // Check 4: Minimum interdecile count spread (P95 - P05)
        if let Some(min_spread) = options.min_count_spread {
            sort_buf.copy_from_slice(&temp_buf);
            sort_buf.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let idx_05 = ((num_species - 1) as f64 * 0.05).round() as usize;
            let idx_95 = ((num_species - 1) as f64 * 0.95).round() as usize;
            let spread = sort_buf[idx_95] - sort_buf[idx_05];
            if spread < min_spread {
                discarded_low_spread += 1;
                continue;
            }
        }

        // Check 5: Minimum variance
        let mean = sum / n_f64;
        let mut var_sum = 0.0;
        for &v in &temp_buf {
            let diff = v - mean;
            var_sum += diff * diff;
        }
        let variance = var_sum / (n_f64 - 1.0).max(1.0);

        if variance < options.min_variance {
            discarded_low_variance += 1;
            continue;
        }

        // Retain candidate
        candidate_ids.push(raw_table.orthogroup_ids[m].clone());
        candidate_indices.push(m);
        candidate_counts.extend_from_slice(&temp_buf);
        if options.top_ogs.is_some() {
            candidate_variances.push(variance);
        }
    }

    // Step 6: Optional top_ogs selection
    let mut discarded_top_ogs = 0;
    let (orthogroup_ids, original_indices, counts) = if let Some(top_n) = options.top_ogs {
        let n_cands = candidate_ids.len();
        if n_cands > top_n {
            discarded_top_ogs = n_cands - top_n;
            // Rank candidate indices by variance descending
            let mut order: Vec<usize> = (0..n_cands).collect();
            order.sort_by(|&a, &b| {
                candidate_variances[b]
                    .partial_cmp(&candidate_variances[a])
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            order.truncate(top_n);
            // Restore original file order among top OGs for reproducibility
            order.sort_unstable();

            let mut final_ids = Vec::with_capacity(top_n);
            let mut final_indices = Vec::with_capacity(top_n);
            let mut final_counts = Vec::with_capacity(top_n * num_species);

            for idx in order {
                final_ids.push(candidate_ids[idx].clone());
                final_indices.push(candidate_indices[idx]);
                let start = idx * num_species;
                final_counts.extend_from_slice(&candidate_counts[start..start + num_species]);
            }
            (final_ids, final_indices, final_counts)
        } else {
            (candidate_ids, candidate_indices, candidate_counts)
        }
    } else {
        (candidate_ids, candidate_indices, candidate_counts)
    };

    let retained = orthogroup_ids.len();

    FilteredOrthoData {
        orthogroup_ids,
        original_indices,
        species_names: raw_table.species_names.clone(),
        counts,
        filter_report: OrthoFilterReport {
            initial_orthogroups: num_ogs,
            retained_orthogroups: retained,
            discarded_invariant,
            discarded_low_prevalence,
            discarded_low_variance,
            discarded_low_range,
            discarded_low_spread,
            discarded_top_ogs,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_filter_orthogroups_qc() {
        let table = OrthoGeneCountTable {
            orthogroup_ids: vec![
                "OG_all_zeros".into(),
                "OG_invariant_ones".into(),
                "OG_rare_1sp".into(),
                "OG_valid_cnv".into(),
            ],
            species_names: vec!["Sp1".into(), "Sp2".into(), "Sp3".into(), "Sp4".into()],
            counts: vec![
                0.0, 0.0, 0.0, 0.0, // all zeros
                1.0, 1.0, 1.0, 1.0, // all ones (zero variance)
                5.0, 0.0, 0.0, 0.0, // present in only 1 species
                2.0, 0.0, 3.0, 1.0, // valid variable OG
            ],
        };

        let options = OrthoFilterOptions {
            min_species_count: 2,
            min_variance: 1e-5,
            ..Default::default()
        };

        let filtered = filter_orthogroups(&table, &options);

        assert_eq!(filtered.num_orthogroups(), 1);
        assert_eq!(filtered.orthogroup_ids[0], "OG_valid_cnv");
        assert_eq!(filtered.original_indices[0], 3);
        assert_eq!(filtered.filter_report.initial_orthogroups, 4);
        assert_eq!(filtered.filter_report.retained_orthogroups, 1);
        assert_eq!(filtered.filter_report.discarded_invariant, 1);
        assert_eq!(filtered.filter_report.discarded_low_prevalence, 1);
        assert_eq!(filtered.filter_report.discarded_low_variance, 1);
    }

    #[test]
    fn test_filter_min_prevalence_and_range() {
        let table = OrthoGeneCountTable {
            orthogroup_ids: vec![
                "OG_50pct".into(),
                "OG_100pct_narrow".into(),
                "OG_100pct_wide".into(),
            ],
            species_names: vec!["Sp1".into(), "Sp2".into(), "Sp3".into(), "Sp4".into()],
            counts: vec![
                1.0, 1.0, 0.0, 0.0, // 50% prevalence
                1.0, 1.0, 2.0, 1.0, // 100% prevalence, range = 1
                1.0, 5.0, 10.0, 2.0, // 100% prevalence, range = 9
            ],
        };

        // Filter: at least 75% prevalence (3 of 4 species)
        let opts_prev = OrthoFilterOptions {
            min_prevalence: Some(0.75),
            min_variance: 1e-5,
            ..Default::default()
        };
        let res_prev = filter_orthogroups(&table, &opts_prev);
        assert_eq!(res_prev.num_orthogroups(), 2);
        assert_eq!(res_prev.filter_report.discarded_low_prevalence, 1);

        // Filter: range >= 5
        let opts_range = OrthoFilterOptions {
            min_count_range: Some(5.0),
            min_variance: 1e-5,
            ..Default::default()
        };
        let res_range = filter_orthogroups(&table, &opts_range);
        assert_eq!(res_range.num_orthogroups(), 1);
        assert_eq!(res_range.orthogroup_ids[0], "OG_100pct_wide");
        assert_eq!(res_range.filter_report.discarded_low_range, 2);
    }

    #[test]
    fn test_filter_top_ogs() {
        let table = OrthoGeneCountTable {
            orthogroup_ids: vec![
                "OG_low_var".into(),
                "OG_med_var".into(),
                "OG_high_var".into(),
            ],
            species_names: vec!["Sp1".into(), "Sp2".into(), "Sp3".into()],
            counts: vec![
                1.0, 2.0, 1.0, // variance ~ 0.33
                1.0, 4.0, 1.0, // variance ~ 3.0
                1.0, 10.0, 1.0, // variance ~ 27.0
            ],
        };

        let opts_top = OrthoFilterOptions {
            min_species_count: 1,
            top_ogs: Some(2),
            ..Default::default()
        };
        let res_top = filter_orthogroups(&table, &opts_top);
        assert_eq!(res_top.num_orthogroups(), 2);
        assert_eq!(res_top.orthogroup_ids, vec!["OG_med_var", "OG_high_var"]);
        assert_eq!(res_top.filter_report.discarded_top_ogs, 1);
    }
}
