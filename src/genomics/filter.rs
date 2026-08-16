//! Quality Control (QC) filtering and transformations for OrthoFinder orthogroups.

use crate::genomics::ortho_reader::OrthoGeneCountTable;

/// Configuration options for filtering and transforming orthogroup gene counts.
#[derive(Debug, Clone)]
pub struct OrthoFilterOptions {
    /// Minimum number of species that must possess at least 1 gene copy.
    pub min_species_count: usize,
    /// Minimum variance of gene counts across species.
    pub min_variance: f64,
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
            min_species_count: 5,
            min_variance: 1e-7,
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

    let mut orthogroup_ids = Vec::with_capacity(num_ogs / 2);
    let mut original_indices = Vec::with_capacity(num_ogs / 2);
    let mut counts = Vec::with_capacity((num_ogs / 2) * num_species);

    let mut discarded_invariant = 0;
    let mut discarded_low_prevalence = 0;
    let mut discarded_low_variance = 0;

    let mut temp_buf = vec![0.0; num_species];

    for m in 0..num_ogs {
        let raw_slice = raw_table.get_og_counts(m);

        // Apply optional transformation to temp buffer
        let mut species_with_gene = 0;
        let mut sum = 0.0;

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
            sum += val;
        }

        // Check 1: Invariant check (all 0s or all identical)
        if species_with_gene == 0 {
            discarded_invariant += 1;
            continue;
        }

        // Check 2: Minimum prevalence
        if species_with_gene < options.min_species_count {
            discarded_low_prevalence += 1;
            continue;
        }

        // Check 3: Minimum variance
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

        // Retain orthogroup
        orthogroup_ids.push(raw_table.orthogroup_ids[m].clone());
        original_indices.push(m);
        counts.extend_from_slice(&temp_buf);
    }

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
            binary_mode: false,
            log_count_mode: false,
            species_total_genes: None,
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
}
