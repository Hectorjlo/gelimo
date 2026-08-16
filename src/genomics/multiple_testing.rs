//! Multiple testing corrections (Benjamini-Hochberg FDR, Bonferroni) and TSV export.

use std::fs::File;
use std::io::{BufWriter, Write};

use crate::genomics::engine::OrthoPglsResult;

/// Fully annotated orthogroup PGLS result with adjusted multiple-testing significance statistics.
#[derive(Debug, Clone)]
pub struct AnnotatedOrthoResult {
    pub orthogroup_id: String,
    pub beta_og: f64,
    pub se_og: f64,
    pub t_stat: f64,
    pub p_value: f64,
    pub q_value_fdr: f64,
    pub p_bonferroni: f64,
    pub partial_r: f64,
    pub r_squared: f64,
    pub sigma2: f64,
    pub mean_copies: f64,
    pub species_with_gene: usize,
}

/// Applies False Discovery Rate (Benjamini-Hochberg) and Bonferroni adjustments and sorts results by significance.
// ponytail: O(M log M) in-place index sort and backward cumulative minimum pass
pub fn adjust_p_values(raw_results: Vec<OrthoPglsResult>) -> Vec<AnnotatedOrthoResult> {
    let m = raw_results.len();
    if m == 0 {
        return Vec::new();
    }

    let m_f64 = m as f64;

    // 1. Sort indices by p-value ascending
    let mut indices: Vec<usize> = (0..m).collect();
    indices.sort_unstable_by(|&i, &j| {
        raw_results[i]
            .p_value
            .partial_cmp(&raw_results[j].p_value)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut q_values = vec![0.0; m];
    let mut bonferroni_values = vec![0.0; m];

    // 2. Bonferroni adjustment
    for &idx in &indices {
        bonferroni_values[idx] = (raw_results[idx].p_value * m_f64).min(1.0);
    }

    // 3. Benjamini-Hochberg FDR adjustment: q_i = min_{j >= i} ( (M / rank_j) * p_j )
    let mut running_min = 1.0_f64;
    for rank in (1..=m).rev() {
        let idx = indices[rank - 1];
        let p_val = raw_results[idx].p_value;
        let candidate = (p_val * m_f64 / (rank as f64)).min(1.0);
        running_min = running_min.min(candidate);
        q_values[idx] = running_min;
    }

    // 4. Assemble sorted AnnotatedOrthoResult list
    let mut annotated = Vec::with_capacity(m);
    for idx in indices {
        let res = raw_results[idx].clone();
        annotated.push(AnnotatedOrthoResult {
            orthogroup_id: res.orthogroup_id,
            beta_og: res.beta_og,
            se_og: res.se_og,
            t_stat: res.t_stat,
            p_value: res.p_value,
            q_value_fdr: q_values[idx],
            p_bonferroni: bonferroni_values[idx],
            partial_r: res.partial_r,
            r_squared: res.r_squared,
            sigma2: res.sigma2,
            mean_copies: res.mean_copies,
            species_with_gene: res.species_with_gene,
        });
    }

    annotated
}

/// Exports annotated orthogroup PGLS results to a TSV file using buffered I/O.
// ponytail: streaming buffered TSV writer
pub fn export_ortho_results_tsv(
    results: &[AnnotatedOrthoResult],
    output_path: &str,
) -> Result<(), String> {
    let file = File::create(output_path)
        .map_err(|e| format!("Failed to create output file '{}': {}", output_path, e))?;
    let mut writer = BufWriter::with_capacity(64 * 1024, file);

    // Header line
    writeln!(
        writer,
        "Orthogroup\tBeta_OG\tStdError\tt_stat\tp_value\tq_value_fdr\tp_bonferroni\tPartial_r\tR_squared\tSigma2\tMean_Copies\tSpecies_With_Gene"
    )
    .map_err(|e| e.to_string())?;

    for res in results {
        writeln!(
            writer,
            "{}\t{:.6}\t{:.6}\t{:.4}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.4}\t{:.4}\t{:.4}\t{:.3}\t{}",
            res.orthogroup_id,
            res.beta_og,
            res.se_og,
            res.t_stat,
            res.p_value,
            res.q_value_fdr,
            res.p_bonferroni,
            res.partial_r,
            res.r_squared,
            res.sigma2,
            res.mean_copies,
            res.species_with_gene
        )
        .map_err(|e| e.to_string())?;
    }

    writer.flush().map_err(|e| e.to_string())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_adjust_p_values_fdr_monotonicity() {
        let raw = vec![
            OrthoPglsResult {
                orthogroup_id: "OG1".into(),
                beta_og: 1.0,
                se_og: 0.1,
                t_stat: 10.0,
                p_value: 0.001,
                partial_r: 0.5,
                r_squared: 0.5,
                sigma2: 1.0,
                mean_copies: 2.0,
                species_with_gene: 10,
            },
            OrthoPglsResult {
                orthogroup_id: "OG2".into(),
                beta_og: 0.5,
                se_og: 0.2,
                t_stat: 2.5,
                p_value: 0.05,
                partial_r: 0.2,
                r_squared: 0.2,
                sigma2: 1.0,
                mean_copies: 1.0,
                species_with_gene: 8,
            },
            OrthoPglsResult {
                orthogroup_id: "OG3".into(),
                beta_og: 0.1,
                se_og: 0.2,
                t_stat: 0.5,
                p_value: 0.6,
                partial_r: 0.05,
                r_squared: 0.01,
                sigma2: 1.0,
                mean_copies: 1.0,
                species_with_gene: 5,
            },
        ];

        let adjusted = adjust_p_values(raw);

        assert_eq!(adjusted.len(), 3);
        assert_eq!(adjusted[0].orthogroup_id, "OG1");
        assert!(adjusted[0].q_value_fdr <= adjusted[1].q_value_fdr);
        assert!(adjusted[1].q_value_fdr <= adjusted[2].q_value_fdr);
        assert!(adjusted[0].p_bonferroni <= adjusted[1].p_bonferroni);
        assert_eq!(adjusted[0].partial_r, 0.5);
    }
}
