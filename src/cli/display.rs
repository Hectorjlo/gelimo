//! Terminal display, table formatting, and colorized reports for `gelimo`.

use colored::Colorize;
use crate::genomics::pipeline::OrthoPipelineSummary;
use crate::pipeline::{PglsConfig, PglsSummary};

pub fn format_p_value(p: f64) -> colored::ColoredString {
    let s = format!("{:.4e}", p);
    if p < 0.001 {
        format!("{} ***", s).bright_green().bold()
    } else if p < 0.01 {
        format!("{} **", s).green().bold()
    } else if p < 0.05 {
        format!("{} *", s).green()
    } else {
        s.dimmed()
    }
}

pub fn print_fit_summary(config: &PglsConfig, summary: &PglsSummary) {
    println!("{}", "======================================================================".cyan().bold());
    println!("{}", "gelimo v0.1.0 - Generalized Phylogenetic Linear Models Engine".bright_white().bold());
    println!("{}", "======================================================================".cyan().bold());
    println!("{:<16} {}", "Tree file:".bright_cyan().bold(), config.tree_path);
    println!("{:<16} {}", "Data file:".bright_cyan().bold(), config.data_path);
    println!("{:<16} {}", "Response (Y):".bright_cyan().bold(), config.response_col.yellow().bold());
    println!("{:<16} {}", "Predictors:".bright_cyan().bold(), config.predictor_cols.join(", ").yellow());
    if config.optimize_lambda {
        println!(
            "{:<16} {}",
            "Pagel Lambda:".bright_cyan().bold(),
            format!("{:.4} (Optimized by ML)", summary.lambda_used).bright_green().bold()
        );
    } else {
        println!(
            "{:<16} {}",
            "Pagel Lambda:".bright_cyan().bold(),
            format!("{:.4} (Fixed)", summary.lambda_used).cyan()
        );
    }
    println!();
    println!("{}", "--- Data Alignment & Pruning Summary ---".bright_yellow().bold());
    println!("{:<24} {}", "Matched species count:".bright_white(), summary.comp_data.report.matched_species_count.to_string().green().bold());
    if !summary.comp_data.report.taxonomy_matches.is_empty() {
        println!("{:<24} {}", "Taxonomy matches:".bright_white(), summary.comp_data.report.taxonomy_matches.len().to_string().green());
        for m in &summary.comp_data.report.taxonomy_matches {
            println!("  * {}", m.green());
        }
    }
    if !summary.comp_data.report.pruned_tree_tips.is_empty() {
        println!("{:<24} {}", "Pruned tree tips:".bright_white(), summary.comp_data.report.pruned_tree_tips.len().to_string().yellow());
    }
    if !summary.comp_data.report.pruned_table_rows.is_empty() {
        println!("{:<24} {}", "Pruned table rows:".bright_white(), summary.comp_data.report.pruned_table_rows.len().to_string().yellow());
    }
    println!();
    println!("{}", "--- PGLS Model Coefficients ---".bright_yellow().bold());
    println!(
        "{}",
        format!(
            "{:<22} {:>12} {:>12} {:>12} {:>16}",
            "Parameter", "Estimate", "Std. Error", "t value", "Pr(>|t|)"
        )
        .bright_cyan()
        .bold()
    );
    println!("{}", "-".repeat(78).bright_black());

    let param_names = std::iter::once("(Intercept)".to_string())
        .chain(config.predictor_cols.clone());

    for (i, p_name) in param_names.enumerate() {
        let b = summary.gls_result.beta[(i, 0)];
        let se = summary.gls_result.se[i];
        let t_val = if se > 0.0 { b / se } else { 0.0 };
        let p_val = summary.p_values[i];

        println!(
            "{:<22} {:>12} {:>12} {:>12} {:>18}",
            p_name.bright_white().bold(),
            format!("{:>12.4}", b).cyan().bold(),
            format!("{:>12.4}", se).white(),
            format!("{:>12.4}", t_val).bright_cyan(),
            format_p_value(p_val)
        );
    }

    println!();
    println!("{}", "--- Model Fit Statistics ---".bright_yellow().bold());
    println!("{:<32} {:>12.4}", "Residual Variance (sigma2 ML):".bright_white(), summary.gls_result.sigma2_ml);
    println!("{:<32} {:>12.4}", "Log-Likelihood:".bright_white(), summary.gls_result.log_likelihood);
    println!("{:<32} {}", "Multiple R-squared:".bright_white(), format!("{:>12.4}", summary.r_squared).bright_magenta().bold());
    println!("{:<32} {:>12.4}", "AIC:".bright_white(), summary.aic);
    println!("{:<32} {:>12?}", "Solver Algorithm:".bright_white(), config.solver_mode);

    if let Some(ref opt) = summary.lambda_opt {
        println!();
        println!("{}", "--- Pagel Lambda Likelihood Ratio Tests ---".bright_yellow().bold());
        println!("{:<32} {}", "Estimated ML Lambda:".bright_white(), format!("{:>12.4}", opt.optimal_lambda).bright_green().bold());
        println!("{:<32} {:>12.4}", "Log-Likelihood (ML Lambda):".bright_white(), opt.log_likelihood_ml);
        println!("{:<32} {:>12.4}", "Log-Likelihood (Lambda = 0):".bright_white(), opt.log_likelihood_lambda0);
        println!("{:<32} {:>12.4}", "Log-Likelihood (Lambda = 1):".bright_white(), opt.log_likelihood_lambda1);
        println!("{:<32} {:>20}", "LRT p-value vs Lambda = 0.0:".bright_white(), format_p_value(opt.lrt_p_value_vs_0));
        println!("{:<32} {:>20}", "LRT p-value vs Lambda = 1.0:".bright_white(), format_p_value(opt.lrt_p_value_vs_1));
    }
    println!("{}", "======================================================================".cyan().bold());
}

pub fn print_ortho_summary(summary: &OrthoPipelineSummary, elapsed: std::time::Duration) {
    let raw_ogs = summary.filter_report.initial_orthogroups;
    let retained = summary.filter_report.retained_orthogroups;
    let pct = |count: usize| -> String {
        if raw_ogs > 0 {
            format!("{:>6.2}%", (count as f64 / raw_ogs as f64) * 100.0)
        } else {
            "  0.00%".to_string()
        }
    };

    println!();
    println!("{}", "--- Alignment & QC Filtering Summary ---".bright_yellow().bold());
    println!("{:<32} {}", "Input Phenotype Species:".bright_white(), summary.tri_report.pheno_rows_count.to_string().cyan());
    println!("{:<32} {}", "Input Tree Tips:".bright_white(), summary.tri_report.tree_tips_count.to_string().cyan());
    println!("{:<32} {}", "Input OrthoFinder Species:".bright_white(), summary.tri_report.orthogroup_species_count.to_string().cyan());
    println!("{:<32} {}", "Matched / Aligned Species:".bright_white(), format!("{} (Fuzzy matches: {})", summary.tri_report.matched_species_count, summary.tri_report.taxonomy_matches.len()).green().bold());
    println!();
    println!("{:<32} {}", "Total Raw Orthogroups:".bright_white(), raw_ogs.to_string().cyan().bold());
    println!("{:<32} {}", "Retained QC Orthogroups:".bright_white(), format!("{} ({})", retained, pct(retained)).bright_green().bold());
    println!("{:<32} {}", "  - Discarded invariant:".bright_white(), format!("{} ({})", summary.filter_report.discarded_invariant, pct(summary.filter_report.discarded_invariant)).yellow());
    println!("{:<32} {}", "  - Discarded low prevalence:".bright_white(), format!("{} ({})", summary.filter_report.discarded_low_prevalence, pct(summary.filter_report.discarded_low_prevalence)).yellow());
    println!("{:<32} {}", "  - Discarded low variance:".bright_white(), format!("{} ({})", summary.filter_report.discarded_low_variance, pct(summary.filter_report.discarded_low_variance)).yellow());
    if summary.filter_report.discarded_low_range > 0 {
        println!("{:<32} {}", "  - Discarded low range:".bright_white(), format!("{} ({})", summary.filter_report.discarded_low_range, pct(summary.filter_report.discarded_low_range)).yellow());
    }
    if summary.filter_report.discarded_low_spread > 0 {
        println!("{:<32} {}", "  - Discarded low spread:".bright_white(), format!("{} ({})", summary.filter_report.discarded_low_spread, pct(summary.filter_report.discarded_low_spread)).yellow());
    }
    if summary.filter_report.discarded_top_ogs > 0 {
        println!("{:<32} {}", "  - Discarded by top-N limit:".bright_white(), format!("{} ({})", summary.filter_report.discarded_top_ogs, pct(summary.filter_report.discarded_top_ogs)).yellow());
    }
    println!("{:<32} {}", "Pagel's lambda applied:".bright_white(), format!("{:.4}", summary.lambda_used).cyan().bold());
    println!("{:<32} {}", "Total execution time:".bright_white(), format!("{:.2?}", elapsed).bright_green().bold());

    if !summary.eval_ranges_summary.is_empty() {
        println!();
        println!("{}", "--- Multi-Threshold Power & FDR Summary ---".bright_yellow().bold());
        println!(
            "{}",
            format!(
                "{:<22} {:>14} {:>24} {:>24}",
                "Threshold (Spread >=)", "Retained OGs", "Signif (FDR q < 0.05)", "Signif (FDR q < 0.01)"
            )
            .bright_cyan()
            .bold()
        );
        println!("{}", "-".repeat(88).bright_black());
        for entry in &summary.eval_ranges_summary {
            println!(
                "{:<22} {:>14} {:>24} {:>24}",
                format!("Spread >= {:.1}", entry.threshold_spread).bright_white().bold(),
                entry.retained_ogs.to_string().cyan(),
                entry.significant_fdr_05.to_string().bright_green().bold(),
                entry.significant_fdr_01.to_string().green()
            );
        }
    }

    println!();
    println!("{}", "--- Multiple Testing Significance Summary ---".bright_yellow().bold());
    println!("{:<40} {}", "Significant OGs (FDR q < 0.05):".bright_white(), summary.significant_fdr_05.to_string().bright_green().bold());
    println!("{:<40} {}", "Significant OGs (FDR q < 0.01):".bright_white(), summary.significant_fdr_01.to_string().bright_green().bold());
    println!("{:<40} {}", "Significant OGs (Bonferroni p < 0.05):".bright_white(), summary.significant_bonferroni.to_string().bright_green().bold());
    println!();
    println!("{}", "--- Top 10 Genomic Hits (Sorted by FDR q-value) ---".bright_yellow().bold());
    println!(
        "{}",
        format!(
            "{:<12} {:>10} {:>10} {:>8} {:>12} {:>14} {:>10} {:>8} {:>10}",
            "Orthogroup", "Beta_OG", "Std. Error", "t value", "p-value", "FDR q-val", "Partial_r", "R2", "Species"
        )
        .bright_cyan()
        .bold()
    );
    println!("{}", "-".repeat(104).bright_black());

    for hit in &summary.top_hits {
        println!(
            "{:<12} {:>10.4} {:>10.4} {:>8.3} {:>12} {:>14} {:>10.4} {:>8.4} {:>10}",
            hit.orthogroup_id.bright_white().bold(),
            hit.beta_og,
            hit.se_og,
            hit.t_stat,
            format_p_value(hit.p_value),
            format_p_value(hit.q_value_fdr),
            hit.partial_r,
            hit.r_squared,
            format!("{}/{}", hit.species_with_gene, summary.tri_report.matched_species_count).cyan()
        );
    }
    println!("{}", "======================================================================".cyan().bold());
}
