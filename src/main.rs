use clap::builder::styling::{AnsiColor, Effects, Styles};
use clap::{Args, Parser, Subcommand};
use colored::Colorize;
use gelimo::genomics::{run_ortho_pipeline, OrthoPipelineConfig, OrthoPipelineSummary};
use gelimo::math::GlsSolverMode;
use gelimo::pipeline::{run_pgls_pipeline, PglsConfig, PglsSummary};

fn cli_styles() -> Styles {
    Styles::styled()
        .header(AnsiColor::BrightCyan.on_default() | Effects::BOLD)
        .usage(AnsiColor::BrightYellow.on_default() | Effects::BOLD)
        .literal(AnsiColor::BrightGreen.on_default() | Effects::BOLD)
        .placeholder(AnsiColor::Cyan.on_default())
        .valid(AnsiColor::BrightGreen.on_default())
        .invalid(AnsiColor::Red.on_default() | Effects::BOLD)
        .error(AnsiColor::Red.on_default() | Effects::BOLD)
}

/// gelimo: Generalized Phylogenetic Linear Models Engine in Rust
#[derive(Parser, Debug)]
#[command(
    name = "gelimo",
    author = "Hector Jesus Lopez Ordaz <hectorjl@lcg.unam.mx>",
    version = "0.1.0",
    about = "Generalized Phylogenetic Linear Models and High-Throughput Genomic Screening Engine in Rust",
    long_about = "High-performance PGLS engine compatible with R caper and optimized for genome-wide orthogroup analysis with Rayon.",
    styles = cli_styles()
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Fit a single PGLS model on comparative species traits (caper compatible)
    #[command(alias = "pgls", styles = cli_styles())]
    Fit(FitArgs),

    /// Run high-throughput parallel PGLS screening on OrthoFinder GeneCount tables
    #[command(alias = "gwas", styles = cli_styles())]
    Orthogroups(OrthoArgs),
}

#[derive(Args, Debug)]
struct FitArgs {
    /// Path to Newick phylogenetic tree file
    #[arg(short, long)]
    tree: String,

    /// Path to CSV/TSV species trait data file
    #[arg(short, long)]
    data: String,

    /// Response variable column name (Y)
    #[arg(short = 'y', long)]
    response: String,

    /// Predictor variable column names (X), comma-separated
    #[arg(short = 'x', long, value_delimiter = ',')]
    predictors: Vec<String>,

    /// Pagel's lambda transformation parameter (default: 1.0)
    #[arg(short = 'l', long, default_value_t = 1.0)]
    lambda: f64,

    /// Optimize Pagel's lambda by Maximum Likelihood (ML)
    #[arg(long, default_value_t = false)]
    opt_lambda: bool,

    /// Enable automatic tree and table pruning for non-matching species
    #[arg(short = 'p', long, default_value_t = false)]
    prune: bool,

    /// Solver algorithm mode: "dense" (O(N^3) Cholesky) or "linear" (O(N) Ho & Ane 2014)
    #[arg(short = 's', long, default_value = "dense")]
    solver: String,
}

#[derive(Args, Debug)]
struct OrthoArgs {
    /// Path to Newick phylogenetic tree file
    #[arg(short, long)]
    tree: String,

    /// Path to CSV/TSV species trait data file (phenotype + covariates)
    #[arg(short, long)]
    data: String,

    /// Path to OrthoFinder Orthogroups.GeneCount.tsv file
    #[arg(short = 'g', long)]
    orthogroups: String,

    /// Response variable column name in phenotype table (Y)
    #[arg(short = 'y', long)]
    response: String,

    /// Optional covariate column names (e.g. Log Adult weight), comma-separated
    #[arg(short = 'x', long, value_delimiter = ',')]
    covariates: Vec<String>,

    /// Output TSV file path for full genome-wide PGLS results
    #[arg(short = 'o', long, default_value = "orthogroups_pgls_results.tsv")]
    output: String,

    /// Pagel's lambda transformation parameter (default: 1.0)
    #[arg(short = 'l', long, default_value_t = 1.0)]
    lambda: f64,

    /// Optimize Pagel's lambda by Maximum Likelihood (ML) on the baseline trait model
    #[arg(long, default_value_t = false)]
    opt_lambda: bool,

    /// Minimum species count possessing at least 1 gene copy (QC filter)
    #[arg(long, default_value_t = 5)]
    min_species: usize,

    /// Minimum variance of gene counts across species (QC filter)
    #[arg(long, default_value_t = 1e-6)]
    min_variance: f64,

    /// Convert gene counts to binary presence/absence (0 vs 1)
    #[arg(long, default_value_t = false)]
    binary: bool,

    /// Apply log10(count + 1) transformation to gene counts
    #[arg(long, default_value_t = false)]
    log_count: bool,

    /// Enable automatic tree and table pruning for non-matching species
    #[arg(short = 'p', long, default_value_t = true)]
    prune: bool,
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Fit(args) => run_fit_command(args),
        Commands::Orthogroups(args) => run_ortho_command(args),
    }
}

fn run_fit_command(args: FitArgs) {
    let solver_mode = match args.solver.to_lowercase().as_str() {
        "linear" | "ho_ane" => GlsSolverMode::LinearHoAne,
        _ => GlsSolverMode::DenseCholesky,
    };

    let config = PglsConfig {
        tree_path: args.tree,
        data_path: args.data,
        response_col: args.response,
        predictor_cols: args.predictors,
        lambda: args.lambda,
        optimize_lambda: args.opt_lambda,
        prune: args.prune,
        solver_mode,
    };

    match run_pgls_pipeline(&config) {
        Ok(summary) => print_fit_summary(&config, &summary),
        Err(err) => {
            eprintln!("{} {}", "Error running PGLS model:".red().bold(), err);
            std::process::exit(1);
        }
    }
}

fn run_ortho_command(args: OrthoArgs) {
    let config = OrthoPipelineConfig {
        tree_path: args.tree,
        data_path: args.data,
        ortho_path: args.orthogroups,
        response_col: args.response,
        covariate_cols: args.covariates,
        output_path: args.output,
        lambda: args.lambda,
        optimize_lambda: args.opt_lambda,
        min_species_count: args.min_species,
        min_variance: args.min_variance,
        binary_mode: args.binary,
        log_count_mode: args.log_count,
        prune: args.prune,
        fuzzy_taxonomy: true,
    };

    println!("{}", "======================================================================".cyan().bold());
    println!("{}", "gelimo v0.1.0 - High-Throughput Genomic Orthogroups PGLS Engine".bright_white().bold());
    println!("{}", "======================================================================".cyan().bold());
    println!("{:<18} {}", "Tree file:".bright_cyan().bold(), config.tree_path);
    println!("{:<18} {}", "Phenotype file:".bright_cyan().bold(), config.data_path);
    println!("{:<18} {}", "OrthoFinder file:".bright_cyan().bold(), config.ortho_path);
    println!("{:<18} {}", "Response (Y):".bright_cyan().bold(), config.response_col.yellow().bold());
    if !config.covariate_cols.is_empty() {
        println!("{:<18} {}", "Covariates:".bright_cyan().bold(), config.covariate_cols.join(", ").yellow());
    }
    println!("{:<18} {}", "Output TSV:".bright_cyan().bold(), config.output_path.green().bold());
    println!(
        "{:<18} {}",
        "Mode:".bright_cyan().bold(),
        if config.log_count_mode {
            "Log10 Count (log10(count + 1))".bright_magenta().bold()
        } else if config.binary_mode {
            "Binary (Presence/Absence)".magenta().bold()
        } else {
            "Continuous (Gene Count / CNV)".bright_magenta().bold()
        }
    );
    println!();
    println!("{}", "Starting tripartite alignment, QC filtering and parallel PGLS...".bright_black().italic());

    let start = std::time::Instant::now();
    match run_ortho_pipeline(&config) {
        Ok(summary) => {
            let elapsed = start.elapsed();
            print_ortho_summary(&summary, elapsed);
        }
        Err(err) => {
            eprintln!("{} {}", "Error in genomic PGLS pipeline:".red().bold(), err);
            std::process::exit(1);
        }
    }
}

fn format_p_value(p: f64) -> colored::ColoredString {
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

fn print_fit_summary(config: &PglsConfig, summary: &PglsSummary) {
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

fn print_ortho_summary(summary: &OrthoPipelineSummary, elapsed: std::time::Duration) {
    println!();
    println!("{}", "--- Alignment & QC Filtering Summary ---".bright_yellow().bold());
    println!("{:<28} {}", "Matched species count:".bright_white(), summary.tri_report.matched_species_count.to_string().green().bold());
    println!("{:<28} {}", "Taxonomy fuzzy matches:".bright_white(), summary.tri_report.taxonomy_matches.len().to_string().green());
    println!("{:<28} {}", "Total raw orthogroups:".bright_white(), summary.filter_report.initial_orthogroups);
    println!("{:<28} {}", "Retained QC orthogroups:".bright_white(), summary.filter_report.retained_orthogroups.to_string().bright_green().bold());
    println!("{:<28} {}", "Discarded invariant:".bright_white(), summary.filter_report.discarded_invariant.to_string().yellow());
    println!("{:<28} {}", "Discarded low prevalence:".bright_white(), summary.filter_report.discarded_low_prevalence.to_string().yellow());
    println!("{:<28} {}", "Discarded low variance:".bright_white(), summary.filter_report.discarded_low_variance.to_string().yellow());
    println!("{:<28} {}", "Pagel's lambda applied:".bright_white(), format!("{:.4}", summary.lambda_used).cyan().bold());
    println!("{:<28} {}", "Total execution time:".bright_white(), format!("{:.2?}", elapsed).bright_green().bold());
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
            "{:<12} {:>10} {:>10} {:>8} {:>12} {:>14} {:>8} {:>10}",
            "Orthogroup", "Beta_OG", "Std. Error", "t value", "p-value", "FDR q-val", "R2", "Species"
        )
        .bright_cyan()
        .bold()
    );
    println!("{}", "-".repeat(92).bright_black());

    for hit in &summary.top_hits {
        let q_str = format!("{:>14.4e}", hit.q_value_fdr);
        let q_colored = if hit.q_value_fdr < 0.01 {
            q_str.bright_green().bold()
        } else if hit.q_value_fdr < 0.05 {
            q_str.green()
        } else {
            q_str.yellow()
        };

        println!(
            "{:<12} {:>10} {:>10} {:>8} {:>12.4e} {:>14} {:>8} {:>6}/{}",
            hit.orthogroup_id.bright_yellow().bold(),
            format!("{:>10.4}", hit.beta_og).cyan().bold(),
            format!("{:>10.4}", hit.se_og).white(),
            format!("{:>8.3}", hit.t_stat).bright_cyan(),
            hit.p_value,
            q_colored,
            format!("{:>8.4}", hit.r_squared).bright_magenta().bold(),
            hit.species_with_gene.to_string().green(),
            summary.tri_report.matched_species_count
        );
    }
    println!("{}", "======================================================================".cyan().bold());
}
