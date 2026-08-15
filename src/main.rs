use clap::{Args, Parser, Subcommand};
use gelimo::genomics::{run_ortho_pipeline, OrthoPipelineConfig, OrthoPipelineSummary};
use gelimo::math::GlsSolverMode;
use gelimo::pipeline::{run_pgls_pipeline, PglsConfig, PglsSummary};

/// gelimo: Generalized Phylogenetic Linear Models Engine in Rust
#[derive(Parser, Debug)]
#[command(
    name = "gelimo",
    author = "Hector Rodriguez",
    version = "0.1.0",
    about = "Generalized Phylogenetic Linear Models and High-Throughput Genomic Screening Engine in Rust",
    long_about = "High-performance PGLS engine compatible with R caper and optimized for genome-wide orthogroup analysis with Rayon."
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Fit a single PGLS model on comparative species traits (caper compatible)
    #[command(alias = "pgls")]
    Fit(FitArgs),

    /// Run high-throughput parallel PGLS screening on OrthoFinder GeneCount tables
    #[command(alias = "gwas")]
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
            eprintln!("Error running PGLS model: {}", err);
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
        prune: args.prune,
        fuzzy_taxonomy: true,
    };

    println!("======================================================================");
    println!("gelimo v0.1.0 - High-Throughput Genomic Orthogroups PGLS Engine");
    println!("======================================================================");
    println!("Tree file:        {}", config.tree_path);
    println!("Phenotype file:   {}", config.data_path);
    println!("OrthoFinder file: {}", config.ortho_path);
    println!("Response (Y):     {}", config.response_col);
    if !config.covariate_cols.is_empty() {
        println!("Covariates:       {}", config.covariate_cols.join(", "));
    }
    println!("Output TSV:       {}", config.output_path);
    println!("Mode:             {}", if config.binary_mode { "Binary (Presence/Absence)" } else { "Continuous (Gene Count / CNV)" });
    println!();
    println!("Starting tripartite alignment, QC filtering and parallel PGLS...");

    let start = std::time::Instant::now();
    match run_ortho_pipeline(&config) {
        Ok(summary) => {
            let elapsed = start.elapsed();
            print_ortho_summary(&summary, elapsed);
        }
        Err(err) => {
            eprintln!("Error in genomic PGLS pipeline: {}", err);
            std::process::exit(1);
        }
    }
}

fn print_fit_summary(config: &PglsConfig, summary: &PglsSummary) {
    println!("======================================================================");
    println!("gelimo v0.1.0 - Generalized Phylogenetic Linear Models Engine");
    println!("======================================================================");
    println!("Tree file:  {}", config.tree_path);
    println!("Data file:  {}", config.data_path);
    println!("Response:   {}", config.response_col);
    println!("Predictors: {}", config.predictor_cols.join(", "));
    if config.optimize_lambda {
        println!("Pagel Lambda: {:.4} (Optimized by Maximum Likelihood)", summary.lambda_used);
    } else {
        println!("Pagel Lambda: {:.4} (Fixed)", summary.lambda_used);
    }
    println!();
    println!("--- Data Alignment & Pruning Summary ---");
    println!("Matched species count: {}", summary.comp_data.report.matched_species_count);
    if !summary.comp_data.report.taxonomy_matches.is_empty() {
        println!("Taxonomy matches:      {}", summary.comp_data.report.taxonomy_matches.len());
        for m in &summary.comp_data.report.taxonomy_matches {
            println!("  * {}", m);
        }
    }
    if !summary.comp_data.report.pruned_tree_tips.is_empty() {
        println!("Pruned tree tips:      {}", summary.comp_data.report.pruned_tree_tips.len());
    }
    if !summary.comp_data.report.pruned_table_rows.is_empty() {
        println!("Pruned table rows:     {}", summary.comp_data.report.pruned_table_rows.len());
    }
    println!();
    println!("--- PGLS Model Coefficients ---");
    println!(
        "{:<22} {:>12} {:>12} {:>12} {:>12}",
        "Parameter", "Estimate", "Std. Error", "t value", "Pr(>|t|)"
    );
    println!("{}", "-".repeat(70));

    let param_names = std::iter::once("(Intercept)".to_string())
        .chain(config.predictor_cols.clone());

    for (i, p_name) in param_names.enumerate() {
        let b = summary.gls_result.beta[(i, 0)];
        let se = summary.gls_result.se[i];
        let t_val = if se > 0.0 { b / se } else { 0.0 };
        let p_val = summary.p_values[i];

        println!(
            "{:<22} {:>12.4} {:>12.4} {:>12.4} {:>12.4e}",
            p_name, b, se, t_val, p_val
        );
    }

    println!();
    println!("--- Model Fit Statistics ---");
    println!("Residual Variance (sigma2 ML): {:>12.4}", summary.gls_result.sigma2_ml);
    println!("Log-Likelihood:                {:>12.4}", summary.gls_result.log_likelihood);
    println!("Multiple R-squared:            {:>12.4}", summary.r_squared);
    println!("AIC:                           {:>12.4}", summary.aic);
    println!("Solver Algorithm:              {:?}", config.solver_mode);

    if let Some(ref opt) = summary.lambda_opt {
        println!();
        println!("--- Pagel Lambda Likelihood Ratio Tests ---");
        println!("Estimated ML Lambda:           {:>12.4}", opt.optimal_lambda);
        println!("Log-Likelihood (ML Lambda):    {:>12.4}", opt.log_likelihood_ml);
        println!("Log-Likelihood (Lambda = 0):   {:>12.4}", opt.log_likelihood_lambda0);
        println!("Log-Likelihood (Lambda = 1):   {:>12.4}", opt.log_likelihood_lambda1);
        println!("LRT p-value vs Lambda = 0.0:   {:>12.4e}", opt.lrt_p_value_vs_0);
        println!("LRT p-value vs Lambda = 1.0:   {:>12.4e}", opt.lrt_p_value_vs_1);
    }
    println!("======================================================================");
}

fn print_ortho_summary(summary: &OrthoPipelineSummary, elapsed: std::time::Duration) {
    println!();
    println!("--- Alignment & QC Filtering Summary ---");
    println!("Matched species count:     {}", summary.tri_report.matched_species_count);
    println!("Taxonomy fuzzy matches:    {}", summary.tri_report.taxonomy_matches.len());
    println!("Total raw orthogroups:     {}", summary.filter_report.initial_orthogroups);
    println!("Retained QC orthogroups:   {}", summary.filter_report.retained_orthogroups);
    println!("Discarded invariant:       {}", summary.filter_report.discarded_invariant);
    println!("Discarded low prevalence:  {}", summary.filter_report.discarded_low_prevalence);
    println!("Discarded low variance:    {}", summary.filter_report.discarded_low_variance);
    println!("Pagel's lambda applied:    {:.4}", summary.lambda_used);
    println!("Total execution time:      {:.2?}", elapsed);
    println!();
    println!("--- Multiple Testing Significance Summary ---");
    println!("Significant OGs (FDR q < 0.05):        {:>6}", summary.significant_fdr_05);
    println!("Significant OGs (FDR q < 0.01):        {:>6}", summary.significant_fdr_01);
    println!("Significant OGs (Bonferroni p < 0.05): {:>6}", summary.significant_bonferroni);
    println!();
    println!("--- Top 10 Genomic Hits (Sorted by FDR q-value) ---");
    println!(
        "{:<12} {:>10} {:>10} {:>8} {:>12} {:>12} {:>8} {:>10}",
        "Orthogroup", "Beta_OG", "Std. Error", "t value", "p-value", "FDR q-val", "R2", "Species"
    );
    println!("{}", "-".repeat(88));

    for hit in &summary.top_hits {
        println!(
            "{:<12} {:>10.4} {:>10.4} {:>8.3} {:>12.4e} {:>12.4e} {:>8.4} {:>7}/{}",
            hit.orthogroup_id,
            hit.beta_og,
            hit.se_og,
            hit.t_stat,
            hit.p_value,
            hit.q_value_fdr,
            hit.r_squared,
            hit.species_with_gene,
            summary.tri_report.matched_species_count
        );
    }
    println!("======================================================================");
}
