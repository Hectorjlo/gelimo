use clap::Parser;
use gelimo::math::GlsSolverMode;
use gelimo::pipeline::{run_pgls_pipeline, PglsConfig, PglsSummary};

/// gelimo: Generalized Phylogenetic Linear Models Engine in Rust
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
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

fn main() {
    let args = Args::parse();

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
        Ok(summary) => print_summary(&config, &summary),
        Err(err) => {
            eprintln!("Error running PGLS model: {}", err);
            std::process::exit(1);
        }
    }
}

fn print_summary(config: &PglsConfig, summary: &PglsSummary) {
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
