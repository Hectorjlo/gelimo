//! Command-line interface orchestration and runners for `gelimo`.

pub mod args;
pub mod display;

use clap::Parser;
use colored::Colorize;
use crate::genomics::{run_ortho_pipeline, OrthoPipelineConfig};
use crate::math::GlsSolverMode;
use crate::pipeline::{run_pgls_pipeline, PglsConfig};

pub use args::{Cli, Commands, FitArgs, OrthoArgs};
pub use display::{format_p_value, print_fit_summary, print_ortho_summary};

/// Main CLI entrypoint parser and dispatcher.
pub fn run() {
    let cli = Cli::parse();

    if let Some(threads) = cli.threads
        && threads > 0
    {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();
    }

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
        Ok(summary) => {
            print_fit_summary(&config, &summary);
        }
        Err(err) => {
            eprintln!("{} {}", "Error fitting PGLS model:".red().bold(), err);
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
        min_prevalence: args.min_prevalence,
        min_variance: args.min_variance,
        min_count_spread: args.min_count_spread,
        min_count_range: args.min_count_range,
        top_ogs: args.top_ogs,
        save_tree_path: args.save_tree,
        eval_ranges: args.eval_ranges,
        binary_mode: args.binary,
        log_count_mode: args.log_count,
        gfs_as_response: args.gfs_as_response,
        prune: args.prune,
        fuzzy_taxonomy: true,
    };

    println!("{}", "======================================================================".cyan().bold());
    println!("{}", "gelimo v0.1.0 - High-Throughput Genomic Orthogroups PGLS Engine".bright_white().bold());
    println!("{}", "======================================================================".cyan().bold());
    println!("{:<18} {}", "Tree file:".bright_cyan().bold(), config.tree_path);
    println!("{:<18} {}", "Phenotype file:".bright_cyan().bold(), config.data_path);
    println!("{:<18} {}", "OrthoFinder file:".bright_cyan().bold(), config.ortho_path);
    println!(
        "{:<18} {}",
        "Model Formula:".bright_cyan().bold(),
        if config.gfs_as_response {
            format!("GFS ~ {} (Dr. Benjamin's Inverse Mode)", if config.covariate_cols.is_empty() { config.response_col.clone() } else { format!("{} + {}", config.response_col, config.covariate_cols.join(" + ")) }).bright_yellow().bold()
        } else {
            format!("{} ~ GFS{}", config.response_col, if config.covariate_cols.is_empty() { "".to_string() } else { format!(" + {}", config.covariate_cols.join(" + ")) }).yellow().bold()
        }
    );
    if !config.covariate_cols.is_empty() {
        println!("{:<18} {}", "Covariates:".bright_cyan().bold(), config.covariate_cols.join(", ").yellow());
    }
    println!("{:<18} {}", "Output TSV:".bright_cyan().bold(), config.output_path.green().bold());
    println!(
        "{:<18} {}",
        "Count Mode:".bright_cyan().bold(),
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
