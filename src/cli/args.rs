//! Command-line argument definitions and styling for `gelimo`.

use clap::builder::styling::{AnsiColor, Effects, Styles};
use clap::{Args, Parser, Subcommand};

pub fn cli_styles() -> Styles {
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
    author = "Héctor Jesús López-Ordaz <hectorjl@lcg.unam.mx>",
    version = "0.1.0",
    about = "Generalized Phylogenetic Linear Models and High-Throughput Genomic Screening Engine in Rust",
    long_about = "High-performance PGLS engine compatible with R caper and optimized for genome-wide orthogroup analysis with Rayon.",
    styles = cli_styles()
)]
pub struct Cli {
    /// Number of parallel CPU threads (default: all available logical cores)
    #[arg(short = 'T', long, global = true)]
    pub threads: Option<usize>,

    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Fit a single PGLS model on comparative species traits (caper compatible)
    #[command(alias = "pgls", styles = cli_styles())]
    Fit(FitArgs),

    /// Run high-throughput parallel PGLS screening on OrthoFinder GeneCount tables
    #[command(alias = "gwas", styles = cli_styles())]
    Orthogroups(OrthoArgs),
}

#[derive(Args, Debug)]
pub struct FitArgs {
    /// Path to Newick phylogenetic tree file
    #[arg(short, long)]
    pub tree: String,

    /// Path to CSV/TSV species trait data file
    #[arg(short, long)]
    pub data: String,

    /// Response variable column name (Y)
    #[arg(short = 'y', long)]
    pub response: String,

    /// Predictor variable column names (X), comma-separated
    #[arg(short = 'x', long, value_delimiter = ',')]
    pub predictors: Vec<String>,

    /// Pagel's lambda transformation parameter (default: 1.0)
    #[arg(short = 'l', long, default_value_t = 1.0)]
    pub lambda: f64,

    /// Optimize Pagel's lambda by Maximum Likelihood (ML)
    #[arg(long, default_value_t = false)]
    pub opt_lambda: bool,

    /// Enable automatic tree and table pruning for non-matching species
    #[arg(short = 'p', long, default_value_t = false)]
    pub prune: bool,

    /// Solver algorithm mode: "dense" (O(N^3) Cholesky) or "linear" (O(N) Ho & Ane 2014)
    #[arg(short = 's', long, default_value = "dense")]
    pub solver: String,
}

#[derive(Args, Debug)]
pub struct OrthoArgs {
    /// Path to Newick phylogenetic tree file
    #[arg(short, long)]
    pub tree: String,

    /// Path to CSV/TSV species trait data file (phenotype + covariates)
    #[arg(short, long)]
    pub data: String,

    /// Path to OrthoFinder Orthogroups.GeneCount.tsv file
    #[arg(short = 'g', long)]
    pub orthogroups: String,

    /// Response variable column name in phenotype table (Y)
    #[arg(short = 'y', long)]
    pub response: String,

    /// Optional covariate column names (e.g. Log Adult weight), comma-separated
    #[arg(short = 'x', long, value_delimiter = ',')]
    pub covariates: Vec<String>,

    /// Output TSV file path for full genome-wide PGLS results
    #[arg(short = 'o', long, default_value = "orthogroups_pgls_results.tsv")]
    pub output: String,

    /// Pagel's lambda transformation parameter (default: 1.0)
    #[arg(short = 'l', long, default_value_t = 1.0)]
    pub lambda: f64,

    /// Optimize Pagel's lambda by Maximum Likelihood (ML) on the baseline trait model
    #[arg(long, default_value_t = false)]
    pub opt_lambda: bool,

    /// Minimum species count possessing at least 1 gene copy (QC filter)
    #[arg(long, default_value_t = 5)]
    pub min_species: usize,

    /// Minimum variance of gene counts across species (QC filter)
    #[arg(long, default_value_t = 1e-6)]
    pub min_variance: f64,

    /// Minimum species prevalence proportion (e.g. 0.8) or percentage (e.g. 80)
    #[arg(long)]
    pub min_prevalence: Option<f64>,

    /// Minimum interdecile count spread (P95 - P05) across species
    #[arg(long)]
    pub min_count_spread: Option<f64>,

    /// Minimum absolute count range (max - min) across species
    #[arg(long)]
    pub min_count_range: Option<f64>,

    /// Maximum number of top most variable orthogroups to retain
    #[arg(long)]
    pub top_ogs: Option<usize>,

    /// Save the pruned phylogenetic tree containing only matched species to Newick file
    #[arg(long)]
    pub save_tree: Option<String>,

    /// Evaluate multiple count spread thresholds (e.g. "3,5,10") and report power comparison
    #[arg(long, value_delimiter = ',')]
    pub eval_ranges: Option<Vec<f64>>,

    /// Convert gene counts to binary presence/absence (0 vs 1)
    #[arg(long, default_value_t = false)]
    pub binary: bool,

    /// Apply log10(count + 1) transformation to gene counts
    #[arg(long, default_value_t = false)]
    pub log_count: bool,

    /// Dr. Benjamin's inverse mode: GFS as response variable (GFS ~ Phenotypes)
    #[arg(long, default_value_t = false)]
    pub gfs_as_response: bool,

    /// Enable automatic tree and table pruning for non-matching species
    #[arg(short = 'p', long, default_value_t = true)]
    pub prune: bool,
}
