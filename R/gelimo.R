#' High-Throughput Genomic PGLS on OrthoFinder Orthogroups
#'
#' Runs ultra-fast parallel PGLS screening across tens of thousands of orthogroups
#' from OrthoFinder gene count tables using the high-performance Rust engine.
#'
#' @param tree Character. Path to Newick phylogenetic tree file.
#' @param data Character. Path to CSV/TSV species trait data file.
#' @param orthogroups Character. Path to OrthoFinder Orthogroups.GeneCount.tsv file.
#' @param response Character. Name of response variable column in trait data (Y).
#' @param covariates Character vector. Optional covariate column names (e.g. "Log Adult weight (g)").
#' @param output_tsv Character. Optional path to save full results as TSV file.
#' @param log_count Logical. If TRUE, applies log10(count + 1) transformation (default: FALSE).
#' @param gfs_as_response Logical. If TRUE, fits Dr. Benjamin's inverse model (GFS ~ Phenotypes). Default: FALSE.
#' @param lambda Numeric. Pagel's lambda parameter (default: 1.0).
#' @param opt_lambda Logical. If TRUE, optimizes Pagel's lambda by Maximum Likelihood on baseline trait model.
#' @param min_species Integer. Minimum species count possessing at least 1 gene copy (default: 5).
#' @param min_prevalence Numeric. Minimum species prevalence proportion (e.g. 0.80) or percentage (80).
#' @param min_variance Numeric. Minimum variance threshold for QC filtering (default: 1e-6).
#' @param min_count_spread Numeric. Minimum interdecile count spread (P95 - P05) across species.
#' @param min_count_range Numeric. Minimum absolute count range (max - min) across species.
#' @param top_ogs Integer. Maximum number of top most variable orthogroups to retain.
#' @param save_tree Character. Optional path to save the pruned phylogenetic tree as Newick file.
#' @param prune Logical. If TRUE, automatically prunes non-matching tree tips and table rows (default: TRUE).
#'
#' @return A data.frame sorted by FDR q-value containing PGLS coefficients, p-values, FDR q-values, Bonferroni corrections, and Partial r effect sizes.
#' @export
#'
#' @examples
#' \dontrun{
#' results <- pgls_orthogroups(
#'   tree = "tree.txt",
#'   data = "phenotypes.tsv",
#'   orthogroups = "Orthogroups.GeneCount.tsv",
#'   response = "Maximum longevity (yrs)",
#'   covariates = "Log Adult weight (g)",
#'   log_count = TRUE,
#'   opt_lambda = TRUE,
#'   min_prevalence = 0.80,
#'   save_tree = "pruned_tree.nwk"
#' )
#' head(results)
#' }
pgls_orthogroups <- function(
    tree,
    data,
    orthogroups,
    response,
    covariates = NULL,
    output_tsv = NULL,
    log_count = FALSE,
    gfs_as_response = FALSE,
    lambda = 1.0,
    opt_lambda = FALSE,
    min_species = 5L,
    min_prevalence = NULL,
    min_variance = 1e-6,
    min_count_spread = NULL,
    min_count_range = NULL,
    top_ogs = NULL,
    save_tree = NULL,
    prune = TRUE) {
  
  if (!file.exists(tree)) stop("Phylogenetic tree file not found: ", tree, call. = FALSE)
  if (!file.exists(data)) stop("Trait data file not found: ", data, call. = FALSE)
  if (!file.exists(orthogroups)) stop("OrthoFinder file not found: ", orthogroups, call. = FALSE)
  
  covs <- if (!is.null(covariates) && length(covariates) > 0) as.character(covariates) else NULL
  
  r_pgls_orthogroups(
    tree = as.character(tree),
    data = as.character(data),
    orthogroups = as.character(orthogroups),
    response = as.character(response),
    covariates = covs,
    output_tsv = if (!is.null(output_tsv)) as.character(output_tsv) else NULL,
    log_count = as.logical(log_count),
    gfs_as_response = as.logical(gfs_as_response),
    lambda = as.numeric(lambda),
    opt_lambda = as.logical(opt_lambda),
    min_species = as.integer(min_species),
    min_prevalence = if (!is.null(min_prevalence)) as.numeric(min_prevalence) else NULL,
    min_variance = as.numeric(min_variance),
    min_count_spread = if (!is.null(min_count_spread)) as.numeric(min_count_spread) else NULL,
    min_count_range = if (!is.null(min_count_range)) as.numeric(min_count_range) else NULL,
    top_ogs = if (!is.null(top_ogs)) as.integer(top_ogs) else NULL,
    save_tree = if (!is.null(save_tree)) as.character(save_tree) else NULL,
    prune = as.logical(prune)
  )
}

#' Fit Single Generalized Phylogenetic Linear Model (PGLS)
#'
#' Fits a comparative PGLS model on species traits with Pagel's lambda optimization.
#'
#' @param tree Character. Path to Newick phylogenetic tree file.
#' @param data Character. Path to CSV/TSV species trait data file.
#' @param response Character. Name of response variable column (Y).
#' @param predictors Character vector. Predictor variable column names (X).
#' @param lambda Numeric. Pagel's lambda parameter (default: 1.0).
#' @param opt_lambda Logical. If TRUE, optimizes lambda by Maximum Likelihood (ML).
#' @param prune Logical. If TRUE, automatically prunes non-matching species (default: TRUE).
#' @param solver Character. Solver mode: "dense" (Cholesky O(N^3)) or "linear" (Ho & Ane O(N)).
#'
#' @return A list containing coefficients data.frame, r_squared, aic, log_likelihood, sigma2, lambda, and matched_species.
#' @export
pgls_fit <- function(
    tree,
    data,
    response,
    predictors = NULL,
    lambda = 1.0,
    opt_lambda = FALSE,
    prune = TRUE,
    solver = "dense") {
  
  if (!file.exists(tree)) stop("Phylogenetic tree file not found: ", tree, call. = FALSE)
  if (!file.exists(data)) stop("Trait data file not found: ", data, call. = FALSE)
  
  preds <- if (!is.null(predictors) && length(predictors) > 0) as.character(predictors) else NULL
  
  r_pgls_fit(
    tree = as.character(tree),
    data = as.character(data),
    response = as.character(response),
    predictors = preds,
    lambda = as.numeric(lambda),
    opt_lambda = as.logical(opt_lambda),
    prune = as.logical(prune),
    solver = as.character(solver)
  )
}
