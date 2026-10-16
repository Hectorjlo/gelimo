# gelimo

[![CI](https://github.com/Hectorjlo/gelimo/actions/workflows/ci.yml/badge.svg)](https://github.com/Hectorjlo/gelimo/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](./LICENSE)
[![Rust](https://img.shields.io/badge/Rust-2024%20Edition-orange.svg)](https://www.rust-lang.org/)

**High-Throughput Generalized Phylogenetic Linear Models in Rust and R**

`gelimo` is an ultra-fast, memory-efficient engine for **Phylogenetic Generalized Least Squares (PGLS)**, designed for comparative genomics and evolutionary trait association. It offers full compatibility with the statistical specifications of R's `caper` package while accelerating whole-genome orthogroup screening by **over $25,000\times$ to $2,000,000\times$** through parallelization, zero-copy matrix factorizations, and pre-computed linear projections.

Available as both a **standalone CLI binary** and a **native R package** powered by `extendr`.

---

## Key Features

- **Ultra-Fast Genome-Wide PGLS**: Fit $>16,000$ PGLS models across 200+ species in **under 0.6 seconds** (compared to hours or days in traditional R loops).
- **OrthoFinder Integration**: Direct streaming ingestion of `Orthogroups.GeneCount.tsv` tables with automatic species matching, QC filtering, and log transformations (`log10(count + 1)`).
- **Multiple Statistical Modes**:
  - **Traditional Trait-as-Response Mode**: $\text{Trait} \sim \text{GFS} + \text{Covariates}$
  - **Inverse Mode**: $\text{GFS} \sim \text{Trait}_1 + \text{Trait}_2 + \dots$ (utilizing single pre-inverted projection matrices).
  - **Partial $r$ Effect Sizes**: Exact partial correlation coefficients controlling for allometric covariates.
- **Continuous Pagel's $\lambda$ Optimization**: Fast 1D Maximum Likelihood (ML) golden-section search.
- **Fuzzy Taxonomy Resolution**: Handles underscore/space mismatches, subspecies strings, and automatic tree/table pruning.
- **Dual Interface**: Full-featured Rust CLI and native R functions returning standard R `data.frame` objects.

---

## Empirical Benchmarks

Evaluated on a mammalian dataset (**216 species**, **16,916 QC-passed OrthoFinder orthogroups**) testing gene family size evolution against Maximum Longevity and Adult Body Weight:

| Implementation / Engine | Methodology | CPU Cores | Time per Model | Total Runtime (16,916 OGs) | Throughput | Speedup vs `gelimo` |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **`gelimo` (Rust / R)** | **Pre-factored Cholesky + Projection** | **Multi-thread** | **33 µs** | **0.57 seconds** | **29,677 mod/s** | **Baseline ($1\times$)** |
| **R Multicore Loop (`nlme::gls` + `mclapply`)** | Fixed $\lambda=1.0$ (Brownian Motion) | 16 cores | **0.863 s** | **~4.05 hours** | 1.16 mod/s | **$25,580\times$ slower** |
| **Traditional R Loop (`caper::pgls`)** | Sequential per-model ML $\lambda$ | 1 core | **74.93 s** | **~352 hours (~14.6 days)** | 0.013 mod/s | **$2,223,500\times$ slower** |

> All implementations yield analytically identical $\beta$ coefficients, $t$-statistics, and $p$-values up to 64-bit floating-point precision (`f64`).

---

## Installation

### 1. Standalone Rust CLI
Requires a working [Rust toolchain](https://rustup.rs/) ($\ge 1.85$):
```bash
git clone https://github.com/Hectorjlo/gelimo.git
cd gelimo
cargo build --release

# The binary will be available at target/release/gelimo
./target/release/gelimo --help
```

### 2. Native R Package
Install directly in R via `devtools` or from local source:
```r
# From GitHub
devtools::install_github("Hectorjlo/gelimo")

# Or from local source directory
devtools::install(".")
```

---

## Quickstart: Command Line (CLI)

### 1. Single Comparative Trait Model (`gelimo fit`)
Fit a single PGLS model with Pagel's $\lambda$ ML optimization:
```bash
gelimo fit \
  --tree species_tree.tree \
  --data traits.tsv \
  --response "Maximum longevity (yrs)" \
  --predictors "Log Adult weight (g)" \
  --opt-lambda \
  --prune
```

### 2. Whole-Genome Orthogroup Screening (`gelimo orthogroups`)
Screen all orthogroups from OrthoFinder with comprehensive QC filters and tree export:
```bash
gelimo orthogroups \
  --tree species_tree.tree \
  --data traits.tsv \
  --orthogroups Orthogroups.GeneCount.tsv \
  --response "Maximum longevity (yrs)" \
  --covariates "Log Adult weight (g)" \
  --output pgls_results.tsv \
  --log-count \
  --gfs-as-response \
  --opt-lambda \
  --min-prevalence 0.80 \
  --min-count-spread 2.0 \
  --save-tree pruned_species_tree.nwk \
  --eval-ranges 0.5,1.0,2.0,3.0 \
  --threads 8 \
  --prune
```

---

## Quickstart: R Package

```r
library(gelimo)

# 1. Genome-Wide Orthogroup Screening (returns native R data.frame in < 0.6s)
results <- pgls_orthogroups(
  tree = "species_tree.tree",
  data = "traits.tsv",
  orthogroups = "Orthogroups.GeneCount.tsv",
  response = "Maximum longevity (yrs)",
  covariates = c("Log Adult weight (g)"),
  log_count = TRUE,
  gfs_as_response = TRUE,       # Inverse model mode: GFS ~ Phenotypes
  opt_lambda = TRUE,            # Optimize Pagel's lambda on baseline trait
  min_prevalence = 0.80,        # Require gene in >= 80% of species
  min_count_spread = 2.0,       # Interdecile count spread P95 - P05 >= 2.0
  save_tree = "pruned_tree.nwk", # Save aligned pruned Newick tree
  output_tsv = "results.tsv",   # Optional TSV export
  prune = TRUE
)

# View top genomic hits sorted by FDR q-value
head(results)
cat("Significant hits (FDR q < 0.05):", sum(results$q_value_fdr < 0.05), "\n")

# 2. Fit a Single Comparative Trait Model
fit <- pgls_fit(
  tree = "species_tree.tree",
  data = "traits.tsv",
  response = "Maximum longevity (yrs)",
  predictors = c("Log Adult weight (g)"),
  opt_lambda = TRUE,
  prune = TRUE
)

print(fit$coefficients)
cat("R-squared:", fit$r_squared, "AIC:", fit$aic, "Lambda:", fit$optimal_lambda, "\n")
```

---

## Quality Control (QC) & Comparative Parameters

| Flag | R Argument | Default | Description |
| :--- | :--- | :--- | :--- |
| `-T, --threads` | — | All cores | Number of CPU threads for Rayon parallel execution |
| `--min-prevalence` | `min_prevalence` | `None` | Minimum prevalence proportion (`0.80`) or percentage (`80%`) of species sharing $\ge 1$ copy |
| `--min-species` | `min_species` | `5` | Absolute minimum number of species possessing at least 1 copy of the gene |
| `--min-count-spread` | `min_count_spread` | `None` | Minimum interdecile count spread ($P_{95} - P_{05}$) across species (robust to outliers) |
| `--min-count-range` | `min_count_range` | `None` | Minimum absolute count range ($\max - \min$) across species |
| `--top-ogs` | `top_ogs` | `None` | Retain only top $N$ most variable orthogroups sorted by count variance |
| `--save-tree` | `save_tree` | `None` | File path to export the pruned phylogenetic tree in Newick format |
| `--eval-ranges` | — | `None` | Comma-separated thresholds (e.g. `0.5,1.0,2.0`) to print a multi-threshold power & FDR summary table |
| `--log-count` | `log_count` | `false` | Apply $\log_{10}(\text{count} + 1)$ transformation to orthogroup gene counts |
| `--binary` | `binary_mode` | `false` | Convert gene counts to binary presence/absence ($0$ vs $1$) |
| `--gfs-as-response`| `gfs_as_response` | `false` | Fit inverse PGLS model ($\text{GFS} \sim \text{Trait} + \text{Covariates}$) |
| `--opt-lambda` | `opt_lambda` | `false` | Optimize Pagel's $\lambda$ via Maximum Likelihood on baseline trait |

---

## Input File Formats & Parameter Specifications

### 1. Phylogenetic Tree (`--tree` / `-t`)
- **Format**: Newick format (e.g. `SpeciesTree_rooted.txt`, `.nwk`, `.tree`).
- **Structure**: Tips must contain the species identifiers (e.g. `Homo_sapiens`, `Mus_musculus`).

### 2. Phenotypic Trait Table (`--data` / `-d`)
- **Format**: TSV (tab-separated) or CSV (comma-separated).
- **Rows**: Each row represents one species.
- **Species Identification**: `gelimo` automatically detects species names formatted as:
  - **Single column**: `species`, `Species`, `species_name`, or `Spp.Name` (e.g. `Homo_sapiens` or `Homo sapiens`).
  - **Two columns**: Separate `Genus` and `Species` columns (e.g. `Homo` and `sapiens`), which are merged and normalized automatically with fuzzy matching.
- **Trait Columns**: Numerical columns containing continuous variables (phenotypic traits and covariates).

*Example `--data` TSV:*
```tsv
Genus	Species	Maximum longevity (yrs)	Adult weight (g)	Log Adult weight (g)
Homo	sapiens	122.5	70000	4.845
Mus	musculus	4.0	20	1.301
Canis	lupus	20.0	35000	4.544
```

### 3. OrthoFinder Gene Count Table (`--orthogroups` / `-g`)
- **Format**: Direct, unmodified `Orthogroups.GeneCount.tsv` file from OrthoFinder.
- **Structure**:
  - First column: `Orthogroup` identifier (e.g. `OG0000001`).
  - Columns 2 to $N$: Species names with integer gene copy counts.
  - Final column: `Total` (detected and ignored automatically).

### 4. Phenotypic Response (`-y`) & Covariates (`-x`)
- **`-y, --response`**: The exact column name in `--data` to be used as the primary phenotype of interest (e.g. `-y "Maximum longevity (yrs)"`).
- **`-x, --covariates`**: *(Optional)* Comma-separated list of column names in `--data` used to control for confounding effects such as body mass or allometry (e.g. `-x "Log Adult weight (g)"`).

### 5. Statistical Model Formulations

`gelimo` supports two complementary modeling frameworks:

- **Mode A: Traditional Trait Response Mode (Default)**
  $$\text{Phenotype } (-y) \sim \text{GFS}_{\text{OG}} + \text{Covariate}_1 (-x) + \dots$$
  Tests whether gene family expansions directly predict changes in the phenotypic trait controlling for covariates.

- **Mode B: Inverse Mode (`--gfs-as-response` / `gfs_as_response = TRUE`)**
  $$\text{GFS}_{\text{OG}} \sim \text{Phenotype } (-y) + \text{Covariate}_1 (-x) + \dots$$
  Tests whether phenotypic divergence drives gene family expansions/contractions across the phylogeny, computing the exact **Partial $r$** effect size for the target trait.

---

## Columns in Orthogroups Results

| Column | Description |
| :--- | :--- |
| `Orthogroup` | OrthoFinder orthogroup identifier (e.g. `OG0010602`) |
| `Beta_OG` | Estimated effect size (slope) for the orthogroup gene count |
| `StdError` | Standard error of the $\beta$ coefficient |
| `t_stat` | Student's $t$-statistic |
| `p_value` | Raw two-tailed $p$-value |
| `q_value_fdr` | False Discovery Rate adjusted $q$-value (Benjamini-Hochberg) |
| `p_bonferroni` | Family-wise error rate adjusted $p$-value (Bonferroni) |
| `Partial_r` | Partial correlation coefficient controlling for covariates |
| `R_squared` | Coefficient of determination ($R^2$) with phylogenetic projection |
| `Sigma2` | Maximum Likelihood evolutionary rate ($\sigma^2$) |
| `Mean_Copies` | Mean gene copy count across aligned species |
| `Species_With_Gene` | Number of aligned species possessing at least 1 copy of the gene |

---

## Architecture & Design Principles

`gelimo` is engineered under the following core principles:
1. **Zero-Copy Matrix Decomposition**: The phylogenetic covariance matrix ($V$) and phenotype design matrices are decomposed only once ($1 \times O(N^3)$), eliminating redundant inversions across 30,000 iterations.
2. **Pre-Computed Linear Projections**: In inverse mode, the projection operator $P = (X^{*T} X^*)^{-1} X^{*T}$ reduces each orthogroup evaluation to a single BLAS matrix-vector product in nanoseconds.
3. **Decoupled Architecture**: `extendr-api` is isolated behind the optional `r-bindings` Cargo feature, keeping standalone CLI builds lightweight and free of R runtime dependencies.

---

## Acknowledgments

- Special thanks to **Dr. Benjamin** for invaluable methodological guidance, statistical insights into high-throughput phylogenetic linear modeling, and the formulation of the inverse model and Type II framework.
- Honorable mention to **Carlos García González** for initial workflows and early exploratory scripts on comparative orthogroup analysis that inspired the development and optimization of this engine.

---

## Citation

If you use `gelimo` in your research, please cite:

```bibtex
@software{lopez_ordaz_2026_gelimo,
  author       = {López-Ordaz, Héctor Jesús},
  title        = {gelimo: High-Throughput Generalized Phylogenetic Linear Models in Rust and R},
  year         = {2026},
  version      = {0.1.0},
  publisher    = {GitHub},
  url          = {https://github.com/Hectorjlo/gelimo}
}
```

A machine-readable citation format is also available in [`CITATION.cff`](./CITATION.cff).

---

## License & Author

- **Author**: Héctor Jesús López-Ordaz (`<hectorjl@lcg.unam.mx>`) — Licenciatura en Ciencias Genómicas, UNAM.
- **License**: [MIT License](./LICENSE).