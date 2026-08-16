//! Linear algebra, VCV mathematical transformations, and PGLS solvers module.

pub mod linear_solver;
pub mod pagel;
pub mod solver;
pub mod vcv;

pub use linear_solver::{solve_gls_linear, LinearSolverWorkspace};
pub use pagel::{apply_pagel_kappa, apply_pagel_lambda, apply_pagel_lambda_tree};
pub use solver::{solve_gls, solve_gls_with_mode, GlsResult, GlsSolverMode};
pub use vcv::compute_vcv;
