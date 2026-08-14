//! Statistical inference module for PGLS metrics.

pub mod metrics;

pub use metrics::{compute_aic, compute_log_likelihood, compute_r_squared};
