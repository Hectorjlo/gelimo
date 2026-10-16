//! Maximum Likelihood (ML) evolutionary parameter optimization and parallel grid search.

pub mod lambda;

pub use lambda::{optimize_pagel_lambda, LambdaOptResult};
