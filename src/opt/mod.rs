//! Maximum Likelihood (ML) evolutionary parameter optimization and parallel grid search.

pub mod lambda;

pub use lambda::{
    grid_search_lambda_parallel, optimize_pagel_lambda, LambdaOptResult,
};
