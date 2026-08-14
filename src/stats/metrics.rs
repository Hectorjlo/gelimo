//! Statistical metrics and likelihood evaluation for PGLS.

use faer::Mat;
use crate::math::solver::{solve_gls, GlsResult};

/// Computes Log-Likelihood for a given VCV matrix, design matrix X, and response Y.
pub fn compute_log_likelihood(vcv: &Mat<f64>, x: &Mat<f64>, y: &Mat<f64>) -> Result<f64, String> {
    let res = solve_gls(vcv, x, y)?;
    Ok(res.log_likelihood)
}

/// Computes Akaike Information Criterion (AIC) given Log-Likelihood and number of estimated parameters $k$.
///
/// $\text{AIC} = 2k - 2\ln L$
// ponytail: one-line AIC metric computation
pub fn compute_aic(log_likelihood: f64, num_params: usize) -> f64 {
    2.0 * (num_params as f64) - 2.0 * log_likelihood
}

/// Computes $R^2$ (Phylogenetic Coefficient of Determination) from transformed system residuals.
///
/// $R^2 = 1 - \frac{\text{RSS}^*}{\text{TSS}^*}$
// ponytail: transformed residual coefficient of determination R^2
pub fn compute_r_squared(gls_result: &GlsResult, y_star: &Mat<f64>) -> f64 {
    let n = y_star.nrows();
    let mut sum_y = 0.0;
    for i in 0..n {
        sum_y += y_star[(i, 0)];
    }
    let mean_y = sum_y / (n as f64);

    let mut tss_star = 0.0;
    for i in 0..n {
        let diff = y_star[(i, 0)] - mean_y;
        tss_star += diff * diff;
    }

    let mut rss_star = 0.0;
    for i in 0..n {
        let e = gls_result.transformed_residuals[(i, 0)];
        rss_star += e * e;
    }

    if tss_star.abs() < 1e-15 {
        0.0
    } else {
        (1.0 - rss_star / tss_star).clamp(0.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_compute_aic() {
        // k = 2, logL = -10.0 -> AIC = 4.0 - (-20.0) = 24.0
        assert_relative_eq!(compute_aic(-10.0, 2), 24.0, epsilon = 1e-12);
    }
}
