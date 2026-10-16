//! Statistical metrics and likelihood evaluation for PGLS.

use faer::Mat;
use crate::math::solver::GlsResult;

/// Computes Akaike Information Criterion (AIC) given Log-Likelihood and number of estimated parameters $k$.
///
/// $\text{AIC} = 2k - 2\ln L$
pub fn compute_aic(log_likelihood: f64, num_params: usize) -> f64 {
    2.0 * (num_params as f64) - 2.0 * log_likelihood
}

/// Computes phylogenetic $R^2$ (Coefficient of Determination) using exact projection on transformed intercept $X_0^* = L^{-1} \mathbf{1}$.
///
/// $\text{TSS}^* = (Y^* - X_0^* \hat{\mu}_{\text{null}})^T (Y^* - X_0^* \hat{\mu}_{\text{null}})$
/// where $\hat{\mu}_{\text{null}} = \frac{(X_0^*)^T Y^*}{(X_0^*)^T X_0^*}$.
/// $R^2 = 1 - \frac{\text{RSS}^*}{\text{TSS}^*}$.
pub fn compute_r_squared(gls_result: &GlsResult, x0_star: &Mat<f64>, y_star: &Mat<f64>) -> f64 {
    let n = y_star.nrows();
    let mut dot_x0_y = 0.0;
    let mut dot_x0_x0 = 0.0;

    for i in 0..n {
        let x0 = x0_star[(i, 0)];
        let y = y_star[(i, 0)];
        dot_x0_y += x0 * y;
        dot_x0_x0 += x0 * x0;
    }

    let mu_null = if dot_x0_x0 > 1e-15 {
        dot_x0_y / dot_x0_x0
    } else {
        0.0
    };

    let tss_star: f64 = (0..n)
        .map(|i| {
            let diff = y_star[(i, 0)] - x0_star[(i, 0)] * mu_null;
            diff * diff
        })
        .sum();

    let rss_star = gls_result.sigma2_ml * (n as f64);

    if tss_star.abs() < 1e-15 {
        0.0
    } else {
        (1.0 - rss_star / tss_star).clamp(0.0, 1.0)
    }
}

/// Standard normal cumulative distribution function $\Phi(x)$ approximation.
pub fn approx_normal_cdf(x: f64) -> f64 {
    0.5 * (1.0 + erff(x / std::f64::consts::SQRT_2))
}

fn erff(x: f64) -> f64 {
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;
    let p = 0.3275911;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x_abs = x.abs();
    let t = 1.0 / (1.0 + p * x_abs);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x_abs * x_abs).exp();

    sign * y
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
