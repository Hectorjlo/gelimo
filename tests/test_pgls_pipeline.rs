use gelimo::data::ComparativeData;
use gelimo::math::{apply_pagel_lambda, compute_vcv, solve_gls};
use gelimo::stats::{compute_aic, compute_r_squared};

#[test]
fn test_full_pgls_pipeline_3sp() {
    // 1. Ingestion and ComparativeData alignment
    let comp_data = ComparativeData::from_files("data/benchmark_3sp.tree", "data/benchmark_3sp.tsv")
        .expect("Failed to load ComparativeData");

    assert_eq!(comp_data.species_names.len(), 3);

    // 2. Extract Response Y (Longevity) and Design Matrix X (Intercept + Adult weight)
    let y = comp_data
        .to_response_vector("Maximum longevity (yrs)")
        .expect("Failed to extract Y");
    let x = comp_data
        .to_design_matrix(&["Adult weight (g)"])
        .expect("Failed to extract X");

    assert_eq!(y.nrows(), 3);
    assert_eq!(x.nrows(), 3);
    assert_eq!(x.ncols(), 2);

    // 3. Compute VCV matrix
    let vcv = compute_vcv(&comp_data.tree);

    // 4. Apply Pagel's Lambda transformation (lambda = 1.0 -> standard Brownian Motion)
    let vcv_lambda1 = apply_pagel_lambda(&vcv, 1.0);

    // 5. Solve GLS / PGLS system
    let gls_res = solve_gls(&vcv_lambda1, &x, &y).expect("PGLS solve failed");

    // Check coefficient dimensions and values
    assert_eq!(gls_res.beta.nrows(), 2);
    assert_eq!(gls_res.se.len(), 2);
    assert_eq!(gls_res.residuals.nrows(), 3);

    // Residual variance must be positive
    assert!(gls_res.sigma2_ml > 0.0);
    assert!(gls_res.sigma2_unbiased > gls_res.sigma2_ml);

    // Log-likelihood must be finite
    assert!(gls_res.log_likelihood.is_finite());

    // 6. Compute R^2 and AIC metrics
    let l_mat = vcv_lambda1.as_ref().llt(faer::Side::Lower).unwrap().L().to_owned();
    let y_star = gelimo::math::solver::forward_substitute(&l_mat, &y);
    let x_star = gelimo::math::solver::forward_substitute(&l_mat, &x);

    let r2 = compute_r_squared(&gls_res, &x_star, &y_star);
    let aic = compute_aic(gls_res.log_likelihood, 3); // 2 betas + 1 sigma2

    assert!((0.0..=1.0).contains(&r2));
    assert!(aic.is_finite());
}
