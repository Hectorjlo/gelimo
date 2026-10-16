#[test]
fn test_3sp_files_exist() {
    let tree_path = std::path::Path::new("data/benchmark_3sp.tree");
    let tsv_path = std::path::Path::new("data/benchmark_3sp.tsv");

    assert!(tree_path.exists(), "benchmark_3sp.tree file should exist");
    assert!(tsv_path.exists(), "benchmark_3sp.tsv file should exist");
}

#[test]
fn test_3sp_comparative_data_loading() {
    use gelimo::data::ComparativeData;

    let comp_data = ComparativeData::from_files("data/benchmark_3sp.tree", "data/benchmark_3sp.tsv")
        .expect("Failed to load ComparativeData from 3sp benchmark files");

    assert_eq!(comp_data.species_names.len(), 3);
    assert_eq!(
        comp_data.species_names,
        vec!["Canis_lupus_familiaris", "Homo_sapiens", "Macaca_mulatta"]
    );
    assert_eq!(comp_data.data.height(), 3);
}

#[test]
fn test_3sp_vcv_and_pagel_lambda_phase_1_and_2() {
    use approx::assert_relative_eq;
    use gelimo::data::ComparativeData;
    use gelimo::math::{apply_pagel_lambda, compute_vcv};

    let comp_data = ComparativeData::from_files("data/benchmark_3sp.tree", "data/benchmark_3sp.tsv")
        .expect("Failed to load ComparativeData from 3sp benchmark files");

    // Phase 1: Compute VCV Matrix (V)
    let vcv = compute_vcv(&comp_data.tree);

    // Tips order: Canis_lupus_familiaris (idx 0), Homo_sapiens (idx 1), Macaca_mulatta (idx 2)
    // Tree: (Canis_lupus_familiaris_protein:0.5,(Homo_sapiens_protein:1,Macaca_mulatta_protein.faa:1)1:0.5);
    // Root to Canis = 0.5
    // Root to Homo = 0.5 + 1.0 = 1.5
    // Root to Macaca = 0.5 + 1.0 = 1.5
    // MRCA(Homo, Macaca) = Node 1 -> dist = 0.5

    assert_relative_eq!(vcv[(0, 0)], 0.5, epsilon = 1e-12);
    assert_relative_eq!(vcv[(0, 1)], 0.0, epsilon = 1e-12);
    assert_relative_eq!(vcv[(0, 2)], 0.0, epsilon = 1e-12);

    assert_relative_eq!(vcv[(1, 0)], 0.0, epsilon = 1e-12);
    assert_relative_eq!(vcv[(1, 1)], 1.5, epsilon = 1e-12);
    assert_relative_eq!(vcv[(1, 2)], 0.5, epsilon = 1e-12);

    assert_relative_eq!(vcv[(2, 0)], 0.0, epsilon = 1e-12);
    assert_relative_eq!(vcv[(2, 1)], 0.5, epsilon = 1e-12);
    assert_relative_eq!(vcv[(2, 2)], 1.5, epsilon = 1e-12);

    // Phase 2: Apply Pagel's Lambda (lambda = 0.5)
    let v_lambda = apply_pagel_lambda(&vcv, 0.5);

    assert_relative_eq!(v_lambda[(0, 0)], 0.5, epsilon = 1e-12);
    assert_relative_eq!(v_lambda[(1, 1)], 1.5, epsilon = 1e-12);
    assert_relative_eq!(v_lambda[(2, 2)], 1.5, epsilon = 1e-12);

    // Off-diagonal MRCA(Homo, Macaca) rescaled: 0.5 * 0.5 = 0.25
    assert_relative_eq!(v_lambda[(1, 2)], 0.25, epsilon = 1e-12);
    assert_relative_eq!(v_lambda[(2, 1)], 0.25, epsilon = 1e-12);
}

#[test]
fn test_3sp_solver_modes_comparison_phase_4() {
    use approx::assert_relative_eq;
    use gelimo::data::ComparativeData;
    use gelimo::math::{
        compute_vcv, solve_gls_linear, solve_gls_with_mode, GlsSolverMode,
        LinearSolverWorkspace,
    };

    let comp_data = ComparativeData::from_files("data/benchmark_3sp.tree", "data/benchmark_3sp.tsv")
        .expect("Failed to load ComparativeData");

    let y = comp_data.to_response_vector("Maximum longevity (yrs)").unwrap();
    let x = comp_data.to_design_matrix(&["Adult weight (g)"]).unwrap();
    let vcv = compute_vcv(&comp_data.tree);

    // Mode 1: Dense Cholesky O(N^3)
    let res_dense = solve_gls_with_mode(
        GlsSolverMode::DenseCholesky,
        Some(&comp_data.tree),
        Some(&vcv),
        &x,
        &y,
    )
    .expect("Dense Cholesky solve failed");

    // Mode 2: Linear Ho & Ané (2014) O(N) with pre-allocated zero-allocation workspace
    let mut ws = LinearSolverWorkspace::new(comp_data.tree.nodes.len(), 2);
    let res_linear = solve_gls_linear(&comp_data.tree, &x, &y, Some(&mut ws))
        .expect("Linear Ho & Ané solve failed");

    // Floating-point precision validation (tolerance <= 1e-8)
    assert_relative_eq!(
        res_dense.beta[(0, 0)],
        res_linear.beta[(0, 0)],
        epsilon = 1e-8
    );
    assert_relative_eq!(
        res_dense.beta[(1, 0)],
        res_linear.beta[(1, 0)],
        epsilon = 1e-8
    );

    assert_relative_eq!(res_dense.sigma2_ml, res_linear.sigma2_ml, epsilon = 1e-8);
    assert_relative_eq!(
        res_dense.log_likelihood,
        res_linear.log_likelihood,
        epsilon = 1e-8
    );
}
