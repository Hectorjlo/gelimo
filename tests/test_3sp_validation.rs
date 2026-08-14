mod common;

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
