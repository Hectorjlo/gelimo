use gelimo::genomics::{
    filter_orthogroups, run_ortho_pipeline, OrthoFilterOptions, OrthoPipelineConfig,
};
use gelimo::genomics::ortho_reader::OrthoGeneCountTable;
use gelimo::tree::parse_newick;

#[test]
fn test_min_prevalence_fraction_vs_percentage() {
    let table = OrthoGeneCountTable {
        orthogroup_ids: vec!["OG1".into(), "OG2".into(), "OG3".into()],
        species_names: vec![
            "Sp1".into(), "Sp2".into(), "Sp3".into(), "Sp4".into(), "Sp5".into(),
            "Sp6".into(), "Sp7".into(), "Sp8".into(), "Sp9".into(), "Sp10".into(),
        ],
        // OG1: in 10/10 species (100%)
        // OG2: in 8/10 species (80%)
        // OG3: in 5/10 species (50%)
        counts: vec![
            1.0, 2.0, 1.0, 3.0, 1.0, 2.0, 1.0, 1.0, 2.0, 1.0, // OG1
            1.0, 0.0, 1.0, 2.0, 0.0, 1.0, 3.0, 1.0, 1.0, 2.0, // OG2
            1.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2.0, // OG3
        ],
    };

    // Test fraction 0.80 (requires >= 8 species)
    let opts_frac = OrthoFilterOptions {
        min_prevalence: Some(0.80),
        min_variance: 1e-5,
        ..Default::default()
    };
    let res_frac = filter_orthogroups(&table, &opts_frac);
    assert_eq!(res_frac.num_orthogroups(), 2);
    assert_eq!(res_frac.orthogroup_ids, vec!["OG1", "OG2"]);

    // Test percentage 80.0 (requires >= 8 species) -> same result
    let opts_pct = OrthoFilterOptions {
        min_prevalence: Some(80.0),
        min_variance: 1e-5,
        ..Default::default()
    };
    let res_pct = filter_orthogroups(&table, &opts_pct);
    assert_eq!(res_pct.num_orthogroups(), 2);
    assert_eq!(res_pct.orthogroup_ids, vec!["OG1", "OG2"]);

    // Test fraction 0.90 (requires >= 9 species) -> only OG1
    let opts_strict = OrthoFilterOptions {
        min_prevalence: Some(0.90),
        min_variance: 1e-5,
        ..Default::default()
    };
    let res_strict = filter_orthogroups(&table, &opts_strict);
    assert_eq!(res_strict.num_orthogroups(), 1);
    assert_eq!(res_strict.orthogroup_ids, vec!["OG1"]);
}

#[test]
fn test_min_count_spread_interdecile() {
    let table = OrthoGeneCountTable {
        orthogroup_ids: vec!["OG_narrow".into(), "OG_wide".into()],
        species_names: (1..=20).map(|i| format!("Sp{}", i)).collect(),
        // 20 species: P05 is index 1, P95 is index 18
        // OG_narrow: counts from 1 to 2 -> spread ~ 1
        // OG_wide: counts from 1 to 10 -> spread ~ 8
        counts: {
            let mut c = Vec::new();
            for i in 0..20 {
                c.push(if i % 2 == 0 { 1.0 } else { 2.0 });
            }
            for i in 0..20 {
                c.push(i as f64);
            }
            c
        },
    };

    let opts = OrthoFilterOptions {
        min_count_spread: Some(5.0),
        min_variance: 1e-5,
        ..Default::default()
    };
    let res = filter_orthogroups(&table, &opts);
    assert_eq!(res.num_orthogroups(), 1);
    assert_eq!(res.orthogroup_ids[0], "OG_wide");
    assert_eq!(res.filter_report.discarded_low_spread, 1);
}

#[test]
fn test_end_to_end_pipeline_with_save_tree_and_eval_ranges() {
    let tree_path = "user/test/mammals_ng/tree/SpeciesTree_rooted.txt";
    let pheno_path = "user/test/AnAge/anage_data_log_mammals_filtered.tsv";
    let ortho_path = "user/test/mammals_ng/Orthogroups/Orthogroups.GeneCount.tsv";
    let temp_tree_path = std::env::temp_dir().join("gelimo_output_pruned_tree.nwk");
    let temp_tree_out = temp_tree_path.to_str().unwrap();

    if !std::path::Path::new(tree_path).exists()
        || !std::path::Path::new(pheno_path).exists()
        || !std::path::Path::new(ortho_path).exists()
    {
        println!("Skipping real data pipeline test: files not found");
        return;
    }

    let config = OrthoPipelineConfig {
        tree_path: tree_path.to_string(),
        data_path: pheno_path.to_string(),
        ortho_path: ortho_path.to_string(),
        response_col: "Maximum longevity (yrs)".to_string(),
        covariate_cols: vec!["Log Adult weight (g)".to_string()],
        output_path: "".to_string(), // In-memory only
        lambda: 1.0,
        optimize_lambda: false,
        min_species_count: 5,
        min_prevalence: Some(0.80),
        min_variance: 1e-6,
        min_count_spread: None,
        min_count_range: None,
        top_ogs: Some(50),
        save_tree_path: Some(temp_tree_out.to_string()),
        eval_ranges: Some(vec![1.0, 2.0]),
        binary_mode: false,
        log_count_mode: true,
        gfs_as_response: true,
        prune: true,
        fuzzy_taxonomy: true,
    };

    let summary = run_ortho_pipeline(&config).expect("Pipeline execution failed");

    // Check top_ogs constraint
    assert_eq!(summary.total_analyzed, 50);
    assert_eq!(summary.all_results.len(), 50);

    // Check saved tree
    assert!(temp_tree_path.exists());
    let saved_content = std::fs::read_to_string(&temp_tree_path).expect("Read saved tree");
    let parsed_tree = parse_newick(&saved_content).expect("Parse saved tree");
    assert_eq!(parsed_tree.tips.len(), summary.tri_report.matched_species_count);

    // Check multi-threshold summary
    assert_eq!(summary.eval_ranges_summary.len(), 2);
    assert_eq!(summary.eval_ranges_summary[0].threshold_spread, 1.0);
    assert_eq!(summary.eval_ranges_summary[1].threshold_spread, 2.0);

    // Clean up temporary tree file
    let _ = std::fs::remove_file(&temp_tree_path);
}
