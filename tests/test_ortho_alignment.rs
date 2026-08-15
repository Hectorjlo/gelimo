use gelimo::data::ComparativeDataOptions;
use gelimo::genomics::{
    align_tripartite, filter_orthogroups, run_ortho_pgls_parallel, OrthoFilterOptions,
    OrthoPglsConfig,
};

#[test]
fn test_tripartite_alignment_mammals_real_data() {
    let tree_path = "user/test/mammals_ng/tree/SpeciesTree_rooted.txt";
    let pheno_path = "user/test/AnAge/anage_data_log_mammals_filtered.tsv";
    let ortho_path = "user/test/mammals_ng/Orthogroups/Orthogroups.GeneCount.tsv";

    if !std::path::Path::new(tree_path).exists()
        || !std::path::Path::new(pheno_path).exists()
        || !std::path::Path::new(ortho_path).exists()
    {
        println!("Skipping real data test: files not found");
        return;
    }

    let options = ComparativeDataOptions {
        prune: true,
        fuzzy_taxonomy: true,
    };

    // Phase A: Tripartite Alignment
    let start = std::time::Instant::now();
    let tri_data = align_tripartite(tree_path, pheno_path, &options, ortho_path)
        .expect("Tripartite alignment failed");
    let elapsed = start.elapsed();

    println!("Tripartite alignment completed in: {:?}", elapsed);
    println!("Matched species count: {}", tri_data.species_names.len());
    println!("Total orthogroups: {}", tri_data.orthogroups.num_orthogroups());

    // Phase B: QC Filtering
    let filter_opts = OrthoFilterOptions {
        min_species_count: 10,
        min_variance: 1e-5,
        binary_mode: false,
        species_total_genes: None,
    };

    let filter_start = std::time::Instant::now();
    let filtered_ortho = filter_orthogroups(&tri_data.orthogroups, &filter_opts);
    let filter_elapsed = filter_start.elapsed();

    println!("QC Filtering completed in: {:?}", filter_elapsed);
    println!("Retained orthogroups: {}", filtered_ortho.num_orthogroups());

    // Phase C: Parallel PGLS Engine
    // Wrap aligned pheno DataFrame into ComparativeData helper to extract Y and covariates
    let comp_helper = gelimo::data::ComparativeData {
        tree: tri_data.tree.clone(),
        species_names: tri_data.species_names.clone(),
        data: tri_data.pheno_data.clone(),
        report: Default::default(),
    };

    let y = comp_helper
        .to_response_vector("Maximum longevity (yrs)")
        .expect("Failed to extract Y");

    let cov = comp_helper
        .to_response_vector("Log Adult weight (g)")
        .ok();

    let pgls_config = OrthoPglsConfig { lambda: 0.8723 };

    let pgls_start = std::time::Instant::now();
    let results = run_ortho_pgls_parallel(
        &tri_data.tree,
        &y,
        cov.as_ref(),
        &filtered_ortho,
        &pgls_config,
    )
    .expect("Parallel PGLS execution failed");
    let pgls_elapsed = pgls_start.elapsed();

    println!(
        "Parallel PGLS across {} orthogroups completed in: {:?}",
        results.len(),
        pgls_elapsed
    );

    assert_eq!(results.len(), filtered_ortho.num_orthogroups());

    // Inspect top hits with lowest p-values
    let mut sorted_results = results.clone();
    sorted_results.sort_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap());

    println!("--- Top 5 Longevity-Associated Orthogroups ---");
    for (rank, hit) in sorted_results.iter().take(5).enumerate() {
        println!(
            "{}. {} | Beta: {:+.4} | SE: {:.4} | t: {:+.3} | p-val: {:.4e} | R2: {:.4} | Species: {}/{}",
            rank + 1,
            hit.orthogroup_id,
            hit.beta_og,
            hit.se_og,
            hit.t_stat,
            hit.p_value,
            hit.r_squared,
            hit.species_with_gene,
            tri_data.species_names.len()
        );
    }
}
