use gelimo::data::ComparativeDataOptions;
use gelimo::genomics::{align_tripartite, filter_orthogroups, OrthoFilterOptions};

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

    let start = std::time::Instant::now();
    let tri_data = align_tripartite(tree_path, pheno_path, &options, ortho_path)
        .expect("Tripartite alignment failed");
    let elapsed = start.elapsed();

    println!("Tripartite alignment completed in: {:?}", elapsed);
    println!("Matched species count: {}", tri_data.species_names.len());
    println!("Total orthogroups: {}", tri_data.orthogroups.num_orthogroups());
    println!("Taxonomy matches: {}", tri_data.report.taxonomy_matches.len());

    assert!(tri_data.species_names.len() > 100);
    assert_eq!(tri_data.tree.tips.len(), tri_data.species_names.len());
    assert_eq!(tri_data.pheno_data.height(), tri_data.species_names.len());
    assert_eq!(tri_data.orthogroups.num_species(), tri_data.species_names.len());
    assert_eq!(tri_data.orthogroups.num_orthogroups(), 29965);

    // Verify first orthogroup slice has exact length N
    let og0_counts = tri_data.orthogroups.get_og_counts(0);
    assert_eq!(og0_counts.len(), tri_data.species_names.len());

    // Phase B Test: QC Filtering
    let filter_opts = OrthoFilterOptions {
        min_species_count: 10, // at least 10 species possess the gene
        min_variance: 1e-5,
        binary_mode: false,
        species_total_genes: None,
    };

    let filter_start = std::time::Instant::now();
    let filtered_ortho = filter_orthogroups(&tri_data.orthogroups, &filter_opts);
    let filter_elapsed = filter_start.elapsed();

    println!("QC Filtering completed in: {:?}", filter_elapsed);
    println!("Initial orthogroups:       {}", filtered_ortho.filter_report.initial_orthogroups);
    println!("Retained orthogroups:      {}", filtered_ortho.filter_report.retained_orthogroups);
    println!("Discarded invariant:       {}", filtered_ortho.filter_report.discarded_invariant);
    println!("Discarded low prevalence:  {}", filtered_ortho.filter_report.discarded_low_prevalence);
    println!("Discarded low variance:    {}", filtered_ortho.filter_report.discarded_low_variance);

    assert!(filtered_ortho.num_orthogroups() > 10000);
    assert_eq!(filtered_ortho.num_species(), tri_data.species_names.len());
}
