use gelimo::data::ComparativeDataOptions;
use gelimo::genomics::align_tripartite;

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
}
