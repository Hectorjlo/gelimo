use clap::Parser;
use gelimo::data::{ComparativeData, ComparativeDataOptions};
use gelimo::tree::parser::parse_newick;
use std::fs;

/// phylo_check: Inspection tool for phylogenetic Newick trees and species data tables
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to Newick phylogenetic tree file
    #[arg(short, long)]
    tree: Option<String>,

    /// Path to CSV/TSV species trait data file
    #[arg(short, long)]
    data: Option<String>,
}

fn main() {
    let args = Args::parse();

    if args.tree.is_none() && args.data.is_none() {
        eprintln!("Please specify at least --tree <FILE> or --data <FILE>");
        std::process::exit(1);
    }

    if let Some(tree_path) = &args.tree {
        println!("=== Tree Inspection: {} ===", tree_path);
        match fs::read_to_string(tree_path) {
            Ok(content) => match parse_newick(&content) {
                Ok(tree) => {
                    println!("Total nodes:    {}", tree.nodes.len());
                    println!("Tips count:     {}", tree.tips.len());
                    println!("Root ID:        {}", tree.root);
                    println!("Tip labels:");
                    for (i, &tip_idx) in tree.tips.iter().enumerate() {
                        let name = tree.nodes[tip_idx].name.as_deref().unwrap_or("<unnamed>");
                        println!("  {:3}. {}", i + 1, name);
                    }
                }
                Err(err) => eprintln!("Error parsing Newick tree: {}", err),
            },
            Err(err) => eprintln!("Error reading file: {}", err),
        }
        println!();
    }

    if let (Some(tree_path), Some(data_path)) = (&args.tree, &args.data) {
        println!("=== Comparative Alignment Test ===");
        let options = ComparativeDataOptions {
            prune: true,
            fuzzy_taxonomy: true,
        };
        match ComparativeData::from_files_with_options(tree_path, data_path, &options) {
            Ok(comp) => {
                println!("Matched species count: {}", comp.report.matched_species_count);
                println!("Data table height:     {}", comp.data.height());
                println!("Data table columns:    {:?}", comp.data.get_column_names());
            }
            Err(err) => eprintln!("Alignment error: {}", err),
        }
    }
}
