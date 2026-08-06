//! Comparative data module for aligning species traits and trees.

use polars::prelude::*;
use std::fs::File;
use std::io::Read;
use std::path::Path;

use crate::tree::node::normalize_species_name;
use crate::tree::parser::parse_newick;
use crate::tree::PhyloTree;

/// Alignment container for species traits and phylogenetic tree matching `caper::comparative.data`.
#[derive(Debug, Clone)]
pub struct ComparativeData {
    pub tree: PhyloTree,
    pub species_names: Vec<String>,
    pub data: DataFrame,
}

#[derive(Debug, thiserror::Error)]
pub enum ComparativeDataError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Tree parse error: {0}")]
    TreeParse(String),

    #[error("Polars data error: {0}")]
    Polars(#[from] PolarsError),

    #[error("Species column not found in trait table")]
    MissingSpeciesColumn,

    #[error("Species mismatch between tree and trait table. Tree species: {tree_species:?}, Table species: {table_species:?}")]
    SpeciesMismatch {
        tree_species: Vec<String>,
        table_species: Vec<String>,
    },
}

impl ComparativeData {
    /// Constructs `ComparativeData` from Newick string and TSV string input with strict matching.
    pub fn from_newick_and_tsv(
        newick_str: &str,
        tsv_str: &str
    ) -> Result<Self, ComparativeDataError> {
        let tree = parse_newick(
            newick_str
        ).map_err(ComparativeDataError::TreeParse)?;

        let cursor = std::io::Cursor::new(tsv_str);
        let df = CsvReadOptions::default()
            .with_has_header(true)
            .with_parse_options(
                CsvParseOptions::default().with_separator(b'\t')
            )
            .into_reader_with_file_handle(cursor)
            .finish()?;

        Self::align(tree, df)
    }

    /// Reads `ComparativeData` directly from Newick tree file and TSV trait data file.
    pub fn from_files<P1: AsRef<Path>, P2: AsRef<Path>>(
        tree_path: P1,
        tsv_path: P2,
    ) -> Result<Self, ComparativeDataError> {
        let mut tree_file = File::open(tree_path)?;
        let mut newick_str = String::new();
        tree_file.read_to_string(&mut newick_str)?;

        let df = CsvReadOptions::default()
            .with_has_header(true)
            .with_parse_options(
                CsvParseOptions::default().with_separator(b'\t')
            )
            .try_into_reader_with_file_path(
                Some(tsv_path.as_ref().to_path_buf())
            )?
            .finish()?;

        let tree = parse_newick(
            &newick_str
        ).map_err(ComparativeDataError::TreeParse)?;
        Self::align(tree, df)
    }

    /// Aligns tree tip species with table species, ensuring strict match and proper row ordering.
    fn align(
        tree: PhyloTree,
        mut df: DataFrame
    ) -> Result<Self, ComparativeDataError> {
        let tree_normalized_tips = tree.normalized_tip_names();

        // Unified species_id extraction supporting Genus+Species or single species column
        let table_species: Vec<String> = if df.column("Genus").is_ok() && df.column("Species").is_ok() {
            let genus_ca = df.column("Genus")?.str()?;
            let species_ca = df.column("Species")?.str()?;
            (0..df.height())
                .map(|i| {
                    let g = genus_ca.get(i).unwrap_or_default();
                    let s = species_ca.get(i).unwrap_or_default();
                    normalize_species_name(&format!("{}_{}", g, s))
                })
                .collect()
        } else if let Some(col) = ["species", "Species", "SPECIES"].into_iter().find_map(|name| df.column(name).ok()) {
            col.str()?
                .into_no_null_iter()
                .map(normalize_species_name)
                .collect()
        } else {
            return Err(ComparativeDataError::MissingSpeciesColumn);
        };

        // Always ensure unified "species_id" column is added to DataFrame for alignment lookup
        let s = StringChunked::new(
            PlSmallStr::from_static("species_id"),
            table_species.clone(),
        )
        .into_series();
        df.with_column(s)?;

        // Strict species matching check
        let mut sorted_tree_species = tree_normalized_tips.clone();
        let mut sorted_table_species = table_species.clone();
        sorted_tree_species.sort();
        sorted_table_species.sort();

        if sorted_tree_species != sorted_table_species {
            return Err(ComparativeDataError::SpeciesMismatch {
                tree_species: tree_normalized_tips,
                table_species,
            });
        }

        // Reorder DataFrame rows to match exact order of tree tips
        let species_col = df.column("species_id")?.str()?.clone();
        let mut row_indices = Vec::with_capacity(tree_normalized_tips.len());

        for target in &tree_normalized_tips {
            let idx = species_col
                .into_iter()
                .position(
                    |opt_val| opt_val == Some(target.as_str())
                )
                .ok_or_else(|| ComparativeDataError::SpeciesMismatch {
                    tree_species: tree_normalized_tips.clone(),
                    table_species: table_species.clone(),
                })?;
            row_indices.push(idx as u32);
        }

        let idx_ca = UInt32Chunked::new(PlSmallStr::from_static("idx"), row_indices);
        let aligned_df = df.take(&idx_ca)?;

        Ok(Self {
            tree,
            species_names: tree_normalized_tips,
            data: aligned_df,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_comparative_data_alignment_success() {
        let newick = "(Canis_lupus_familiaris_protein:0.5,(Homo_sapiens_protein:1,Macaca_mulatta_protein.faa:1)1:0.5);";
        let tsv = "Genus\tSpecies\tCommon name\tBirth weight (g)\tAdult weight (g)\tMaximum longevity (yrs)\n\
                   Canis\tlupus_familiaris\tDog\t250\t20000\t24\n\
                   Homo\tsapiens\tHuman\t3312.5\t62035\t122.5\n\
                   Macaca\tmulatta\tRhesus monkey\t464\t8235\t40\n";

        let comp_data = ComparativeData::from_newick_and_tsv(newick, tsv).expect("ComparativeData alignment failed");
        assert_eq!(comp_data.species_names.len(), 3);
        assert_eq!(comp_data.data.height(), 3);

        // Verify ordering matches tree tips: Canis_lupus_familiaris, Homo_sapiens, Macaca_mulatta
        let species_col = comp_data.data.column("species_id").unwrap().str().unwrap();
        let species_vec: Vec<&str> = species_col.into_iter().map(|s| s.unwrap()).collect();
        assert_eq!(species_vec, vec!["Canis_lupus_familiaris", "Homo_sapiens", "Macaca_mulatta"]);
    }

    #[test]
    fn test_comparative_data_mismatch_error() {
        let newick = "((Homo_sapiens_protein:1,Macaca_mulatta_protein.faa:1):0.5);";
        let tsv = "Genus\tSpecies\tCommon name\tBirth weight (g)\n\
                   Homo\tsapiens\tHuman\t3312.5\n";

        let res = ComparativeData::from_newick_and_tsv(newick, tsv);
        assert!(res.is_err(), "Expected error on species mismatch");
    }

    #[test]
    fn test_comparative_data_single_species_column() {
        let newick = "((Homo_sapiens_protein:1,Macaca_mulatta_protein.faa:1):0.5);";
        let tsv = "species\ttrait_val\n\
                   Homo_sapiens\t3312.5\n\
                   Macaca_mulatta\t464.0\n";

        let comp_data = ComparativeData::from_newick_and_tsv(newick, tsv).expect("Single species column alignment failed");
        assert_eq!(comp_data.species_names.len(), 2);
        assert_eq!(comp_data.data.height(), 2);
        assert!(comp_data.data.column("species_id").is_ok());
    }
}
