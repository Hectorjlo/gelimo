//! Comparative data module for aligning species traits and trees with pruning and taxonomy reports.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use crate::tree::node::normalize_species_name;
use crate::tree::parser::parse_newick;
use crate::tree::PhyloTree;

/// Lightweight in-memory trait table supporting TSV and CSV formats.
#[derive(Debug, Clone, Default)]
pub struct TraitTable {
    pub headers: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl TraitTable {
    /// Parses a TSV or CSV trait table from a string or reader.
    pub fn from_reader<R: BufRead>(mut reader: R) -> Result<Self, ComparativeDataError> {
        let mut header_line = String::new();
        if reader.read_line(&mut header_line)? == 0 {
            return Err(ComparativeDataError::MissingSpeciesColumn);
        }
        let trimmed_header = header_line.trim_end_matches(&['\r', '\n'][..]);
        let sep = if trimmed_header.contains('\t') { '\t' } else { ',' };
        let headers: Vec<String> = trimmed_header
            .split(sep)
            .map(|s| s.trim().trim_matches('"').to_string())
            .collect();

        let mut rows = Vec::new();
        let mut line = String::new();
        while reader.read_line(&mut line)? > 0 {
            let trimmed = line.trim_end_matches(&['\r', '\n'][..]);
            if !trimmed.trim().is_empty() {
                let row: Vec<String> = trimmed
                    .split(sep)
                    .map(|s| s.trim().trim_matches('"').to_string())
                    .collect();
                rows.push(row);
            }
            line.clear();
        }

        Ok(Self { headers, rows })
    }

    /// Returns the number of rows (species) in the table.
    pub fn height(&self) -> usize {
        self.rows.len()
    }

    /// Finds the index of a column by name.
    pub fn column_index(&self, name: &str) -> Option<usize> {
        self.headers.iter().position(|h| h == name)
    }

    /// Returns a new `TraitTable` with rows reordered according to `indices`.
    pub fn reorder(&self, indices: &[usize]) -> Self {
        let rows = indices
            .iter()
            .filter_map(|&idx| self.rows.get(idx).cloned())
            .collect();
        Self {
            headers: self.headers.clone(),
            rows,
        }
    }

    /// Extracts a column of `f64` values by column name.
    pub fn column_f64(&self, name: &str) -> Result<Vec<f64>, ComparativeDataError> {
        let col_idx = self
            .column_index(name)
            .ok_or_else(|| ComparativeDataError::ColumnError(format!("Column '{}' not found", name)))?;
        Ok(self
            .rows
            .iter()
            .map(|r| {
                r.get(col_idx)
                    .and_then(|s| s.trim().parse::<f64>().ok())
                    .unwrap_or(0.0)
            })
            .collect())
    }
}

/// Options for configuring `ComparativeData` alignment, pruning, and taxonomy matching.
#[derive(Debug, Clone)]
pub struct ComparativeDataOptions {
    /// If `true`, prunes non-matching species between tree and table. Defaults to `false` (returns error on mismatch).
    pub prune: bool,
    /// If `true`, attempts Genus + Subspecies matching (e.g. `Canis_lupus_familiaris` -> `Canis_familiaris`). Defaults to `true`.
    pub fuzzy_taxonomy: bool,
}

impl Default for ComparativeDataOptions {
    fn default() -> Self {
        Self {
            prune: false,
            fuzzy_taxonomy: true,
        }
    }
}

/// Report providing full transparency on matched species count, pruned species, and taxonomy adjustments.
#[derive(Debug, Clone, Default)]
pub struct PruneReport {
    /// Total number of species that successfully matched and were included in analysis.
    pub matched_species_count: usize,
    /// Tips removed from tree because they were missing in the trait table.
    pub pruned_tree_tips: Vec<String>,
    /// Rows removed from trait table because they were missing in the tree.
    pub pruned_table_rows: Vec<String>,
    /// Logs of automatic taxonomy adjustments (e.g. trinomial to binomial matches).
    pub taxonomy_matches: Vec<String>,
}

/// Alignment container for species traits and phylogenetic tree matching `caper::comparative.data`.
#[derive(Debug, Clone)]
pub struct ComparativeData {
    pub tree: PhyloTree,
    pub species_names: Vec<String>,
    pub data: TraitTable,
    pub report: PruneReport,
}

#[derive(Debug)]
pub enum ComparativeDataError {
    Io(std::io::Error),
    TreeParse(String),
    ColumnError(String),
    MissingSpeciesColumn,
    PruningError(String),
    SpeciesMismatch {
        matched_count: usize,
        tree_species: Vec<String>,
        table_species: Vec<String>,
    },
}

impl std::fmt::Display for ComparativeDataError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "IO error: {}", e),
            Self::TreeParse(e) => write!(f, "Tree parse error: {}", e),
            Self::ColumnError(e) => write!(f, "Data column error: {}", e),
            Self::MissingSpeciesColumn => write!(f, "Species column not found in trait table"),
            Self::PruningError(e) => write!(f, "Tree pruning error: {}", e),
            Self::SpeciesMismatch {
                matched_count,
                tree_species,
                table_species,
            } => write!(
                f,
                "Species mismatch between tree and trait table. Matched count: {}. Unmatched tree species: {:?}, Unmatched table species: {:?}",
                matched_count, tree_species, table_species
            ),
        }
    }
}

impl std::error::Error for ComparativeDataError {}

impl From<std::io::Error> for ComparativeDataError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

impl ComparativeData {
    /// Constructs `ComparativeData` from Newick string and TSV string with default options.
    pub fn from_newick_and_tsv(newick_str: &str, tsv_str: &str) -> Result<Self, ComparativeDataError> {
        Self::from_newick_and_tsv_with_options(newick_str, tsv_str, &ComparativeDataOptions::default())
    }

    /// Constructs `ComparativeData` from Newick string and TSV string with explicit options.
    pub fn from_newick_and_tsv_with_options(
        newick_str: &str,
        tsv_str: &str,
        options: &ComparativeDataOptions,
    ) -> Result<Self, ComparativeDataError> {
        let tree = parse_newick(newick_str).map_err(ComparativeDataError::TreeParse)?;
        let df = TraitTable::from_reader(std::io::Cursor::new(tsv_str))?;
        Self::align_with_options(tree, df, options)
    }

    /// Reads `ComparativeData` directly from files with default options.
    pub fn from_files<P1: AsRef<Path>, P2: AsRef<Path>>(
        tree_path: P1,
        tsv_path: P2,
    ) -> Result<Self, ComparativeDataError> {
        Self::from_files_with_options(tree_path, tsv_path, &ComparativeDataOptions::default())
    }

    /// Reads `ComparativeData` directly from files with explicit options.
    pub fn from_files_with_options<P1: AsRef<Path>, P2: AsRef<Path>>(
        tree_path: P1,
        tsv_path: P2,
        options: &ComparativeDataOptions,
    ) -> Result<Self, ComparativeDataError> {
        let mut tree_file = File::open(tree_path)?;
        let mut newick_str = String::new();
        tree_file.read_to_string(&mut newick_str)?;

        let tsv_file = File::open(tsv_path)?;
        let df = TraitTable::from_reader(BufReader::with_capacity(64 * 1024, tsv_file))?;

        let tree = parse_newick(&newick_str).map_err(ComparativeDataError::TreeParse)?;
        Self::align_with_options(tree, df, options)
    }

    /// Core alignment engine handling matching, fuzzy taxonomy, pruning, and reporting.
    fn align_with_options(
        tree: PhyloTree,
        df: TraitTable,
        options: &ComparativeDataOptions,
    ) -> Result<Self, ComparativeDataError> {
        let tree_normalized_tips = tree.normalized_tip_names();

        // Extract normalized species list from TraitTable
        let table_species: Vec<String> = if let (Some(g_idx), Some(s_idx)) =
            (df.column_index("Genus"), df.column_index("Species"))
        {
            df.rows
                .iter()
                .map(|r| {
                    let g = r.get(g_idx).map(|s| s.as_str()).unwrap_or_default();
                    let s = r.get(s_idx).map(|s| s.as_str()).unwrap_or_default();
                    normalize_species_name(&format!("{}_{}", g, s))
                })
                .collect()
        } else if let Some(sp_idx) = ["species", "Species", "SPECIES", "species_name", "Spp.Name"]
            .into_iter()
            .find_map(|name| df.column_index(name))
        {
            df.rows
                .iter()
                .map(|r| normalize_species_name(r.get(sp_idx).map(|s| s.as_str()).unwrap_or_default()))
                .collect()
        } else {
            return Err(ComparativeDataError::MissingSpeciesColumn);
        };

        let mut report = PruneReport::default();
        let table_set: HashSet<String> = table_species.iter().cloned().collect();
        let mut matched_tree_to_table: HashMap<String, String> = HashMap::new();
        let mut matched_table_set: HashSet<String> = HashSet::new();

        // Step 1: Direct exact matching
        for tip in &tree_normalized_tips {
            if table_set.contains(tip) {
                matched_tree_to_table.insert(tip.clone(), tip.clone());
                matched_table_set.insert(tip.clone());
            }
        }

        // Step 2: Genus + Subspecies fuzzy taxonomy matching (if enabled)
        if options.fuzzy_taxonomy {
            for tip in &tree_normalized_tips {
                if matched_tree_to_table.contains_key(tip) {
                    continue;
                }

                let tokens: Vec<&str> = tip.split('_').collect();
                if tokens.len() >= 2 {
                    let genus = tokens[0];
                    let sub_sp = tokens[tokens.len() - 1];

                    for tbl_sp in &table_species {
                        if matched_table_set.contains(tbl_sp) {
                            continue;
                        }
                        let tbl_tokens: Vec<&str> = tbl_sp.split('_').collect();
                        if tbl_tokens.len() >= 2
                            && tbl_tokens[0] == genus
                            && tbl_tokens.contains(&sub_sp)
                        {
                            matched_tree_to_table.insert(tip.clone(), tbl_sp.clone());
                            matched_table_set.insert(tbl_sp.clone());
                            report.taxonomy_matches.push(format!(
                                "Fuzzy taxonomy match: tree tip '{}' matched with table row '{}'",
                                tip, tbl_sp
                            ));
                            break;
                        }
                    }
                }
            }
        }

        // Step 3: Identify unmatched species
        for tip in &tree_normalized_tips {
            if !matched_tree_to_table.contains_key(tip) {
                report.pruned_tree_tips.push(tip.clone());
            }
        }
        for tbl_sp in &table_species {
            if !matched_table_set.contains(tbl_sp) {
                report.pruned_table_rows.push(tbl_sp.clone());
            }
        }

        report.matched_species_count = matched_tree_to_table.len();

        // Step 4: Strict validation or Pruning execution
        if !options.prune
            && (!report.pruned_tree_tips.is_empty() || !report.pruned_table_rows.is_empty())
        {
            return Err(ComparativeDataError::SpeciesMismatch {
                matched_count: report.matched_species_count,
                tree_species: report.pruned_tree_tips,
                table_species: report.pruned_table_rows,
            });
        }

        // Apply tree pruning if needed
        let final_tree = if !report.pruned_tree_tips.is_empty() {
            let keep_set: HashSet<String> = matched_tree_to_table.keys().cloned().collect();
            tree.prune(&keep_set)
                .map_err(ComparativeDataError::PruningError)?
        } else {
            tree
        };

        let final_tree_tips = final_tree.normalized_tip_names();

        // Reorder rows to match final pruned tree tip order
        let mut row_indices = Vec::with_capacity(final_tree_tips.len());
        for target_tip in &final_tree_tips {
            let target_table_sp = matched_tree_to_table.get(target_tip).unwrap_or(target_tip);
            let idx = table_species
                .iter()
                .position(|sp| sp == target_table_sp)
                .ok_or_else(|| {
                    ComparativeDataError::PruningError(format!(
                        "Failed to find row for '{}'",
                        target_table_sp
                    ))
                })?;
            row_indices.push(idx);
        }

        let aligned_df = df.reorder(&row_indices);

        Ok(Self {
            tree: final_tree,
            species_names: final_tree_tips,
            data: aligned_df,
            report,
        })
    }

    /// Extracts a response vector $Y$ ($N \times 1$) from a numeric column in `data`.
    pub fn to_response_vector(&self, col_name: &str) -> Result<faer::Mat<f64>, ComparativeDataError> {
        let vals = self.data.column_f64(col_name)?;
        let n = vals.len();
        let mut y = faer::Mat::zeros(n, 1);
        for (i, v) in vals.into_iter().enumerate() {
            y[(i, 0)] = v;
        }
        Ok(y)
    }

    /// Extracts a design matrix $X$ ($N \times p$) with an intercept column of 1s and specified predictor columns.
    pub fn to_design_matrix(&self, predictor_cols: &[&str]) -> Result<faer::Mat<f64>, ComparativeDataError> {
        let n = self.data.height();
        let p = predictor_cols.len() + 1;
        let mut x = faer::Mat::zeros(n, p);

        for i in 0..n {
            x[(i, 0)] = 1.0;
        }

        for (j, &col_name) in predictor_cols.iter().enumerate() {
            let vals = self.data.column_f64(col_name)?;
            for (i, v) in vals.into_iter().enumerate() {
                x[(i, j + 1)] = v;
            }
        }

        Ok(x)
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
        assert_eq!(comp_data.report.matched_species_count, 3);
        assert!(comp_data.report.pruned_tree_tips.is_empty());
    }

    #[test]
    fn test_fuzzy_taxonomy_matching_canis() {
        let newick = "(Canis_lupus_familiaris_protein:0.5,(Homo_sapiens_protein:1,Macaca_mulatta_protein.faa:1)1:0.5);";
        let tsv = "Genus\tSpecies\tCommon name\tBirth weight (g)\n\
                   Canis\tfamiliaris\tDog\t250\n\
                   Homo\tsapiens\tHuman\t3312.5\n\
                   Macaca\tmulatta\tRhesus monkey\t464\n";

        let comp_data = ComparativeData::from_newick_and_tsv(newick, tsv).expect("Fuzzy taxonomy alignment failed");
        assert_eq!(comp_data.report.matched_species_count, 3);
        assert_eq!(comp_data.report.taxonomy_matches.len(), 1);
        assert!(comp_data.report.taxonomy_matches[0].contains("Canis_lupus_familiaris"));
    }

    #[test]
    fn test_pruning_with_report() {
        let newick = "(Canis_lupus_familiaris_protein:0.5,(Homo_sapiens_protein:1,Macaca_mulatta_protein.faa:1)1:0.5);";
        let tsv = "Genus\tSpecies\tCommon name\tBirth weight (g)\n\
                   Homo\tsapiens\tHuman\t3312.5\n\
                   Macaca\tmulatta\tRhesus monkey\t464\n";

        let opts = ComparativeDataOptions {
            prune: true,
            fuzzy_taxonomy: true,
        };

        let comp_data = ComparativeData::from_newick_and_tsv_with_options(newick, tsv, &opts).unwrap();
        assert_eq!(comp_data.report.matched_species_count, 2);
        assert_eq!(comp_data.report.pruned_tree_tips, vec!["Canis_lupus_familiaris"]);
        assert_eq!(comp_data.species_names.len(), 2);
    }
}
