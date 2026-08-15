//! Fast streaming reader and tripartite aligner for OrthoFinder GeneCount tables.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::data::comparative::{ComparativeData, ComparativeDataOptions};
use crate::tree::node::normalize_species_name;
use crate::tree::PhyloTree;

/// In-memory table representation of OrthoFinder `Orthogroups.GeneCount.tsv`.
#[derive(Debug, Clone)]
pub struct OrthoGeneCountTable {
    /// Identifiers for all orthogroups (e.g. `OG0000000`, `OG0000001`, ...).
    pub orthogroup_ids: Vec<String>,
    /// Normalized species names from table columns.
    pub species_names: Vec<String>,
    /// Flat contiguous matrix of shape (M_orthogroups x K_species).
    /// Access count for orthogroup `m` and species `k` as `counts[m * K + k]`.
    pub counts: Vec<f64>,
}

impl OrthoGeneCountTable {
    /// Number of orthogroups (rows).
    pub fn num_orthogroups(&self) -> usize {
        self.orthogroup_ids.len()
    }

    /// Number of species (columns).
    pub fn num_species(&self) -> usize {
        self.species_names.len()
    }

    /// Returns a slice of gene counts across all species for a given orthogroup index.
    // ponytail: zero-copy row slice access
    pub fn get_og_counts(&self, og_idx: usize) -> &[f64] {
        let k = self.species_names.len();
        let start = og_idx * k;
        &self.counts[start..start + k]
    }
}

/// Fast streaming parser for OrthoFinder `Orthogroups.GeneCount.tsv`.
// ponytail: fast zero-allocation line reader for OrthoFinder TSVs
pub fn read_orthogroup_gene_counts(path: &str) -> Result<OrthoGeneCountTable, String> {
    let file = File::open(path).map_err(
        |e| format!("Failed to open OrthoFinder file '{}': {}", path, e)
    )?;
    let mut reader = BufReader::with_capacity(
        64 * 1024, file
    );

    let mut header_line = String::new();
    if reader.read_line(&mut header_line).map_err(
        |e| e.to_string()
    )? == 0 {
        return Err("OrthoFinder file is empty".to_string());
    }

    let header_cols: Vec<&str> = header_line.trim_end().split('\t').collect();
    if header_cols.len() < 2 {
        return Err("OrthoFinder file header must have at least 2 columns".to_string());
    }

    // Determine species column indices (skip 'Orthogroup' at col 0 and 'Total' at last col if present)
    let end_idx = if header_cols.last().map(
        |s| s.eq_ignore_ascii_case("Total")
    ).unwrap_or(false) {
        header_cols.len() - 1
    } else {
        header_cols.len()
    };

    let species_names: Vec<String> = header_cols[1..end_idx]
        .iter()
        .map(|s| normalize_species_name(s))
        .collect();

    let num_species = species_names.len();
    let mut orthogroup_ids = Vec::with_capacity(30000);
    let mut counts = Vec::with_capacity(30000 * num_species);

    let mut line = String::new();
    while reader.read_line(&mut line).map_err(|e| e.to_string())? > 0 {
        let trimmed = line.trim_end();
        if !trimmed.is_empty() {
            let mut parts = trimmed.split('\t');
            if let Some(og_id) = parts.next() {
                orthogroup_ids.push(og_id.to_string());
                for (sp_idx, val_str) in parts.enumerate() {
                    if sp_idx < num_species {
                        let count = val_str.parse::<f64>().unwrap_or(0.0);
                        counts.push(count);
                    }
                }
            }
        }
        line.clear();
    }

    Ok(OrthoGeneCountTable {
        orthogroup_ids,
        species_names,
        counts,
    })
}

/// Summary report of the tripartite alignment between Tree, Phenotype, and OrthoFinder.
#[derive(Debug, Clone, Default)]
pub struct TripartiteReport {
    pub tree_tips_count: usize,
    pub pheno_rows_count: usize,
    pub orthogroup_species_count: usize,
    pub total_orthogroups: usize,
    pub matched_species_count: usize,
    pub pruned_tree_tips: Vec<String>,
    pub pruned_pheno_species: Vec<String>,
    pub pruned_ortho_species: Vec<String>,
    pub taxonomy_matches: Vec<String>,
}

/// Aligned comparative dataset containing Tree, Phenotype table, and OrthoFinder counts for common species.
#[derive(Debug, Clone)]
pub struct TripartiteComparativeData {
    pub tree: PhyloTree,
    pub species_names: Vec<String>,
    pub pheno_data: polars::prelude::DataFrame,
    pub orthogroups: OrthoGeneCountTable,
    pub report: TripartiteReport,
}

/// Aligns Tree, Phenotype table (e.g. AnAge), and OrthoFinder GeneCount table into a unified structure.
// ponytail: single-pass tripartite intersection and tree pruning
pub fn align_tripartite(
    tree_path: &str,
    pheno_path: &str,
    pheno_options: &ComparativeDataOptions,
    ortho_path: &str,
) -> Result<TripartiteComparativeData, String> {
    // 1. Ingest ComparativeData (Tree + Phenotype)
    let comp_data = ComparativeData::from_files_with_options(tree_path, pheno_path, pheno_options)
        .map_err(|e| format!("Phenotype alignment error: {}", e))?;

    // 2. Ingest OrthoFinder GeneCount table
    let raw_ortho = read_orthogroup_gene_counts(ortho_path)?;

    // 3. Match species between ComparativeData and OrthoFinder
    let mut ortho_sp_map: HashMap<String, usize> = HashMap::new();
    for (idx, sp) in raw_ortho.species_names.iter().enumerate() {
        ortho_sp_map.insert(sp.clone(), idx);
    }

    let mut matched_species = Vec::new();
    let mut matched_tree_indices = Vec::new();
    let mut matched_ortho_indices = Vec::new();
    let mut taxonomy_matches = comp_data.report.taxonomy_matches.clone();

    for (comp_idx, sp) in comp_data.species_names.iter().enumerate() {
        if let Some(&ortho_idx) = ortho_sp_map.get(sp) {
            matched_species.push(sp.clone());
            matched_tree_indices.push(comp_idx);
            matched_ortho_indices.push(ortho_idx);
        } else if pheno_options.fuzzy_taxonomy {
            // Try fuzzy matching: Genus + Subspecies
            let parts: Vec<&str> = sp.split('_').collect();
            let fuzzy_match = if parts.len() >= 3 {
                let candidate = format!("{}_{}", parts[0], parts.last().unwrap());
                ortho_sp_map.get(&candidate).copied().map(|idx| (idx, candidate))
            } else if parts.len() == 2 {
                ortho_sp_map.iter().find_map(|(name, &idx)| {
                    let o_parts: Vec<&str> = name.split('_').collect();
                    if o_parts.len() >= 3 && o_parts[0] == parts[0] && o_parts.last().unwrap() == &parts[1] {
                        Some((idx, name.clone()))
                    } else {
                        None
                    }
                })
            } else {
                None
            };

            if let Some((ortho_idx, matched_name)) = fuzzy_match {
                taxonomy_matches.push(format!(
                    "OrthoFinder fuzzy match: pheno/tree '{}' matched with ortho '{}'",
                    sp, matched_name
                ));
                matched_species.push(sp.clone());
                matched_tree_indices.push(comp_idx);
                matched_ortho_indices.push(ortho_idx);
            }
        }
    }

    if matched_species.is_empty() {
        return Err("No common species found across Tree, Phenotype table, and OrthoFinder data".to_string());
    }

    let n_matched = matched_species.len();

    // 4. Prune Tree to the tripartite common species
    let matched_species_set: HashSet<String> = matched_species.iter().cloned().collect();
    let final_tree = comp_data.tree.prune(&matched_species_set)?;

    // Get tips order from pruned tree
    let final_tree_tips: Vec<String> = final_tree
        .tips
        .iter()
        .map(|&idx| {
            final_tree.nodes[idx]
                .name
                .as_deref()
                .map(normalize_species_name)
                .unwrap_or_default()
        })
        .collect();

    // Map each tree tip to its index in matched_species
    let tree_sp_to_idx: HashMap<String, usize> = matched_species
        .iter()
        .enumerate()
        .map(|(i, s)| (s.clone(), i))
        .collect();

    let mut reordered_ortho_indices = Vec::with_capacity(n_matched);
    let mut reordered_tree_indices = Vec::with_capacity(n_matched);

    for tip in &final_tree_tips {
        if let Some(&match_idx) = tree_sp_to_idx.get(tip) {
            reordered_ortho_indices.push(matched_ortho_indices[match_idx]);
            reordered_tree_indices.push(matched_tree_indices[match_idx]);
        }
    }

    // 5. Re-index OrthoFinder counts matrix aligned to final tree tips order
    let num_ogs = raw_ortho.num_orthogroups();
    let raw_k = raw_ortho.num_species();
    let mut aligned_counts = Vec::with_capacity(num_ogs * n_matched);

    for m in 0..num_ogs {
        let og_offset = m * raw_k;
        for &ortho_sp_idx in &reordered_ortho_indices {
            aligned_counts.push(raw_ortho.counts[og_offset + ortho_sp_idx]);
        }
    }

    let aligned_ortho = OrthoGeneCountTable {
        orthogroup_ids: raw_ortho.orthogroup_ids,
        species_names: final_tree_tips.clone(),
        counts: aligned_counts,
    };

    // 6. Filter Phenotype DataFrame rows aligned to final tree tips order
    use polars::prelude::*;
    let row_indices: Vec<u32> = reordered_tree_indices.iter().map(|&idx| idx as u32).collect();
    let idx_ca = UInt32Chunked::from_slice("idx".into(), &row_indices);
    let aligned_pheno = comp_data.data.take(&idx_ca).map_err(|e| format!("Polars take error: {}", e))?;

    let report = TripartiteReport {
        tree_tips_count: comp_data.tree.tips.len(),
        pheno_rows_count: comp_data.data.height(),
        orthogroup_species_count: raw_ortho.species_names.len(),
        total_orthogroups: num_ogs,
        matched_species_count: n_matched,
        pruned_tree_tips: comp_data.report.pruned_tree_tips,
        pruned_pheno_species: comp_data.report.pruned_table_rows,
        pruned_ortho_species: raw_ortho
            .species_names
            .into_iter()
            .filter(|s| !matched_species_set.contains(s))
            .collect(),
        taxonomy_matches,
    };

    Ok(TripartiteComparativeData {
        tree: final_tree,
        species_names: final_tree_tips,
        pheno_data: aligned_pheno,
        orthogroups: aligned_ortho,
        report,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_orthogroup_gene_counts_mock() {
        let tsv_content = "Orthogroup\tSp_A\tSp_B\tTotal\nOG0001\t2\t5\t7\nOG0002\t0\t3\t3\n";
        let temp_path = "agy/test/mock_ortho.tsv";
        std::fs::write(temp_path, tsv_content).unwrap();

        let table = read_orthogroup_gene_counts(temp_path).unwrap();
        assert_eq!(table.num_orthogroups(), 2);
        assert_eq!(table.num_species(), 2);
        assert_eq!(table.species_names, vec!["Sp_A", "Sp_B"]);
        assert_eq!(table.get_og_counts(0), &[2.0, 5.0]);
        assert_eq!(table.get_og_counts(1), &[0.0, 3.0]);

        let _ = std::fs::remove_file(temp_path);
    }
}
