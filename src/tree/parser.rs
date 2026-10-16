//! Newick format parser module.

use super::node::{Node, PhyloTree};

fn parse_label_and_length(s: &str) -> Result<(Option<String>, Option<f64>), String> {
    let trimmed = s.trim();
    if trimmed.is_empty() {
        return Ok((None, None));
    }
    if let Some((name_part, len_part)) = trimmed.split_once(':') {
        let name = (!name_part.trim().is_empty()).then(|| name_part.trim().to_string());
        let bl = len_part
            .trim()
            .parse::<f64>()
            .map_err(|e| format!("Invalid branch length '{}': {}", len_part, e))?;
        Ok((name, Some(bl)))
    } else {
        Ok((Some(trimmed.to_string()), None))
    }
}

fn parse_subtree(
    input: &str,
    parent: Option<usize>,
    nodes: &mut Vec<Node>,
    tips: &mut Vec<usize>,
) -> Result<usize, String> {
    let s = input.trim();
    if s.is_empty() {
        return Err("Empty Newick node".to_string());
    }

    let id = nodes.len();
    if s.starts_with('(') {
        let mut depth = 0usize;
        let mut close_idx = None;
        let mut split_positions = Vec::new();

        for (i, ch) in s.char_indices() {
            match ch {
                '(' => depth += 1,
                ')' => {
                    depth -= 1;
                    if depth == 0 {
                        close_idx = Some(i);
                        break;
                    }
                }
                ',' if depth == 1 => split_positions.push(i),
                _ => {}
            }
        }

        let close_pos = close_idx.ok_or_else(|| "Unmatched '(' in Newick string".to_string())?;
        let inner = &s[1..close_pos];
        let meta = &s[close_pos + 1..];
        let (name, branch_length) = parse_label_and_length(meta)?;

        nodes.push(Node {
            id,
            name,
            branch_length,
            children: Vec::new(),
            parent,
            is_tip: false,
        });

        let mut children = Vec::with_capacity(split_positions.len() + 1);
        let mut start = 0;
        for pos in split_positions {
            let child_str = &inner[start..pos - 1];
            children.push(parse_subtree(child_str, Some(id), nodes, tips)?);
            start = pos;
        }
        children.push(parse_subtree(&inner[start..], Some(id), nodes, tips)?);
        nodes[id].children = children;
    } else {
        if s.contains(['(', ')', ',', ';']) {
            return Err(format!("Unexpected delimiter in tip '{}'", s));
        }
        let (name, branch_length) = parse_label_and_length(s)?;
        nodes.push(Node {
            id,
            name,
            branch_length,
            children: Vec::new(),
            parent,
            is_tip: true,
        });
        tips.push(id);
    }

    Ok(id)
}

/// Parses a Newick format string into a `PhyloTree`.
pub fn parse_newick(input: &str) -> Result<PhyloTree, String> {
    let trimmed = input.trim();
    let body = trimmed.strip_suffix(';').unwrap_or(trimmed).trim();
    if body.is_empty() {
        return Err("Empty Newick tree".to_string());
    }

    let mut nodes = Vec::new();
    let mut tips = Vec::new();
    let root = parse_subtree(body, None, &mut nodes, &mut tips)?;

    Ok(PhyloTree { nodes, root, tips })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::node::normalize_species_name;

    #[test]
    fn test_normalize_species_name() {
        assert_eq!(normalize_species_name("Homo sapiens"), "Homo_sapiens");
        assert_eq!(normalize_species_name("Homo_sapiens_protein"), "Homo_sapiens");
        assert_eq!(normalize_species_name("Macaca_mulatta_protein.faa"), "Macaca_mulatta");
        assert_eq!(normalize_species_name("Canis_lupus_familiaris"), "Canis_lupus_familiaris");
    }

    #[test]
    fn test_parse_newick_3sp() {
        let input = "((SpA:1.0,SpB:1.0):1.0,SpC:2.0);";
        let tree = parse_newick(input).expect("Failed to parse Newick");
        assert_eq!(tree.tips.len(), 3);
        let tip_names = tree.tip_names();
        assert_eq!(tip_names, vec!["SpA", "SpB", "SpC"]);
    }

    #[test]
    fn test_parse_newick_benchmark_3sp() {
        let input = "(Canis_lupus_familiaris_protein:0.5,(Homo_sapiens_protein:1,Macaca_mulatta_protein.faa:1)1:0.5);";
        let tree = parse_newick(input).expect("Failed to parse benchmark tree");
        assert_eq!(tree.tips.len(), 3);
        let normalized = tree.normalized_tip_names();
        assert_eq!(normalized, vec!["Canis_lupus_familiaris", "Homo_sapiens", "Macaca_mulatta"]);
    }
}
