//! Newick format parser module.

use super::node::{Node, PhyloTree};
use nom::{
    bytes::complete::take_while1,
    character::complete::{char, space0},
    combinator::{map, opt},
    multi::separated_list1,
    number::complete::double,
    sequence::{delimited, preceded, tuple},
    IResult,
};

#[derive(Debug)]
struct RawNode {
    name: Option<String>,
    branch_length: Option<f64>,
    children: Vec<RawNode>,
}

fn parse_name(input: &str) -> IResult<&str, String> {
    map(
        take_while1(|c: char| !c.is_whitespace() && c != ':' && c != ',' && c != ')' && c != '(' && c != ';'),
        |s: &str| s.to_string(),
    )(input)
}

fn parse_branch_length(input: &str) -> IResult<&str, f64> {
    preceded(char(':'), double)(input)
}

fn parse_node(input: &str) -> IResult<&str, RawNode> {
    let (input, _) = space0(input)?;
    if input.starts_with('(') {
        // Internal node: ( child1, child2, ... ) [name] [:branch_length]
        let (input, children) = delimited(
            char('('),
            separated_list1(char(','), parse_node),
            char(')'),
        )(input)?;
        let (input, name) = opt(parse_name)(input)?;
        let (input, branch_length) = opt(parse_branch_length)(input)?;
        let (input, _) = space0(input)?;
        Ok((
            input,
            RawNode {
                name,
                branch_length,
                children,
            },
        ))
    } else {
        // Tip node: name [:branch_length]
        let (input, (name, branch_length)) = tuple((opt(parse_name), opt(parse_branch_length)))(input)?;
        let (input, _) = space0(input)?;
        Ok((
            input,
            RawNode {
                name,
                branch_length,
                children: Vec::new(),
            },
        ))
    }
}

/// Parses a Newick format string into a `PhyloTree`.
pub fn parse_newick(input: &str) -> Result<PhyloTree, String> {
    let trimmed = input.trim().trim_end_matches(';');
    let (remaining, raw_root) = parse_node(trimmed).map_err(|e| format!("Newick parse error: {:?}", e))?;

    if !remaining.trim().is_empty() {
        return Err(format!("Unexpected trailing characters: '{}'", remaining));
    }

    let mut nodes = Vec::new();
    let mut tips = Vec::new();

    fn build_nodes(
        raw: RawNode,
        parent: Option<usize>,
        nodes: &mut Vec<Node>,
        tips: &mut Vec<usize>,
    ) -> usize {
        let id = nodes.len();
        let is_tip = raw.children.is_empty();

        // Push temporary node
        nodes.push(Node {
            id,
            name: raw.name,
            branch_length: raw.branch_length,
            children: Vec::new(),
            parent,
            is_tip,
        });

        if is_tip {
            tips.push(id);
        } else {
            let child_ids: Vec<usize> = raw
                .children
                .into_iter()
                .map(|child| build_nodes(child, Some(id), nodes, tips))
                .collect();
            nodes[id].children = child_ids;
        }

        id
    }

    let root = build_nodes(raw_root, None, &mut nodes, &mut tips);

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
