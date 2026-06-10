use petgraph::graph::{NodeIndex, UnGraph};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use structscope_core::{Residue, Structure};

#[derive(Debug, Clone, Copy)]
pub enum GraphType {
    Residue,
    Atom,
    Interface,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResidueNode {
    pub residue_id: String,
    pub chain_id: String,
    pub residue_name: String,
    pub seq_number: i32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AtomNode {
    pub atom_id: String,
    pub residue_id: String,
    pub chain_id: String,
    pub atom_name: String,
    pub element: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContactEdge {
    pub distance: f64,
    pub kind: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChemicalInteraction {
    pub kind: String,
    pub res_id_a: String,
    pub res_id_b: String,
    pub distance: f64,
}

pub fn atom_id_to_residue_id(atom_id: &str) -> Option<String> {
    let parts: Vec<&str> = atom_id.split(':').collect();
    if parts.len() >= 3 {
        let residue_parts = &parts[..parts.len() - 2];
        Some(residue_parts.join(":"))
    } else {
        None
    }
}

fn interaction_priority(kind: &str) -> i32 {
    match kind {
        "disulfide" => 7,
        "salt_bridge" => 6,
        "hydrogen_bond" => 5,
        "cation_pi" => 4,
        "pi_pi_parallel" => 3,
        "pi_pi_perpendicular" => 3,
        "hydrophobic" => 2,
        "distance_contact" | "interface_contact" => 1,
        _ => 0, // covalent_adjacency or unknown
    }
}

pub type ResidueGraph = UnGraph<ResidueNode, ContactEdge>;
pub type AtomGraph = UnGraph<AtomNode, ContactEdge>;

pub fn build_residue_graph(
    structure: &Structure,
    threshold_angstroms: f64,
    interactions: Option<&[ChemicalInteraction]>,
) -> ResidueGraph {
    let mut graph = ResidueGraph::default();
    let mut entries = Vec::new();

    for chain in &structure.chains {
        for residue in &chain.residues {
            let index = graph.add_node(ResidueNode {
                residue_id: residue.id.clone(),
                chain_id: chain.id.clone(),
                residue_name: residue.name.clone(),
                seq_number: residue.seq_number,
            });
            entries.push((chain.id.as_str(), residue, index));
        }
    }

    for idx in 0..entries.len() {
        if idx + 1 < entries.len() {
            let (left_chain, left_residue, left_index) = entries[idx];
            let (right_chain, right_residue, right_index) = entries[idx + 1];
            if left_chain == right_chain && are_sequential_neighbors(left_residue, right_residue) {
                graph.add_edge(
                    left_index,
                    right_index,
                    ContactEdge {
                        distance: residue_distance(left_residue, right_residue),
                        kind: "covalent_adjacency".to_string(),
                    },
                );
            }
        }
    }

    for left in 0..entries.len() {
        for right in (left + 1)..entries.len() {
            let (_, left_residue, left_index) = entries[left];
            let (_, right_residue, right_index) = entries[right];
            let distance = residue_distance(left_residue, right_residue);
            if distance <= threshold_angstroms && graph.find_edge(left_index, right_index).is_none() {
                graph.add_edge(
                    left_index,
                    right_index,
                    ContactEdge {
                        distance,
                        kind: "distance_contact".to_string(),
                    },
                );
            }
        }
    }

    if let Some(interactions) = interactions {
        let mut id_to_index = HashMap::new();
        for (_, residue, index) in &entries {
            id_to_index.insert(residue.id.clone(), *index);
        }

        for interaction in interactions {
            if let (Some(&a_idx), Some(&b_idx)) = (
                id_to_index.get(&interaction.res_id_a),
                id_to_index.get(&interaction.res_id_b),
            ) {
                if let Some(edge_idx) = graph.find_edge(a_idx, b_idx) {
                    let edge = &mut graph[edge_idx];
                    let existing_prio = interaction_priority(&edge.kind);
                    let new_prio = interaction_priority(&interaction.kind);
                    if edge.kind != "covalent_adjacency"
                        && (new_prio > existing_prio
                            || (new_prio == existing_prio && interaction.distance < edge.distance))
                    {
                        edge.kind = interaction.kind.clone();
                        edge.distance = interaction.distance;
                    }
                } else {
                    graph.add_edge(
                        a_idx,
                        b_idx,
                        ContactEdge {
                            distance: interaction.distance,
                            kind: interaction.kind.clone(),
                        },
                    );
                }
            }
        }
    }

    graph
}

/// Residue graph containing only inter-chain contacts. Residues without an
/// interface contact are omitted so the graph represents the interface itself.
pub fn build_interface_graph(
    structure: &Structure,
    threshold_angstroms: f64,
    interactions: Option<&[ChemicalInteraction]>,
) -> ResidueGraph {
    let residues: Vec<(&str, &Residue)> = structure
        .chains
        .iter()
        .flat_map(|chain| chain.residues.iter().map(move |residue| (chain.id.as_str(), residue)))
        .collect();

    let mut graph = ResidueGraph::default();
    let mut node_index: HashMap<usize, NodeIndex> = HashMap::new();
    let ensure_node = |graph: &mut ResidueGraph, map: &mut HashMap<usize, NodeIndex>, i: usize| -> NodeIndex {
        *map.entry(i).or_insert_with(|| {
            let (chain_id, residue) = residues[i];
            graph.add_node(ResidueNode {
                residue_id: residue.id.clone(),
                chain_id: chain_id.to_string(),
                residue_name: residue.name.clone(),
                seq_number: residue.seq_number,
            })
        })
    };

    for left in 0..residues.len() {
        for right in (left + 1)..residues.len() {
            if residues[left].0 == residues[right].0 {
                continue;
            }
            let distance = residue_distance(residues[left].1, residues[right].1);
            if distance <= threshold_angstroms {
                let a = ensure_node(&mut graph, &mut node_index, left);
                let b = ensure_node(&mut graph, &mut node_index, right);
                graph.add_edge(a, b, ContactEdge { distance, kind: "interface_contact".to_string() });
            }
        }
    }

    if let Some(interactions) = interactions {
        let mut residue_id_to_slice_index = HashMap::new();
        for (idx, (_, residue)) in residues.iter().enumerate() {
            residue_id_to_slice_index.insert(residue.id.as_str(), idx);
        }

        for interaction in interactions {
            if let (Some(&left_idx), Some(&right_idx)) = (
                residue_id_to_slice_index.get(interaction.res_id_a.as_str()),
                residue_id_to_slice_index.get(interaction.res_id_b.as_str()),
            ) {
                if residues[left_idx].0 != residues[right_idx].0 {
                    let a = ensure_node(&mut graph, &mut node_index, left_idx);
                    let b = ensure_node(&mut graph, &mut node_index, right_idx);

                    if let Some(edge_idx) = graph.find_edge(a, b) {
                        let edge = &mut graph[edge_idx];
                        let existing_prio = interaction_priority(&edge.kind);
                        let new_prio = interaction_priority(&interaction.kind);
                        if new_prio > existing_prio
                            || (new_prio == existing_prio && interaction.distance < edge.distance)
                        {
                            edge.kind = interaction.kind.clone();
                            edge.distance = interaction.distance;
                        }
                    } else {
                        graph.add_edge(
                            a,
                            b,
                            ContactEdge {
                                distance: interaction.distance,
                                kind: interaction.kind.clone(),
                            },
                        );
                    }
                }
            }
        }
    }

    graph
}

/// Atom-level contact graph using a spatial grid to find pairs within the threshold.
pub fn build_atom_graph(structure: &Structure, threshold_angstroms: f64) -> AtomGraph {
    let mut graph = AtomGraph::default();
    let mut points = Vec::new();
    let mut indices = Vec::new();

    for chain in &structure.chains {
        for residue in &chain.residues {
            for atom in &residue.atoms {
                let index = graph.add_node(AtomNode {
                    atom_id: atom.id.clone(),
                    residue_id: residue.id.clone(),
                    chain_id: chain.id.clone(),
                    atom_name: atom.name.clone(),
                    element: atom.element.clone().unwrap_or_default(),
                });
                points.push((atom.x, atom.y, atom.z));
                indices.push(index);
            }
        }
    }

    for (i, j, distance) in grid_pairs(&points, threshold_angstroms) {
        graph.add_edge(indices[i], indices[j], ContactEdge { distance, kind: "distance_contact".to_string() });
    }

    graph
}

/// Node types that can be serialized to GraphML.
pub trait GraphmlNode {
    fn graphml_id(&self) -> &str;
    fn graphml_attrs(&self) -> Vec<(&'static str, String)>;
}

impl GraphmlNode for ResidueNode {
    fn graphml_id(&self) -> &str {
        &self.residue_id
    }
    fn graphml_attrs(&self) -> Vec<(&'static str, String)> {
        vec![
            ("residue_name", self.residue_name.clone()),
            ("seq_number", self.seq_number.to_string()),
        ]
    }
}

impl GraphmlNode for AtomNode {
    fn graphml_id(&self) -> &str {
        &self.atom_id
    }
    fn graphml_attrs(&self) -> Vec<(&'static str, String)> {
        vec![
            ("atom_name", self.atom_name.clone()),
            ("element", self.element.clone()),
            ("residue_id", self.residue_id.clone()),
        ]
    }
}

pub fn export_graphml<N: GraphmlNode>(graph: &UnGraph<N, ContactEdge>) -> String {
    let mut out = String::from(
        r#"<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns">
  <graph edgedefault="undirected">
"#,
    );
    for node_idx in graph.node_indices() {
        let node = &graph[node_idx];
        out.push_str(&format!("    <node id=\"{}\">", node.graphml_id()));
        for (key, value) in node.graphml_attrs() {
            out.push_str(&format!("<data key=\"{key}\">{value}</data>"));
        }
        out.push_str("</node>\n");
    }
    for edge_idx in graph.edge_indices() {
        let (source, target) = graph.edge_endpoints(edge_idx).expect("edge endpoints");
        let edge = &graph[edge_idx];
        out.push_str(&format!(
            "    <edge source=\"{}\" target=\"{}\"><data key=\"kind\">{}</data><data key=\"distance\">{:.3}</data></edge>\n",
            graph[source].graphml_id(), graph[target].graphml_id(), edge.kind, edge.distance
        ));
    }
    out.push_str("  </graph>\n</graphml>\n");
    out
}

pub fn export_gml<N: GraphmlNode>(graph: &UnGraph<N, ContactEdge>) -> String {
    let mut out = String::from("graph [\n  directed 0\n");
    for node_idx in graph.node_indices() {
        let node = &graph[node_idx];
        let id = node_idx.index();
        out.push_str("  node [\n");
        out.push_str(&format!("    id {}\n", id));
        out.push_str(&format!("    label \"{}\"\n", node.graphml_id()));
        for (key, value) in node.graphml_attrs() {
            if let Ok(val_i) = value.parse::<i32>() {
                out.push_str(&format!("    {} {}\n", key, val_i));
            } else if let Ok(val_f) = value.parse::<f64>() {
                out.push_str(&format!("    {} {}\n", key, val_f));
            } else {
                out.push_str(&format!("    {} \"{}\"\n", key, value));
            }
        }
        out.push_str("  ]\n");
    }
    for edge_idx in graph.edge_indices() {
        let (source, target) = graph.edge_endpoints(edge_idx).expect("edge endpoints");
        let edge = &graph[edge_idx];
        out.push_str("  edge [\n");
        out.push_str(&format!("    source {}\n", source.index()));
        out.push_str(&format!("    target {}\n", target.index()));
        out.push_str(&format!("    distance {:.3}\n", edge.distance));
        out.push_str(&format!("    kind \"{}\"\n", edge.kind));
        out.push_str("  ]\n");
    }
    out.push_str("]\n");
    out
}

pub fn export_json<N: serde::Serialize>(graph: &UnGraph<N, ContactEdge>) -> String {
    #[derive(Serialize)]
    struct JsonGraph<'a, N> {
        nodes: Vec<JsonNode<'a, N>>,
        links: Vec<JsonLink<'a>>,
    }

    #[derive(Serialize)]
    struct JsonNode<'a, N> {
        id: usize,
        #[serde(flatten)]
        data: &'a N,
    }

    #[derive(Serialize)]
    struct JsonLink<'a> {
        source: usize,
        target: usize,
        #[serde(flatten)]
        data: &'a ContactEdge,
    }

    let mut nodes = Vec::new();
    for node_idx in graph.node_indices() {
        nodes.push(JsonNode {
            id: node_idx.index(),
            data: &graph[node_idx],
        });
    }

    let mut links = Vec::new();
    for edge_idx in graph.edge_indices() {
        let (source, target) = graph.edge_endpoints(edge_idx).expect("edge endpoints");
        links.push(JsonLink {
            source: source.index(),
            target: target.index(),
            data: &graph[edge_idx],
        });
    }

    let json_graph = JsonGraph { nodes, links };
    serde_json::to_string_pretty(&json_graph).unwrap_or_default()
}

/// Return index pairs (i < j) whose points lie within `threshold` using cell hashing.
fn grid_pairs(points: &[(f64, f64, f64)], threshold: f64) -> Vec<(usize, usize, f64)> {
    let cell = threshold.max(f64::MIN_POSITIVE);
    let key = |p: &(f64, f64, f64)| {
        ((p.0 / cell).floor() as i64, (p.1 / cell).floor() as i64, (p.2 / cell).floor() as i64)
    };
    let mut cells: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    for (i, p) in points.iter().enumerate() {
        cells.entry(key(p)).or_default().push(i);
    }

    let t2 = threshold * threshold;
    let mut pairs = Vec::new();
    for (i, p) in points.iter().enumerate() {
        let (cx, cy, cz) = key(p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(bucket) = cells.get(&(cx + dx, cy + dy, cz + dz)) {
                        for &j in bucket {
                            if j <= i {
                                continue;
                            }
                            let q = &points[j];
                            let d2 = (p.0 - q.0).powi(2) + (p.1 - q.1).powi(2) + (p.2 - q.2).powi(2);
                            if d2 <= t2 {
                                pairs.push((i, j, d2.sqrt()));
                            }
                        }
                    }
                }
            }
        }
    }
    pairs
}

fn are_sequential_neighbors(left: &Residue, right: &Residue) -> bool {
    (left.seq_number - right.seq_number).abs() == 1
}

fn residue_distance(left: &Residue, right: &Residue) -> f64 {
    let (lx, ly, lz) = residue_centroid(left);
    let (rx, ry, rz) = residue_centroid(right);
    let dx = lx - rx;
    let dy = ly - ry;
    let dz = lz - rz;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn residue_centroid(residue: &Residue) -> (f64, f64, f64) {
    let count = residue.atoms.len() as f64;
    if count == 0.0 {
        return (0.0, 0.0, 0.0);
    }
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_z = 0.0;
    for atom in &residue.atoms {
        sum_x += atom.x;
        sum_y += atom.y;
        sum_z += atom.z;
    }
    (sum_x / count, sum_y / count, sum_z / count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use structscope_core::{parse_str, InputFormat, ParseOptions};

    const ATOMS: &str = "\
ATOM      1  N   GLY A   1      11.104  13.207   8.292  1.00 20.00           N
ATOM      2  CA  GLY A   1      12.000  12.500   8.000  1.00 20.00           C
ATOM      3  C   GLY A   2      13.100  12.800   8.900  1.00 20.00           C
";

    const TWO_CHAIN: &str = "\
ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 0.00           C
ATOM      2  CA  GLY A   2       1.000   0.000   0.000  1.00 0.00           C
ATOM      3  CA  ALA B   1       2.000   0.000   0.000  1.00 0.00           C
";

    fn parse(pdb: &str) -> Structure {
        parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap()
    }

    #[test]
    fn atom_graph_has_node_per_atom_and_contacts() {
        let graph = build_atom_graph(&parse(ATOMS), 5.0);
        assert_eq!(graph.node_count(), 3);
        assert_eq!(graph.edge_count(), 3); // all three atoms within 5A
    }

    #[test]
    fn interface_graph_only_inter_chain() {
        let graph = build_interface_graph(&parse(TWO_CHAIN), 8.0, None);
        // A:1-B:1 and A:2-B:1 cross the interface; A:1-A:2 is intra-chain and excluded.
        assert_eq!(graph.edge_count(), 2);
        assert_eq!(graph.node_count(), 3);
    }

    #[test]
    fn interface_graph_empty_for_single_chain() {
        let graph = build_interface_graph(&parse(ATOMS), 8.0, None);
        assert_eq!(graph.edge_count(), 0);
        assert_eq!(graph.node_count(), 0);
    }

    #[test]
    fn build_residue_graph_with_chemical_interactions() {
        let structure = parse(TWO_CHAIN);
        let interactions = vec![
            ChemicalInteraction {
                kind: "hydrogen_bond".to_string(),
                res_id_a: format!("{}:A:1:_", structure.id),
                res_id_b: format!("{}:A:2:_", structure.id),
                distance: 2.8,
            }
        ];
        let graph = build_residue_graph(&structure, 8.0, Some(&interactions));
        let edge_idx = graph.find_edge(NodeIndex::new(0), NodeIndex::new(1)).unwrap();
        assert_eq!(graph[edge_idx].kind, "covalent_adjacency");
    }
}
