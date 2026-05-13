use petgraph::graph::UnGraph;
use serde::{Deserialize, Serialize};
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
pub struct ContactEdge {
    pub distance: f64,
    pub kind: String,
}

pub type ResidueGraph = UnGraph<ResidueNode, ContactEdge>;

pub fn build_residue_graph(structure: &Structure, threshold_angstroms: f64) -> ResidueGraph {
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

    graph
}

pub fn export_graphml(graph: &ResidueGraph) -> String {
    let mut out = String::from(
        r#"<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns">
  <graph edgedefault="undirected">
"#,
    );
    for node_idx in graph.node_indices() {
        let node = &graph[node_idx];
        out.push_str(&format!(
            "    <node id=\"{}\"><data key=\"residue_name\">{}</data><data key=\"seq_number\">{}</data></node>\n",
            node.residue_id, node.residue_name, node.seq_number
        ));
    }
    for edge_idx in graph.edge_indices() {
        let (source, target) = graph.edge_endpoints(edge_idx).expect("edge endpoints");
        let edge = &graph[edge_idx];
        out.push_str(&format!(
            "    <edge source=\"{}\" target=\"{}\"><data key=\"kind\">{}</data><data key=\"distance\">{:.3}</data></edge>\n",
            graph[source].residue_id, graph[target].residue_id, edge.kind, edge.distance
        ));
    }
    out.push_str("  </graph>\n</graphml>\n");
    out
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
