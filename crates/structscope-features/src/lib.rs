use serde::{Deserialize, Serialize};
use serde_json::{json, Map, Value};
use structscope_core::Structure;
use structscope_graphs::{build_interface_graph, build_residue_graph};

pub mod dihedral;
pub mod interactions;
pub mod per_residue;
pub mod sasa;
pub mod ss;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FeatureRecord {
    pub structure_id: String,
    pub source_path: Option<String>,
    pub features: Map<String, Value>,
}

pub fn compute_features(structure: &Structure) -> FeatureRecord {
    let summary = structure.summary();
    let residue_graph = build_residue_graph(structure, 8.0);
    let interface_graph = build_interface_graph(structure, 8.0);

    let mut features = Map::new();
    features.insert("atom_count".to_string(), json!(summary.atom_count));
    features.insert("residue_count".to_string(), json!(summary.residue_count));
    features.insert("chain_count".to_string(), json!(summary.chain_count));
    features.insert("ligand_count".to_string(), json!(summary.ligand_count));
    features.insert("heteroatom_count".to_string(), json!(summary.heteroatom_count));
    features.insert("contact_count".to_string(), json!(residue_graph.edge_count()));
    features.insert("interface_contact_count".to_string(), json!(interface_graph.edge_count()));
    features.insert("interface_residue_count".to_string(), json!(interface_graph.node_count()));
    features.insert("radius_of_gyration".to_string(), json!(radius_of_gyration(structure)));
    features.insert("sasa_total".to_string(), json!(sasa::total_sasa(structure)));
    // Buried vs exposed by relative accessibility (standard 25% cutoff); residues
    // without a reference area (non-standard) are excluded from both counts.
    let rsa: Vec<f64> = per_residue::per_residue_features(structure).iter().filter_map(|r| r.rsa).collect();
    features.insert("buried_residue_count".to_string(), json!(rsa.iter().filter(|&&v| v < 0.25).count()));
    features.insert("exposed_residue_count".to_string(), json!(rsa.iter().filter(|&&v| v >= 0.25).count()));
    let ss_all: String = ss::secondary_structure(structure).iter().map(|c| c.ss.clone()).collect();
    features.insert(
        "helix_residue_count".to_string(),
        json!(ss_all.chars().filter(|c| matches!(c, 'H' | 'G' | 'I')).count()),
    );
    features.insert("strand_residue_count".to_string(), json!(ss_all.chars().filter(|&c| c == 'E').count()));
    features.insert("coil_residue_count".to_string(), json!(ss_all.chars().filter(|&c| c == 'C').count()));
    let contacts = interactions::interactions(structure);
    features.insert(
        "disulfide_count".to_string(),
        json!(contacts.iter().filter(|i| i.kind == "disulfide").count()),
    );
    features.insert(
        "salt_bridge_count".to_string(),
        json!(contacts.iter().filter(|i| i.kind == "salt_bridge").count()),
    );
    features.insert(
        "hydrogen_bond_count".to_string(),
        json!(contacts.iter().filter(|i| i.kind == "hydrogen_bond").count()),
    );
    features.insert("centroid".to_string(), json!(structure_centroid(structure)));
    features.insert(
        "graph_density".to_string(),
        json!(graph_density(residue_graph.node_count(), residue_graph.edge_count())),
    );
    features.insert(
        "connected_components".to_string(),
        json!(petgraph::algo::connected_components(&residue_graph)),
    );
    features.insert("clustering_coefficient".to_string(), json!(clustering_coefficient(&residue_graph)));
    features.insert("degree_distribution".to_string(), json!(degree_distribution(&residue_graph)));

    FeatureRecord {
        structure_id: structure.id.clone(),
        source_path: structure.metadata.source_path.clone(),
        features,
    }
}

fn structure_centroid(structure: &Structure) -> [f64; 3] {
    let mut total = 0.0;
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_z = 0.0;

    for atom in structure
        .chains
        .iter()
        .flat_map(|chain| &chain.residues)
        .flat_map(|residue| &residue.atoms)
    {
        total += 1.0;
        sum_x += atom.x;
        sum_y += atom.y;
        sum_z += atom.z;
    }

    if total == 0.0 {
        return [0.0, 0.0, 0.0];
    }

    [sum_x / total, sum_y / total, sum_z / total]
}

fn radius_of_gyration(structure: &Structure) -> f64 {
    let centroid = structure_centroid(structure);
    let mut total = 0.0;
    let mut sum = 0.0;

    for atom in structure
        .chains
        .iter()
        .flat_map(|chain| &chain.residues)
        .flat_map(|residue| &residue.atoms)
    {
        total += 1.0;
        let dx = atom.x - centroid[0];
        let dy = atom.y - centroid[1];
        let dz = atom.z - centroid[2];
        sum += dx * dx + dy * dy + dz * dz;
    }

    if total == 0.0 {
        0.0
    } else {
        (sum / total).sqrt()
    }
}

fn graph_density(nodes: usize, edges: usize) -> f64 {
    if nodes < 2 {
        return 0.0;
    }
    let max_edges = (nodes * (nodes - 1) / 2) as f64;
    edges as f64 / max_edges
}

fn degree_distribution(graph: &structscope_graphs::ResidueGraph) -> Vec<usize> {
    graph.node_indices().map(|idx| graph.neighbors(idx).count()).collect()
}

fn clustering_coefficient(graph: &structscope_graphs::ResidueGraph) -> f64 {
    if graph.node_count() == 0 {
        return 0.0;
    }

    let mut total = 0.0;
    for node in graph.node_indices() {
        let neighbors: Vec<_> = graph.neighbors(node).collect();
        let k = neighbors.len();
        if k < 2 {
            continue;
        }
        let mut triangles = 0usize;
        for left in 0..k {
            for right in (left + 1)..k {
                if graph.find_edge(neighbors[left], neighbors[right]).is_some() {
                    triangles += 1;
                }
            }
        }
        let possible = k * (k - 1) / 2;
        total += triangles as f64 / possible as f64;
    }

    total / graph.node_count() as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use structscope_core::{parse_str, InputFormat, ParseOptions};

    const PDB_SAMPLE: &str = "\
ATOM      1  N   GLY A   1      11.104  13.207   8.292  1.00 20.00           N
ATOM      2  CA  GLY A   1      12.000  12.500   8.000  1.00 20.00           C
ATOM      3  C   GLY A   2      13.100  12.800   8.900  1.00 20.00           C
";

    #[test]
    fn computes_basic_features() {
        let structure = parse_str(PDB_SAMPLE, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let record = compute_features(&structure);
        assert_eq!(record.features["atom_count"].as_u64(), Some(3));
        assert_eq!(record.features["residue_count"].as_u64(), Some(2));
        assert!(record.features["radius_of_gyration"].as_f64().unwrap() >= 0.0);
    }
}
