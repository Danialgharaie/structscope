use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use structscope_core::Structure;
use structscope_graphs::{contacting_chain_pairs, min_heavy_atom_residue_distance};

use crate::lc89::{interface_surface_dots, lc89_score};
use crate::sasa::{atom_buried_sasa_deltas, atom_sasa, atom_sasa_chain_neighbors};

pub struct InterfaceParams {
    pub contact_distance: f64,
    pub area_distance: f64,
    pub sc_distance: f64,
}

#[derive(Debug, Clone, Default)]
pub struct ProteinInterfaceSummary {
    pub interface_pair_count: usize,
    pub interface_bsa_total: f64,
    pub interface_area_total: f64,
    pub interface_sc_mean: f64,
    pub interface_bsa_max: f64,
    pub interface_area_max: f64,
    pub interface_sc_max: f64,
    pub interface_chain_a: String,
    pub interface_chain_b: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InterfaceFeature {
    pub structure_id: String,
    pub chain_a: String,
    pub chain_b: String,
    pub contact_count: usize,
    pub interface_residue_count_a: usize,
    pub interface_residue_count_b: usize,
    pub bsa: f64,
    pub interface_area: f64,
    pub shape_complementarity: f64,
}

pub fn protein_interface_summary(_structure: &Structure, _params: &InterfaceParams) -> ProteinInterfaceSummary {
    let records = per_interface_features(_structure, _params);
    if records.is_empty() {
        return ProteinInterfaceSummary::default();
    }

    let interface_pair_count = records.len();
    let interface_bsa_total: f64 = records.iter().map(|r| r.bsa).sum();
    let interface_area_total: f64 = records.iter().map(|r| r.interface_area).sum();
    let interface_sc_mean = records.iter().map(|r| r.shape_complementarity).sum::<f64>() / interface_pair_count as f64;

    let best = records
        .iter()
        .max_by(|left, right| {
            left.bsa
                .partial_cmp(&right.bsa)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| right.chain_a.cmp(&left.chain_a))
                .then_with(|| right.chain_b.cmp(&left.chain_b))
        })
        .expect("non-empty");

    ProteinInterfaceSummary {
        interface_pair_count,
        interface_bsa_total,
        interface_area_total,
        interface_sc_mean,
        interface_bsa_max: best.bsa,
        interface_area_max: best.interface_area,
        interface_sc_max: best.shape_complementarity,
        interface_chain_a: best.chain_a.clone(),
        interface_chain_b: best.chain_b.clone(),
    }
}

pub fn per_interface_features(_structure: &Structure, _params: &InterfaceParams) -> Vec<InterfaceFeature> {
    if _structure.chains.len() < 2 {
        return Vec::new();
    }

    #[derive(Clone)]
    struct AtomMeta {
        chain_label: String,
        residue_id: String,
        xyz: [f64; 3],
        is_heavy: bool,
    }

    fn is_heavy(element: Option<&str>, atom_name: &str) -> bool {
        if let Some(elem) = element {
            return !elem.trim().eq_ignore_ascii_case("H");
        }
        !atom_name.trim().starts_with('H')
    }

    let mut atom_meta = Vec::new();
    for chain in &_structure.chains {
        for residue in &chain.residues {
            for atom in &residue.atoms {
                atom_meta.push(AtomMeta {
                    chain_label: chain.label.clone(),
                    residue_id: residue.id.clone(),
                    xyz: [atom.x, atom.y, atom.z],
                    is_heavy: is_heavy(atom.element.as_deref(), &atom.name),
                });
            }
        }
    }
    let atom_chain_labels: Vec<String> = atom_meta.iter().map(|a| a.chain_label.clone()).collect();

    let complex_sasa = atom_sasa(_structure);
    let mut records = Vec::new();
    for (chain_a, chain_b) in contacting_chain_pairs(_structure, _params.contact_distance) {
        let mono_a = atom_sasa_chain_neighbors(_structure, &chain_a);
        let mono_b = atom_sasa_chain_neighbors(_structure, &chain_b);
        let deltas = atom_buried_sasa_deltas(
            &complex_sasa,
            &mono_a,
            &mono_b,
            &atom_chain_labels,
            &chain_a,
            &chain_b,
        );
        let bsa: f64 = deltas.iter().sum();

        let area2 = _params.area_distance * _params.area_distance;
        let mut interface_area = 0.0;
        let mut patch_res_a = HashSet::new();
        let mut patch_res_b = HashSet::new();

        for (i, atom) in atom_meta.iter().enumerate() {
            if !atom.is_heavy {
                continue;
            }
            let is_pair_atom = atom.chain_label == chain_a || atom.chain_label == chain_b;
            if !is_pair_atom {
                continue;
            }
            let partner_label = if atom.chain_label == chain_a { &chain_b } else { &chain_a };
            let in_patch = atom_meta.iter().any(|other| {
                if !other.is_heavy || &other.chain_label != partner_label {
                    return false;
                }
                let dx = atom.xyz[0] - other.xyz[0];
                let dy = atom.xyz[1] - other.xyz[1];
                let dz = atom.xyz[2] - other.xyz[2];
                dx * dx + dy * dy + dz * dz <= area2
            });
            if in_patch {
                interface_area += deltas.get(i).copied().unwrap_or(0.0);
                if atom.chain_label == chain_a {
                    patch_res_a.insert(atom.residue_id.clone());
                } else {
                    patch_res_b.insert(atom.residue_id.clone());
                }
            }
        }

        let chain_a_ref = _structure.chains.iter().find(|c| c.label == chain_a);
        let chain_b_ref = _structure.chains.iter().find(|c| c.label == chain_b);
        let mut contact_count = 0;
        if let (Some(a), Some(b)) = (chain_a_ref, chain_b_ref) {
            for ra in &a.residues {
                for rb in &b.residues {
                    if min_heavy_atom_residue_distance(ra, rb) <= _params.contact_distance {
                        contact_count += 1;
                    }
                }
            }
        }

        let dots_a = interface_surface_dots(_structure, &chain_a, &chain_b, _params.sc_distance);
        let dots_b = interface_surface_dots(_structure, &chain_b, &chain_a, _params.sc_distance);
        // We use a single LC89 score over paired directional dot sets.
        let shape_complementarity = lc89_score(&dots_a, &dots_b, _params.sc_distance);

        records.push(InterfaceFeature {
            structure_id: _structure.id.clone(),
            chain_a,
            chain_b,
            contact_count,
            interface_residue_count_a: patch_res_a.len(),
            interface_residue_count_b: patch_res_b.len(),
            bsa,
            interface_area,
            shape_complementarity,
        });
    }
    records
}

#[cfg(test)]
mod tests {
    use super::*;
    use structscope_core::{parse_str, InputFormat, ParseOptions};

    const DIMER: &str = "\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 0.00           C
ATOM      2  CA  ALA B   1       3.500   0.000   0.000  1.00 0.00           C
";

    #[test]
    fn summary_reports_one_contacting_pair() {
        let s = parse_str(DIMER, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let params = InterfaceParams {
            contact_distance: 8.0,
            area_distance: 5.0,
            sc_distance: 5.0,
        };
        let summary = protein_interface_summary(&s, &params);
        assert_eq!(summary.interface_pair_count, 1);
        assert!(summary.interface_bsa_total > 0.0);
        assert!(summary.interface_area_total > 0.0);
        assert!(summary.interface_sc_mean > 0.0);
        assert_eq!(summary.interface_bsa_max, summary.interface_bsa_total);
    }

    #[test]
    fn per_pair_record_shape() {
        let s = parse_str(DIMER, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let params = InterfaceParams {
            contact_distance: 8.0,
            area_distance: 5.0,
            sc_distance: 5.0,
        };
        let records = per_interface_features(&s, &params);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].chain_a, "A");
        assert_eq!(records[0].chain_b, "B");
        assert!(records[0].bsa > 0.0);
    }
}
