//! Per-residue structural features, combining SASA, secondary structure, and
//! backbone dihedrals into one record per residue. All three sources iterate
//! chains -> residues in the same order, so they align positionally.
use crate::{dihedral::backbone_dihedrals, sasa::atom_sasa, ss::secondary_structure};
use serde::{Deserialize, Serialize};
use structscope_core::Structure;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResidueFeature {
    pub residue_id: String,
    pub chain_id: String,
    pub residue_name: String,
    pub seq_number: i32,
    pub sasa: f64,
    pub secondary_structure: char,
    pub phi: Option<f64>,
    pub psi: Option<f64>,
    pub omega: Option<f64>,
}

/// Compute one feature record per residue, in chains -> residues order.
pub fn per_residue_features(structure: &Structure) -> Vec<ResidueFeature> {
    let sasa = atom_sasa(structure);
    let ss = secondary_structure(structure);
    let dihedrals = backbone_dihedrals(structure);

    let mut out = Vec::new();
    let mut atom_i = 0;
    let mut res_i = 0;
    for (ci, chain) in structure.chains.iter().enumerate() {
        let ss_chars: Vec<char> = ss[ci].ss.chars().collect();
        for (ri, residue) in chain.residues.iter().enumerate() {
            let mut res_sasa = 0.0;
            for _ in &residue.atoms {
                res_sasa += sasa[atom_i];
                atom_i += 1;
            }
            let d = &dihedrals[res_i];
            out.push(ResidueFeature {
                residue_id: residue.id.clone(),
                chain_id: chain.id.clone(),
                residue_name: residue.name.clone(),
                seq_number: residue.seq_number,
                sasa: res_sasa,
                secondary_structure: ss_chars[ri],
                phi: d.phi,
                psi: d.psi,
                omega: d.omega,
            });
            res_i += 1;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use structscope_core::{parse_str, InputFormat, ParseOptions};

    #[test]
    fn one_record_per_residue() {
        let pdb = "\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 0.00           C
ATOM      4  N   GLY A   2       3.332   1.540   0.000  1.00 0.00           N
ATOM      5  CA  GLY A   2       3.970   2.840   0.000  1.00 0.00           C
ATOM      6  C   GLY A   2       5.480   2.700   0.000  1.00 0.00           C
";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let feats = per_residue_features(&s);
        assert_eq!(feats.len(), 2);
        assert_eq!(feats[0].residue_name, "ALA");
        assert_eq!(feats[1].seq_number, 2);
        assert!(feats.iter().all(|f| f.sasa >= 0.0));
    }
}
