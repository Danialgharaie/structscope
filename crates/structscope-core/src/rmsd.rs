//! RMSD superposition between two structures with optional sequence-based pairing.

use crate::align::kabsch;
use crate::model::Structure;
use crate::seqalign::{needleman_wunsch, smith_waterman, three_to_one};
use thiserror::Error;

/// Parameters controlling atom selection and residue correspondence for RMSD.
#[derive(Debug, Clone)]
pub struct RmsdParams {
    /// Atom selection label: `ca`, `backbone`, or `all` (any other value selects all atoms).
    pub atoms: String,
    /// Establish residue correspondence by global sequence alignment (CA atoms).
    pub align: bool,
    /// Like `align` but uses local (Smith-Waterman) alignment for partial overlaps.
    pub local: bool,
}

/// Result of optimal-superposition RMSD between two structures.
#[derive(Debug, Clone, PartialEq)]
pub struct RmsdResult {
    pub rmsd: f64,
    pub atom_count: usize,
    pub selection: String,
}

#[derive(Debug, Error)]
pub enum RmsdError {
    #[error("no matching residues found between the two structures")]
    NoMatchingResidues,
    #[error(
        "atom count mismatch under selection '{selection}': reference has {reference}, mobile has {mobile} (use --align for sequence-based correspondence)"
    )]
    AtomCountMismatch {
        selection: String,
        reference: usize,
        mobile: usize,
    },
    #[error("superposition failed (empty selection?)")]
    SuperpositionFailed,
}

fn atom_matches_selection(name: &str, atoms: &str) -> bool {
    match atoms {
        "ca" => name == "CA",
        "backbone" => matches!(name, "N" | "CA" | "C" | "O"),
        _ => true,
    }
}

fn ca_sequence_and_coords(structure: &Structure) -> (Vec<u8>, Vec<[f64; 3]>) {
    let mut seq = Vec::new();
    let mut ca = Vec::new();
    for r in structure.chains.iter().flat_map(|c| &c.residues) {
        if let Some(a) = r.atoms.iter().find(|a| a.name == "CA") {
            seq.push(three_to_one(&r.name));
            ca.push([a.x, a.y, a.z]);
        }
    }
    (seq, ca)
}

fn select_coords(structure: &Structure, atoms: &str) -> Vec<[f64; 3]> {
    structure
        .chains
        .iter()
        .flat_map(|c| &c.residues)
        .flat_map(|r| &r.atoms)
        .filter(|a| atom_matches_selection(&a.name, atoms))
        .map(|a| [a.x, a.y, a.z])
        .collect()
}

fn paired_coords(reference: &Structure, mobile: &Structure, local: bool) -> Result<(Vec<[f64; 3]>, Vec<[f64; 3]>), RmsdError> {
    let (ref_seq, ref_ca) = ca_sequence_and_coords(reference);
    let (mob_seq, mob_ca) = ca_sequence_and_coords(mobile);
    let pairs = if local {
        smith_waterman(&ref_seq, &mob_seq)
    } else {
        needleman_wunsch(&ref_seq, &mob_seq)
    };
    let matched: Vec<(usize, usize)> = pairs
        .into_iter()
        .filter(|&(i, j)| ref_seq[i] == mob_seq[j])
        .collect();
    if matched.is_empty() {
        return Err(RmsdError::NoMatchingResidues);
    }
    Ok((
        matched.iter().map(|&(i, _)| ref_ca[i]).collect(),
        matched.iter().map(|&(_, j)| mob_ca[j]).collect(),
    ))
}

/// Compute optimal-superposition RMSD between `reference` and `mobile`.
pub fn superposition_rmsd(
    reference: &Structure,
    mobile: &Structure,
    params: RmsdParams,
) -> Result<RmsdResult, RmsdError> {
    let (ref_coords, mob_coords) = if params.align || params.local {
        paired_coords(reference, mobile, params.local)?
    } else {
        let r = select_coords(reference, &params.atoms);
        let m = select_coords(mobile, &params.atoms);
        if r.len() != m.len() {
            return Err(RmsdError::AtomCountMismatch {
                selection: params.atoms.clone(),
                reference: r.len(),
                mobile: m.len(),
            });
        }
        (r, m)
    };

    let superposition = kabsch(&mob_coords, &ref_coords).ok_or(RmsdError::SuperpositionFailed)?;
    let selection = if params.local {
        "local-ca".to_string()
    } else if params.align {
        "aligned-ca".to_string()
    } else {
        params.atoms.clone()
    };

    Ok(RmsdResult {
        rmsd: superposition.rmsd,
        atom_count: ref_coords.len(),
        selection,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Atom, Chain, Residue, Structure, StructureMetadata};

    fn ca_atom(x: f64, y: f64, z: f64) -> Atom {
        Atom {
            id: "CA".to_string(),
            serial: Some(1),
            name: "CA".to_string(),
            element: Some("C".to_string()),
            x,
            y,
            z,
            occupancy: None,
            temp_factor: None,
        }
    }

    fn residue(name: &str, ca: [f64; 3]) -> Residue {
        Residue {
            id: name.to_string(),
            name: name.to_string(),
            seq_number: 1,
            insertion_code: None,
            atoms: vec![ca_atom(ca[0], ca[1], ca[2])],
            is_hetero: false,
        }
    }

    fn structure(residues: &[(&str, [f64; 3])]) -> Structure {
        Structure {
            id: "test".to_string(),
            metadata: StructureMetadata {
                source_format: "test".to_string(),
                source_path: None,
                title: None,
            },
            chains: vec![Chain {
                id: "A".to_string(),
                label: "A".to_string(),
                residues: residues
                    .iter()
                    .map(|(name, ca)| residue(name, *ca))
                    .collect(),
            }],
        }
    }

    #[test]
    fn equal_length_ca_zero_rmsd() {
        let coords = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let reference = structure(&[("ALA", coords[0]), ("GLY", coords[1]), ("VAL", coords[2])]);
        let mobile = structure(&[("ALA", coords[0]), ("GLY", coords[1]), ("VAL", coords[2])]);
        let result = superposition_rmsd(
            &reference,
            &mobile,
            RmsdParams {
                atoms: "ca".to_string(),
                align: false,
                local: false,
            },
        )
        .unwrap();
        assert!(result.rmsd < 1e-9, "rmsd {}", result.rmsd);
        assert_eq!(result.atom_count, 3);
        assert_eq!(result.selection, "ca");
    }

    #[test]
    fn aligned_different_length_structures() {
        let reference = structure(&[
            ("ALA", [0.0, 0.0, 0.0]),
            ("GLY", [1.0, 0.0, 0.0]),
            ("VAL", [2.0, 0.0, 0.0]),
            ("LEU", [3.0, 0.0, 0.0]),
            ("ILE", [4.0, 0.0, 0.0]),
        ]);
        let mobile = structure(&[
            ("ALA", [0.0, 0.0, 0.0]),
            ("GLY", [1.0, 0.0, 0.0]),
            ("ILE", [4.0, 0.0, 0.0]),
        ]);
        let result = superposition_rmsd(
            &reference,
            &mobile,
            RmsdParams {
                atoms: "ca".to_string(),
                align: true,
                local: false,
            },
        )
        .unwrap();
        assert!(result.rmsd < 1e-9, "rmsd {}", result.rmsd);
        assert_eq!(result.atom_count, 3);
        assert_eq!(result.selection, "aligned-ca");
    }
}
