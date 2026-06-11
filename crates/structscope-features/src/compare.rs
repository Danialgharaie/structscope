use structscope_core::{superposition_rmsd, RmsdParams, Structure};

#[derive(Debug, Clone)]
pub struct RmsdMatrix {
    pub labels: Vec<String>,
    pub rmsd: Vec<Vec<Option<f64>>>,
}

pub fn rmsd_matrix(structures: &[Structure], params: &RmsdParams) -> RmsdMatrix {
    let n = structures.len();
    let labels: Vec<String> = structures.iter().map(|s| s.id.clone()).collect();
    let mut rmsd = vec![vec![None; n]; n];

    for i in 0..n {
        rmsd[i][i] = Some(0.0);
    }

    for i in 0..n {
        for j in (i + 1)..n {
            let value = superposition_rmsd(&structures[i], &structures[j], params.clone())
                .ok()
                .map(|result| result.rmsd);
            rmsd[i][j] = value;
            rmsd[j][i] = value;
        }
    }

    RmsdMatrix { labels, rmsd }
}

#[cfg(test)]
mod tests {
    use super::*;
    use structscope_core::{Atom, Chain, Residue, Structure, StructureMetadata};

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

    fn structure(id: &str, residues: &[(&str, [f64; 3])]) -> Structure {
        Structure {
            id: id.to_string(),
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

    fn default_params() -> RmsdParams {
        RmsdParams {
            atoms: "ca".to_string(),
            align: false,
            local: false,
        }
    }

    #[test]
    fn identical_structures_have_zero_off_diagonal_rmsd() {
        let coords = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let residues = [("ALA", coords[0]), ("GLY", coords[1]), ("VAL", coords[2])];
        let structures = [
            structure("s1", &residues),
            structure("s2", &residues),
            structure("s3", &residues),
        ];

        let matrix = rmsd_matrix(&structures, &default_params());

        assert_eq!(matrix.labels, vec!["s1", "s2", "s3"]);
        assert_eq!(matrix.rmsd.len(), 3);
        for i in 0..3 {
            assert_eq!(matrix.rmsd[i][i], Some(0.0));
            for j in 0..3 {
                if i != j {
                    let value = matrix.rmsd[i][j].expect("expected rmsd value");
                    assert!(value < 1e-9, "rmsd[{i}][{j}] = {value}");
                    assert_eq!(matrix.rmsd[j][i], matrix.rmsd[i][j]);
                }
            }
        }
    }

    #[test]
    fn offset_structures_have_nonzero_rmsd() {
        let reference = structure(
            "ref",
            &[
                ("ALA", [0.0, 0.0, 0.0]),
                ("GLY", [1.0, 0.0, 0.0]),
                ("VAL", [0.0, 1.0, 0.0]),
            ],
        );
        let mobile = structure(
            "mob",
            &[
                ("ALA", [0.0, 0.0, 0.0]),
                ("GLY", [2.0, 0.0, 0.0]),
                ("VAL", [0.0, 1.0, 0.0]),
            ],
        );

        let matrix = rmsd_matrix(&[reference, mobile], &default_params());

        assert_eq!(matrix.labels, vec!["ref", "mob"]);
        let off_diagonal = matrix.rmsd[0][1].expect("expected rmsd value");
        assert!(off_diagonal > 0.0, "rmsd = {off_diagonal}");
        assert_eq!(matrix.rmsd[1][0], Some(off_diagonal));
    }
}
