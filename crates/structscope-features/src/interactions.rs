//! Typed geometric interactions detected from coordinates and atom identity.
//! Raw distance/geometry rules; emits a typed contact list for the caller to use.
use structscope_core::Structure;

/// One detected interaction between two atoms, with the separating distance (Angstrom).
pub struct Interaction {
    pub kind: &'static str,
    pub atom_id_a: String,
    pub atom_id_b: String,
    pub distance: f64,
}

struct A {
    id: String,
    res_id: String,
    res: String,
    name: String,
    x: f64,
    y: f64,
    z: f64,
}

fn dist(a: &A, b: &A) -> f64 {
    ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt()
}

/// True for sidechain carboxylate oxygens (ASP/GLU) — salt-bridge anions.
pub(crate) fn is_acidic_oxygen(res: &str, name: &str) -> bool {
    matches!((res, name), ("ASP", "OD1") | ("ASP", "OD2") | ("GLU", "OE1") | ("GLU", "OE2"))
}

/// True for sidechain basic nitrogens (LYS/ARG/HIS) — salt-bridge cations.
pub(crate) fn is_basic_nitrogen(res: &str, name: &str) -> bool {
    matches!(
        (res, name),
        ("LYS", "NZ") | ("ARG", "NH1") | ("ARG", "NH2") | ("ARG", "NE") | ("HIS", "ND1") | ("HIS", "NE2")
    )
}

pub(crate) fn is_hydrophobic_carbon(res: &str, name: &str) -> bool {
    let r = res.trim();
    let n = name.trim();
    if matches!(r, "ALA" | "VAL" | "LEU" | "ILE" | "PRO" | "MET" | "PHE" | "TYR" | "TRP") {
        n.starts_with('C') && n != "C" && n != "CA"
    } else {
        false
    }
}

pub(crate) fn is_hydrogen_bond(a_name: &str, b_name: &str, distance: f64) -> bool {
    let elem = |n: &str| n.chars().next().unwrap_or(' ');
    matches!(elem(a_name), 'N' | 'O')
        && matches!(elem(b_name), 'N' | 'O')
        && (2.4..=3.5).contains(&distance)
}

pub(crate) fn is_salt_bridge(a_res: &str, a_name: &str, b_res: &str, b_name: &str, distance: f64) -> bool {
    ((is_acidic_oxygen(a_res, a_name) && is_basic_nitrogen(b_res, b_name))
        || (is_acidic_oxygen(b_res, b_name) && is_basic_nitrogen(a_res, a_name)))
        && distance < 4.0
}

struct AromaticRing {
    res_id: String,
    centroid: [f64; 3],
    normal: [f64; 3],
}

fn cross_product(u: [f64; 3], v: [f64; 3]) -> [f64; 3] {
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}

fn normalize(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt();
    if len > 0.0 {
        [v[0] / len, v[1] / len, v[2] / len]
    } else {
        [0.0, 0.0, 1.0]
    }
}

/// Detect disulfides (CYS SG-SG < 2.5 A), salt bridges (acidic O / basic N < 4.0 A),
/// and polar contacts (N/O donor-acceptor 2.4-3.5 A).
/// Also detects cation-pi (< 6.0 A), parallel/perpendicular pi-pi stacking (< 5.5 A),
/// and hydrophobic carbon-carbon contacts (< 4.5 A).
pub fn interactions(structure: &Structure) -> Vec<Interaction> {
    let atoms: Vec<A> = structure
        .chains
        .iter()
        .flat_map(|c| &c.residues)
        .flat_map(|r| {
            r.atoms.iter().map(move |a| A {
                id: a.id.clone(),
                res_id: r.id.clone(),
                res: r.name.trim().to_string(),
                name: a.name.trim().to_string(),
                x: a.x,
                y: a.y,
                z: a.z,
            })
        })
        .collect();

    let mut out = Vec::new();

    // 1. Classical atom-atom interactions (disulfides, salt bridges, hydrogen bonds)
    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            let (a, b) = (&atoms[i], &atoms[j]);
            let d = dist(a, b);
            if d > 4.0 {
                continue;
            }

            if a.res == "CYS" && a.name == "SG" && b.res == "CYS" && b.name == "SG" && d < 2.5 {
                out.push(Interaction { kind: "disulfide", atom_id_a: a.id.clone(), atom_id_b: b.id.clone(), distance: d });
                continue;
            }

            if is_salt_bridge(&a.res, &a.name, &b.res, &b.name, d) {
                out.push(Interaction { kind: "salt_bridge", atom_id_a: a.id.clone(), atom_id_b: b.id.clone(), distance: d });
                continue;
            }

            // Polar contact: a donor/acceptor pair of N/O atoms on different residues.
            if a.res_id != b.res_id && is_hydrogen_bond(&a.name, &b.name, d) {
                out.push(Interaction { kind: "hydrogen_bond", atom_id_a: a.id.clone(), atom_id_b: b.id.clone(), distance: d });
            }
        }
    }

    // 2. Build Aromatic Rings
    let mut rings = Vec::new();
    for chain in &structure.chains {
        for residue in &chain.residues {
            let res_name = residue.name.trim();
            if res_name == "PHE" || res_name == "TYR" {
                let cg = residue.atoms.iter().find(|a| a.name.trim() == "CG");
                let cd1 = residue.atoms.iter().find(|a| a.name.trim() == "CD1");
                let cd2 = residue.atoms.iter().find(|a| a.name.trim() == "CD2");
                let ce1 = residue.atoms.iter().find(|a| a.name.trim() == "CE1");
                let ce2 = residue.atoms.iter().find(|a| a.name.trim() == "CE2");
                let cz = residue.atoms.iter().find(|a| a.name.trim() == "CZ");
                if let (Some(cg), Some(cd1), Some(cd2), Some(ce1), Some(ce2), Some(cz)) = (cg, cd1, cd2, ce1, ce2, cz) {
                    let centroid = [
                        (cg.x + cd1.x + cd2.x + ce1.x + ce2.x + cz.x) / 6.0,
                        (cg.y + cd1.y + cd2.y + ce1.y + ce2.y + cz.y) / 6.0,
                        (cg.z + cd1.z + cd2.z + ce1.z + ce2.z + cz.z) / 6.0,
                    ];
                    let u = [cd1.x - cg.x, cd1.y - cg.y, cd1.z - cg.z];
                    let v = [ce2.x - cg.x, ce2.y - cg.y, ce2.z - cg.z];
                    let normal = normalize(cross_product(u, v));
                    rings.push(AromaticRing {
                        res_id: residue.id.clone(),
                        centroid,
                        normal,
                    });
                }
            } else if res_name == "TRP" {
                let cg = residue.atoms.iter().find(|a| a.name.trim() == "CG");
                let cd1 = residue.atoms.iter().find(|a| a.name.trim() == "CD1");
                let cd2 = residue.atoms.iter().find(|a| a.name.trim() == "CD2");
                let ne1 = residue.atoms.iter().find(|a| a.name.trim() == "NE1");
                let ce2 = residue.atoms.iter().find(|a| a.name.trim() == "CE2");
                let ce3 = residue.atoms.iter().find(|a| a.name.trim() == "CE3");
                let cz2 = residue.atoms.iter().find(|a| a.name.trim() == "CZ2");
                let cz3 = residue.atoms.iter().find(|a| a.name.trim() == "CZ3");
                let ch2 = residue.atoms.iter().find(|a| a.name.trim() == "CH2");
                if let (Some(cg), Some(cd1), Some(cd2), Some(ne1), Some(ce2), Some(ce3), Some(cz2), Some(cz3), Some(ch2)) =
                    (cg, cd1, cd2, ne1, ce2, ce3, cz2, cz3, ch2)
                {
                    let centroid = [
                        (cg.x + cd1.x + cd2.x + ne1.x + ce2.x + ce3.x + cz2.x + cz3.x + ch2.x) / 9.0,
                        (cg.y + cd1.y + cd2.y + ne1.y + ce2.y + ce3.y + cz2.y + cz3.y + ch2.y) / 9.0,
                        (cg.z + cd1.z + cd2.z + ne1.z + ce2.z + ce3.z + cz2.z + cz3.z + ch2.z) / 9.0,
                    ];
                    let u = [cd1.x - cg.x, cd1.y - cg.y, cd1.z - cg.z];
                    let v = [ce3.x - cg.x, ce3.y - cg.y, ce3.z - cg.z];
                    let normal = normalize(cross_product(u, v));
                    rings.push(AromaticRing {
                        res_id: residue.id.clone(),
                        centroid,
                        normal,
                    });
                }
            }
        }
    }

    // 3. Cation-pi interactions
    for ring in &rings {
        for a in &atoms {
            if a.res_id != ring.res_id && is_basic_nitrogen(&a.res, &a.name) {
                let dx = a.x - ring.centroid[0];
                let dy = a.y - ring.centroid[1];
                let dz = a.z - ring.centroid[2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                if d < 6.0 {
                    out.push(Interaction {
                        kind: "cation_pi",
                        atom_id_a: a.id.clone(),
                        atom_id_b: format!("{}:CG", ring.res_id),
                        distance: d,
                    });
                }
            }
        }
    }

    // 4. Pi-pi stacking
    for i in 0..rings.len() {
        for j in (i + 1)..rings.len() {
            let r1 = &rings[i];
            let r2 = &rings[j];
            if r1.res_id != r2.res_id {
                let dx = r1.centroid[0] - r2.centroid[0];
                let dy = r1.centroid[1] - r2.centroid[1];
                let dz = r1.centroid[2] - r2.centroid[2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                if d < 5.5 {
                    let dot = r1.normal[0] * r2.normal[0] + r1.normal[1] * r2.normal[1] + r1.normal[2] * r2.normal[2];
                    let cos_theta = dot.abs().min(1.0);
                    let theta = cos_theta.acos().to_degrees();
                    if theta < 30.0 {
                        out.push(Interaction {
                            kind: "pi_pi_parallel",
                            atom_id_a: format!("{}:CG", r1.res_id),
                            atom_id_b: format!("{}:CG", r2.res_id),
                            distance: d,
                        });
                    } else if theta >= 60.0 && theta <= 90.0 {
                        out.push(Interaction {
                            kind: "pi_pi_perpendicular",
                            atom_id_a: format!("{}:CG", r1.res_id),
                            atom_id_b: format!("{}:CG", r2.res_id),
                            distance: d,
                        });
                    }
                }
            }
        }
    }

    // 5. Hydrophobic contacts
    let hydrophobic_atoms: Vec<&A> = atoms
        .iter()
        .filter(|a| is_hydrophobic_carbon(&a.res, &a.name))
        .collect();
    for i in 0..hydrophobic_atoms.len() {
        for j in (i + 1)..hydrophobic_atoms.len() {
            let a = hydrophobic_atoms[i];
            let b = hydrophobic_atoms[j];
            if a.res_id != b.res_id {
                let d = dist(a, b);
                if d < 4.5 {
                    out.push(Interaction {
                        kind: "hydrophobic",
                        atom_id_a: a.id.clone(),
                        atom_id_b: b.id.clone(),
                        distance: d,
                    });
                }
            }
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use structscope_core::{parse_str, InputFormat, ParseOptions};

    #[test]
    fn detects_disulfide() {
        // Two CYS SG atoms 2.05 A apart -> one disulfide.
        let pdb = "\
ATOM      1  SG  CYS A   1       0.000   0.000   0.000  1.00 0.00           S
ATOM      2  SG  CYS A  10       2.050   0.000   0.000  1.00 0.00           S
";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let found = interactions(&s);
        assert_eq!(found.iter().filter(|i| i.kind == "disulfide").count(), 1);
    }

    #[test]
    fn detects_salt_bridge() {
        // ASP OD1 to LYS NZ at 3.0 A -> one salt bridge.
        let pdb = "\
ATOM      1  OD1 ASP A   1       0.000   0.000   0.000  1.00 0.00           O
ATOM      2  NZ  LYS A   5       3.000   0.000   0.000  1.00 0.00           N
";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let found = interactions(&s);
        assert_eq!(found.iter().filter(|i| i.kind == "salt_bridge").count(), 1);
    }

    #[test]
    fn far_atoms_yield_nothing() {
        let pdb = "\
ATOM      1  SG  CYS A   1       0.000   0.000   0.000  1.00 0.00           S
ATOM      2  SG  CYS A  10      10.000   0.000   0.000  1.00 0.00           S
";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        assert!(interactions(&s).is_empty());
    }

    #[test]
    fn detects_cation_pi() {
        // PHE CG, CD1, CD2, CE1, CE2, CZ centered around (0,0,0)
        // LYS NZ at (0, 0, 4.0) -> cation-pi interaction
        let pdb = "\
ATOM      1  CG  PHE A   1       0.000   0.000   0.000  1.00 0.00           C
ATOM      2  CD1 PHE A   1       1.000   0.000   0.000  1.00 0.00           C
ATOM      3  CD2 PHE A   1      -1.000   0.000   0.000  1.00 0.00           C
ATOM      4  CE1 PHE A   1       1.000   1.000   0.000  1.00 0.00           C
ATOM      5  CE2 PHE A   1      -1.000   1.000   0.000  1.00 0.00           C
ATOM      6  CZ  PHE A   1       0.000   2.000   0.000  1.00 0.00           C
ATOM      7  NZ  LYS A   2       0.000   0.666   4.000  1.00 0.00           N
";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let found = interactions(&s);
        assert_eq!(found.iter().filter(|i| i.kind == "cation_pi").count(), 1);
    }
}
