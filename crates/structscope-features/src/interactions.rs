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
fn is_acidic_oxygen(res: &str, name: &str) -> bool {
    matches!((res, name), ("ASP", "OD1") | ("ASP", "OD2") | ("GLU", "OE1") | ("GLU", "OE2"))
}

/// True for sidechain basic nitrogens (LYS/ARG/HIS) — salt-bridge cations.
fn is_basic_nitrogen(res: &str, name: &str) -> bool {
    matches!(
        (res, name),
        ("LYS", "NZ") | ("ARG", "NH1") | ("ARG", "NH2") | ("ARG", "NE") | ("HIS", "ND1") | ("HIS", "NE2")
    )
}

/// Detect disulfides (CYS SG-SG < 2.5 A), salt bridges (acidic O / basic N < 4.0 A),
/// and polar contacts (N/O donor-acceptor 2.4-3.5 A). Returns one entry per pair.
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

            let salt = (is_acidic_oxygen(&a.res, &a.name) && is_basic_nitrogen(&b.res, &b.name))
                || (is_acidic_oxygen(&b.res, &b.name) && is_basic_nitrogen(&a.res, &a.name));
            if salt && d < 4.0 {
                out.push(Interaction { kind: "salt_bridge", atom_id_a: a.id.clone(), atom_id_b: b.id.clone(), distance: d });
                continue;
            }

            // Polar contact: a donor/acceptor pair of N/O atoms on different residues.
            let elem = |n: &str| n.chars().next().unwrap_or(' ');
            if a.res_id != b.res_id
                && matches!(elem(&a.name), 'N' | 'O')
                && matches!(elem(&b.name), 'N' | 'O')
                && (2.4..=3.5).contains(&d)
            {
                out.push(Interaction { kind: "hydrogen_bond", atom_id_a: a.id.clone(), atom_id_b: b.id.clone(), distance: d });
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
}
