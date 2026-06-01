//! Backbone dihedral angles (phi, psi, omega) per residue, in degrees.
//! Raw geometry from N/CA/C coordinates; interpretation is the caller's.
use structscope_core::{Residue, Structure};

/// Per-residue backbone torsions in degrees. `None` where an angle is undefined
/// (chain termini) or required backbone atoms are missing.
pub struct ResidueDihedrals {
    pub chain_id: String,
    pub seq_number: i32,
    pub phi: Option<f64>,
    pub psi: Option<f64>,
    pub omega: Option<f64>,
}

fn atom(res: &Residue, name: &str) -> Option<[f64; 3]> {
    res.atoms.iter().find(|a| a.name == name).map(|a| [a.x, a.y, a.z])
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Signed dihedral angle (degrees) defined by four points, IUPAC convention.
fn torsion(p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], p4: [f64; 3]) -> Option<f64> {
    let b1 = sub(p2, p1);
    let b2 = sub(p3, p2);
    let b3 = sub(p4, p3);
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let b2n = norm(b2);
    if b2n < 1e-9 || norm(n1) < 1e-9 || norm(n2) < 1e-9 {
        return None;
    }
    let b2hat = [b2[0] / b2n, b2[1] / b2n, b2[2] / b2n];
    let x = dot(n1, n2);
    let y = dot(cross(n1, n2), b2hat);
    Some(y.atan2(x).to_degrees())
}

/// Compute phi/psi/omega for every residue, chain by chain (termini neighbours within-chain).
pub fn backbone_dihedrals(structure: &Structure) -> Vec<ResidueDihedrals> {
    let mut out = Vec::new();
    for chain in &structure.chains {
        let r = &chain.residues;
        for i in 0..r.len() {
            let (n, ca, c) = (atom(&r[i], "N"), atom(&r[i], "CA"), atom(&r[i], "C"));
            let prev_c = if i > 0 { atom(&r[i - 1], "C") } else { None };
            let prev_ca = if i > 0 { atom(&r[i - 1], "CA") } else { None };
            let next_n = if i + 1 < r.len() { atom(&r[i + 1], "N") } else { None };

            let phi = match (prev_c, n, ca, c) {
                (Some(pc), Some(n), Some(ca), Some(c)) => torsion(pc, n, ca, c),
                _ => None,
            };
            let psi = match (n, ca, c, next_n) {
                (Some(n), Some(ca), Some(c), Some(nn)) => torsion(n, ca, c, nn),
                _ => None,
            };
            let omega = match (prev_ca, prev_c, n, ca) {
                (Some(pca), Some(pc), Some(n), Some(ca)) => torsion(pca, pc, n, ca),
                _ => None,
            };

            out.push(ResidueDihedrals {
                chain_id: chain.id.clone(),
                seq_number: r[i].seq_number,
                phi,
                psi,
                omega,
            });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn torsion_sign_matches_iupac() {
        // Hand-derived: p1=(0,0,1) p2=origin p3=(1,0,0) p4=(1,1,0) -> -90 deg (IUPAC).
        let a = torsion([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]).unwrap();
        assert!((a - -90.0).abs() < 1e-6, "got {a}");
        // Mirror across the b2 axis flips the sign to +90 deg.
        let b = torsion([0.0, 0.0, -1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]).unwrap();
        assert!((b - 90.0).abs() < 1e-6, "got {b}");
    }

    #[test]
    fn torsion_trans_is_180() {
        let a = torsion([0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, -1.0, 0.0]).unwrap();
        assert!((a.abs() - 180.0).abs() < 1e-6, "got {a}");
    }

    #[test]
    fn termini_angles_are_none() {
        use structscope_core::{parse_str, InputFormat, ParseOptions};
        let pdb = "\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 0.00           C
ATOM      4  N   ALA A   2       3.332   1.540   0.000  1.00 0.00           N
ATOM      5  CA  ALA A   2       3.970   2.840   0.000  1.00 0.00           C
ATOM      6  C   ALA A   2       5.480   2.700   0.000  1.00 0.00           C
";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let d = backbone_dihedrals(&s);
        assert_eq!(d.len(), 2);
        assert!(d[0].phi.is_none(), "first residue has no phi");
        assert!(d[1].psi.is_none(), "last residue has no psi");
        // Interior connectivity exists: residue 1 psi and residue 2 phi are defined.
        assert!(d[0].psi.is_some() && d[1].phi.is_some());
    }
}
