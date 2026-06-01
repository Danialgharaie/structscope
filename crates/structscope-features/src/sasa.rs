//! Solvent accessible surface area via the Shrake-Rupley algorithm.
//! Emits raw per-atom areas (Angstrom^2); interpretation is the caller's.
use structscope_core::Structure;

const PROBE_RADIUS: f64 = 1.4;
const N_POINTS: usize = 92;

/// Bondi van der Waals radii (Angstrom); default covers unlisted elements.
fn vdw_radius(element: &str) -> f64 {
    match element.trim().to_ascii_uppercase().as_str() {
        "H" => 1.20,
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "S" => 1.80,
        "P" => 1.80,
        "F" => 1.47,
        "CL" => 1.75,
        "BR" => 1.85,
        "I" => 1.98,
        _ => 1.70,
    }
}

/// Evenly distributed unit-sphere points via the golden spiral.
fn sphere_points(n: usize) -> Vec<[f64; 3]> {
    let golden = std::f64::consts::PI * (1.0 + 5.0_f64.sqrt());
    (0..n)
        .map(|k| {
            let y = 1.0 - 2.0 * (k as f64 + 0.5) / n as f64;
            let r = (1.0 - y * y).max(0.0).sqrt();
            let theta = golden * (k as f64 + 0.5);
            [theta.cos() * r, y, theta.sin() * r]
        })
        .collect()
}

struct Sphere {
    x: f64,
    y: f64,
    z: f64,
    r: f64,
}

/// Per-atom SASA in atom-iteration order (chains -> residues -> atoms).
pub fn atom_sasa(structure: &Structure) -> Vec<f64> {
    let spheres: Vec<Sphere> = structure
        .chains
        .iter()
        .flat_map(|chain| &chain.residues)
        .flat_map(|residue| &residue.atoms)
        .map(|atom| Sphere {
            x: atom.x,
            y: atom.y,
            z: atom.z,
            r: vdw_radius(atom.element.as_deref().unwrap_or("")) + PROBE_RADIUS,
        })
        .collect();

    let points = sphere_points(N_POINTS);

    spheres
        .iter()
        .enumerate()
        .map(|(i, atom)| {
            // Neighbours whose expanded spheres can overlap this atom's surface.
            let neighbours: Vec<&Sphere> = spheres
                .iter()
                .enumerate()
                .filter(|(j, other)| {
                    *j != i && {
                        let (dx, dy, dz) = (other.x - atom.x, other.y - atom.y, other.z - atom.z);
                        let cutoff = atom.r + other.r;
                        dx * dx + dy * dy + dz * dz < cutoff * cutoff
                    }
                })
                .map(|(_, other)| other)
                .collect();

            let accessible = points
                .iter()
                .filter(|p| {
                    let (px, py, pz) = (atom.x + p[0] * atom.r, atom.y + p[1] * atom.r, atom.z + p[2] * atom.r);
                    neighbours.iter().all(|other| {
                        let (dx, dy, dz) = (px - other.x, py - other.y, pz - other.z);
                        dx * dx + dy * dy + dz * dz > other.r * other.r
                    })
                })
                .count();

            4.0 * std::f64::consts::PI * atom.r * atom.r * accessible as f64 / N_POINTS as f64
        })
        .collect()
}

/// Total molecular SASA (Angstrom^2).
pub fn total_sasa(structure: &Structure) -> f64 {
    atom_sasa(structure).iter().sum()
}

/// Maximum accessible surface area per residue type (Angstrom^2), theoretical
/// values from Tien et al. 2013. Used to normalise SASA into relative exposure.
/// Returns `None` for non-standard residues (so RSA is left undefined).
pub fn max_accessible_area(residue_name: &str) -> Option<f64> {
    let v = match residue_name.trim().to_ascii_uppercase().as_str() {
        "ALA" => 129.0,
        "ARG" => 274.0,
        "ASN" => 195.0,
        "ASP" => 193.0,
        "CYS" => 167.0,
        "GLU" => 223.0,
        "GLN" => 225.0,
        "GLY" => 104.0,
        "HIS" => 224.0,
        "ILE" => 197.0,
        "LEU" => 201.0,
        "LYS" => 236.0,
        "MET" => 224.0,
        "PHE" => 240.0,
        "PRO" => 159.0,
        "SER" => 155.0,
        "THR" => 172.0,
        "TRP" => 285.0,
        "TYR" => 263.0,
        "VAL" => 174.0,
        _ => return None,
    };
    Some(v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use structscope_core::{parse_str, InputFormat, ParseOptions};

    #[test]
    fn isolated_atom_equals_full_sphere() {
        // Single carbon: SASA must equal the full expanded-sphere area 4*pi*(r+probe)^2.
        let pdb = "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 0.00           C\n";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let expected = 4.0 * std::f64::consts::PI * (1.70 + PROBE_RADIUS).powi(2);
        assert!((total_sasa(&s) - expected).abs() < 1e-9, "got {}", total_sasa(&s));
    }

    #[test]
    fn buried_atom_has_less_area() {
        // A central atom surrounded by close neighbours must expose less than isolated.
        let pdb = "\
ATOM      1  C   GLY A   1       0.000   0.000   0.000  1.00 0.00           C
ATOM      2  C   GLY A   2       2.000   0.000   0.000  1.00 0.00           C
ATOM      3  C   GLY A   3      -2.000   0.000   0.000  1.00 0.00           C
ATOM      4  C   GLY A   4       0.000   2.000   0.000  1.00 0.00           C
ATOM      5  C   GLY A   5       0.000  -2.000   0.000  1.00 0.00           C
";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let per_atom = atom_sasa(&s);
        let isolated = 4.0 * std::f64::consts::PI * (1.70 + PROBE_RADIUS).powi(2);
        assert!(per_atom[0] < isolated, "central atom should be partly buried");
    }
}
