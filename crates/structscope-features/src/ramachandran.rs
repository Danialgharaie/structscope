//! Ramachandran favored / allowed / outlier classification from phi/psi (degrees).
//! Simplified Lovell et al. (2003)–style regions as embedded polygons.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RamachandranClass {
    Favored,
    Allowed,
    Outlier,
}

type Polygon = &'static [[f64; 2]; 4];

// General (non-Gly, non-Pro) — alpha-helix and beta-sheet favored cores.
const GENERAL_FAVORED: &[Polygon] = &[
    &[[-150.0, -12.0], [-33.0, -12.0], [-33.0, -78.0], [-150.0, -78.0]],
    &[[-180.0, 110.0], [-90.0, 110.0], [-90.0, 180.0], [-180.0, 180.0]],
    &[[90.0, 110.0], [180.0, 110.0], [180.0, 180.0], [90.0, 180.0]],
];

const GENERAL_ALLOWED: &[Polygon] = &[
    &[[-180.0, -180.0], [180.0, -180.0], [180.0, 180.0], [-180.0, 180.0]],
];

const GLY_FAVORED: &[Polygon] = &[
    &[[-150.0, -12.0], [-33.0, -12.0], [-33.0, -78.0], [-150.0, -78.0]],
    &[[-180.0, 30.0], [-60.0, 30.0], [-60.0, 180.0], [-180.0, 180.0]],
];

const GLY_ALLOWED: &[Polygon] = &[
    &[[-180.0, -180.0], [180.0, -180.0], [180.0, 80.0], [-180.0, 80.0]],
    &[[-180.0, 80.0], [-60.0, 80.0], [-60.0, 180.0], [-180.0, 180.0]],
];

const PRO_FAVORED: &[Polygon] = &[
    &[[-100.0, -30.0], [-30.0, -30.0], [-30.0, 180.0], [-100.0, 180.0]],
];

const PRO_ALLOWED: &[Polygon] = &[
    &[[-180.0, -180.0], [180.0, -180.0], [180.0, 180.0], [-180.0, 180.0]],
];

fn point_in_polygon(phi: f64, psi: f64, polygon: &[[f64; 2]]) -> bool {
    let n = polygon.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = (polygon[i][0], polygon[i][1]);
        let (xj, yj) = (polygon[j][0], polygon[j][1]);
        if ((yi > psi) != (yj > psi))
            && (phi < (xj - xi) * (psi - yi) / (yj - yi + f64::EPSILON) + xi)
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}

fn in_any(phi: f64, psi: f64, polygons: &[Polygon]) -> bool {
    polygons.iter().any(|poly| point_in_polygon(phi, psi, *poly))
}

fn regions(residue_name: &str) -> (&'static [Polygon], &'static [Polygon]) {
    match residue_name.trim().to_ascii_uppercase().as_str() {
        "GLY" => (GLY_FAVORED, GLY_ALLOWED),
        "PRO" => (PRO_FAVORED, PRO_ALLOWED),
        _ => (GENERAL_FAVORED, GENERAL_ALLOWED),
    }
}

pub fn classify(phi: f64, psi: f64, residue_name: &str) -> RamachandranClass {
    let (favored, allowed) = regions(residue_name);
    if in_any(phi, psi, favored) {
        RamachandranClass::Favored
    } else if in_any(phi, psi, allowed) {
        RamachandranClass::Allowed
    } else {
        RamachandranClass::Outlier
    }
}

pub fn class_label(class: RamachandranClass) -> &'static str {
    match class {
        RamachandranClass::Favored => "favored",
        RamachandranClass::Allowed => "allowed",
        RamachandranClass::Outlier => "outlier",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn helix_phi_psi_is_favored() {
        assert_eq!(classify(-60.0, -45.0, "ALA"), RamachandranClass::Favored);
    }

    #[test]
    fn gly_outside_general_region_is_outlier() {
        assert_eq!(classify(120.0, 120.0, "GLY"), RamachandranClass::Outlier);
    }

    #[test]
    fn proline_uses_pro_regions() {
        let c = classify(-60.0, 150.0, "PRO");
        assert_eq!(c, RamachandranClass::Favored);
    }
}
