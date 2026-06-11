use structscope_core::{Atom, Structure};

use crate::sasa::{sphere_points, vdw_radius, PROBE_RADIUS};

const N_SURFACE_POINTS: usize = 92;

#[derive(Debug, Clone, Copy)]
pub struct SurfaceDot {
    pub pos: [f64; 3],
    pub normal: [f64; 3],
}

fn is_heavy_atom(atom: &Atom) -> bool {
    if let Some(elem) = &atom.element {
        return !elem.trim().eq_ignore_ascii_case("H");
    }
    !atom.name.trim().starts_with('H')
}

fn weight_distance(d: f64) -> f64 {
    if d < 2.0 {
        1.0 - (d / 2.0).powi(2)
    } else {
        0.0
    }
}

fn dot_pair_score(_a: &SurfaceDot, _b: &SurfaceDot) -> f64 {
    let r = [
        _b.pos[0] - _a.pos[0],
        _b.pos[1] - _a.pos[1],
        _b.pos[2] - _a.pos[2],
    ];
    let d2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    if d2 <= f64::EPSILON {
        return 0.0;
    }

    let d = d2.sqrt();
    let u = [r[0] / d, r[1] / d, r[2] / d];
    // LC89 uses a product of distance and local geometric alignment terms.
    let wd = (-0.5 * d2).exp();
    let to_b = _a.normal[0] * u[0] + _a.normal[1] * u[1] + _a.normal[2] * u[2];
    let to_a = -(_b.normal[0] * u[0] + _b.normal[1] * u[1] + _b.normal[2] * u[2]);
    let wl = (to_b * to_a).max(0.0);
    wd * wl
}

pub fn lc89_score(dots_a: &[SurfaceDot], dots_b: &[SurfaceDot], search_radius: f64) -> f64 {
    if dots_a.is_empty() || dots_b.is_empty() {
        return 0.0;
    }

    let mut v = 0.0;
    let radius2 = search_radius * search_radius;
    for a in dots_a {
        for b in dots_b {
            let dx = a.pos[0] - b.pos[0];
            let dy = a.pos[1] - b.pos[1];
            let dz = a.pos[2] - b.pos[2];
            if dx * dx + dy * dy + dz * dz <= radius2 {
                v += dot_pair_score(a, b);
            }
        }
    }

    fn area_weight(dots: &[SurfaceDot]) -> f64 {
        let mut area = 0.0;
        for (i, a) in dots.iter().enumerate() {
            let mut local = 1.0;
            for (j, b) in dots.iter().enumerate() {
                if i == j {
                    continue;
                }
                let dx = a.pos[0] - b.pos[0];
                let dy = a.pos[1] - b.pos[1];
                let dz = a.pos[2] - b.pos[2];
                local += weight_distance((dx * dx + dy * dy + dz * dz).sqrt());
            }
            area += local;
        }
        area
    }

    let area = 0.5 * (area_weight(dots_a) + area_weight(dots_b));
    if area <= 0.0 {
        0.0
    } else {
        v / area
    }
}

pub fn interface_surface_dots(
    structure: &Structure,
    chain_label: &str,
    partner_label: &str,
    patch_distance: f64,
) -> Vec<SurfaceDot> {
    #[derive(Clone)]
    struct SphereAtom {
        center: [f64; 3],
        radius: f64,
        chain_label: String,
    }

    let mut atoms = Vec::new();
    for chain in &structure.chains {
        for residue in &chain.residues {
            for atom in &residue.atoms {
                if !is_heavy_atom(atom) {
                    continue;
                }
                atoms.push(SphereAtom {
                    center: [atom.x, atom.y, atom.z],
                    radius: vdw_radius(atom.element.as_deref().unwrap_or("")) + PROBE_RADIUS,
                    chain_label: chain.label.clone(),
                });
            }
        }
    }

    let chain_indices: Vec<usize> = atoms
        .iter()
        .enumerate()
        .filter_map(|(i, a)| (a.chain_label == chain_label).then_some(i))
        .collect();
    let partner_indices: Vec<usize> = atoms
        .iter()
        .enumerate()
        .filter_map(|(i, a)| (a.chain_label == partner_label).then_some(i))
        .collect();

    if chain_indices.is_empty() || partner_indices.is_empty() {
        return Vec::new();
    }

    let patch2 = patch_distance * patch_distance;
    let interface_indices: Vec<usize> = chain_indices
        .into_iter()
        .filter(|&i| {
            let a = &atoms[i];
            partner_indices.iter().any(|&j| {
                let b = &atoms[j];
                let dx = a.center[0] - b.center[0];
                let dy = a.center[1] - b.center[1];
                let dz = a.center[2] - b.center[2];
                dx * dx + dy * dy + dz * dz <= patch2
            })
        })
        .collect();

    let points = sphere_points(N_SURFACE_POINTS);
    let mut dots = Vec::new();
    for i in interface_indices {
        let atom = &atoms[i];
        let neighbours: Vec<&SphereAtom> = atoms
            .iter()
            .enumerate()
            .filter_map(|(j, other)| {
                if i == j {
                    return None;
                }
                let dx = other.center[0] - atom.center[0];
                let dy = other.center[1] - atom.center[1];
                let dz = other.center[2] - atom.center[2];
                let cutoff = atom.radius + other.radius;
                (dx * dx + dy * dy + dz * dz < cutoff * cutoff).then_some(other)
            })
            .collect();

        for p in &points {
            let pos = [
                atom.center[0] + p[0] * atom.radius,
                atom.center[1] + p[1] * atom.radius,
                atom.center[2] + p[2] * atom.radius,
            ];
            let accessible = neighbours.iter().all(|other| {
                let dx = pos[0] - other.center[0];
                let dy = pos[1] - other.center[1];
                let dz = pos[2] - other.center[2];
                dx * dx + dy * dy + dz * dz > other.radius * other.radius
            });
            if accessible {
                dots.push(SurfaceDot {
                    pos,
                    normal: *p,
                });
            }
        }
    }

    dots
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aligned_patches_have_higher_sc_than_misaligned() {
        let good = vec![
            SurfaceDot {
                pos: [0.0, 0.0, 0.0],
                normal: [0.0, 0.0, 1.0],
            },
            SurfaceDot {
                pos: [0.0, 0.0, 3.5],
                normal: [0.0, 0.0, -1.0],
            },
        ];
        let poor = vec![
            SurfaceDot {
                pos: [0.0, 0.0, 0.0],
                normal: [0.0, 0.0, 1.0],
            },
            SurfaceDot {
                pos: [0.0, 5.0, 3.5],
                normal: [0.0, 0.0, -1.0],
            },
        ];
        let sc_good = lc89_score(&[good[0]], &[good[1]], 5.0);
        let sc_poor = lc89_score(&[poor[0]], &[poor[1]], 5.0);
        assert!(sc_good > sc_poor);
        assert!(sc_good > 0.0);
    }
}
