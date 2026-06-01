//! Optimal rigid-body superposition via the quaternion method (Horn 1987; QCP).
//! Returns raw RMSD and the transform; the caller decides how to use them.

/// Result of superposing a mobile point set onto a reference set.
/// `aligned = rotation * mobile + translation`.
pub struct Superposition {
    pub rmsd: f64,
    pub rotation: [[f64; 3]; 3],
    pub translation: [f64; 3],
}

fn centroid(p: &[[f64; 3]]) -> [f64; 3] {
    let n = p.len() as f64;
    let mut c = [0.0; 3];
    for v in p {
        c[0] += v[0];
        c[1] += v[1];
        c[2] += v[2];
    }
    [c[0] / n, c[1] / n, c[2] / n]
}

/// Eigenvalues/eigenvectors of a symmetric 4x4 matrix via cyclic Jacobi rotations.
/// Eigenvector j is column j of the returned matrix.
fn jacobi_eigen(mut a: [[f64; 4]; 4]) -> ([f64; 4], [[f64; 4]; 4]) {
    let mut v = [[0.0; 4]; 4];
    for i in 0..4 {
        v[i][i] = 1.0;
    }
    for _ in 0..100 {
        let off = a[0][1].abs() + a[0][2].abs() + a[0][3].abs() + a[1][2].abs() + a[1][3].abs() + a[2][3].abs();
        if off < 1e-15 {
            break;
        }
        for p in 0..3 {
            for q in (p + 1)..4 {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for i in 0..4 {
                    let aip = a[i][p];
                    let aiq = a[i][q];
                    a[i][p] = c * aip - s * aiq;
                    a[i][q] = s * aip + c * aiq;
                }
                for i in 0..4 {
                    let api = a[p][i];
                    let aqi = a[q][i];
                    a[p][i] = c * api - s * aqi;
                    a[q][i] = s * api + c * aqi;
                }
                for i in 0..4 {
                    let vip = v[i][p];
                    let viq = v[i][q];
                    v[i][p] = c * vip - s * viq;
                    v[i][q] = s * vip + c * viq;
                }
            }
        }
    }
    ([a[0][0], a[1][1], a[2][2], a[3][3]], v)
}

/// Superpose `mobile` onto `reference` (equal length, point-to-point correspondence).
pub fn kabsch(mobile: &[[f64; 3]], reference: &[[f64; 3]]) -> Option<Superposition> {
    let n = mobile.len();
    if n == 0 || n != reference.len() {
        return None;
    }
    let cm = centroid(mobile);
    let cr = centroid(reference);

    // Correlation matrix s[i][j] = sum_k mobile_i * reference_j (centered), and E0.
    let mut s = [[0.0; 3]; 3];
    let mut e0 = 0.0;
    for k in 0..n {
        let b = [mobile[k][0] - cm[0], mobile[k][1] - cm[1], mobile[k][2] - cm[2]];
        let a = [reference[k][0] - cr[0], reference[k][1] - cr[1], reference[k][2] - cr[2]];
        e0 += a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
        for i in 0..3 {
            for j in 0..3 {
                s[i][j] += b[i] * a[j];
            }
        }
    }

    let (sxx, sxy, sxz) = (s[0][0], s[0][1], s[0][2]);
    let (syx, syy, syz) = (s[1][0], s[1][1], s[1][2]);
    let (szx, szy, szz) = (s[2][0], s[2][1], s[2][2]);

    // Horn key matrix.
    let k = [
        [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ];

    let (vals, vecs) = jacobi_eigen(k);
    let mut best = 0;
    for i in 1..4 {
        if vals[i] > vals[best] {
            best = i;
        }
    }
    let lambda = vals[best];
    let q = [vecs[0][best], vecs[1][best], vecs[2][best], vecs[3][best]];

    let rmsd = ((e0 - 2.0 * lambda).max(0.0) / n as f64).sqrt();

    // Quaternion (w,x,y,z) -> rotation mapping mobile onto reference. S was built as
    // mobile*reference^T (transpose of Horn's convention), so use the transposed matrix.
    let (w, x, y, z) = (q[0], q[1], q[2], q[3]);
    let rotation = [
        [w * w + x * x - y * y - z * z, 2.0 * (x * y - w * z), 2.0 * (z * x + w * y)],
        [2.0 * (x * y + w * z), w * w - x * x + y * y - z * z, 2.0 * (y * z - w * x)],
        [2.0 * (z * x - w * y), 2.0 * (y * z + w * x), w * w - x * x - y * y + z * z],
    ];

    let translation = [
        cr[0] - (rotation[0][0] * cm[0] + rotation[0][1] * cm[1] + rotation[0][2] * cm[2]),
        cr[1] - (rotation[1][0] * cm[0] + rotation[1][1] * cm[1] + rotation[1][2] * cm[2]),
        cr[2] - (rotation[2][0] * cm[0] + rotation[2][1] * cm[1] + rotation[2][2] * cm[2]),
    ];

    Some(Superposition { rmsd, rotation, translation })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn apply(sp: &Superposition, p: [f64; 3]) -> [f64; 3] {
        [
            sp.rotation[0][0] * p[0] + sp.rotation[0][1] * p[1] + sp.rotation[0][2] * p[2] + sp.translation[0],
            sp.rotation[1][0] * p[0] + sp.rotation[1][1] * p[1] + sp.rotation[1][2] * p[2] + sp.translation[1],
            sp.rotation[2][0] * p[0] + sp.rotation[2][1] * p[1] + sp.rotation[2][2] * p[2] + sp.translation[2],
        ]
    }

    #[test]
    fn identical_sets_zero_rmsd() {
        let p = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        let sp = kabsch(&p, &p).unwrap();
        assert!(sp.rmsd < 1e-9, "rmsd {}", sp.rmsd);
    }

    #[test]
    fn recovers_known_rotation() {
        // reference points; mobile = reference rotated 90 deg about z then translated.
        let reference = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 0.0], [2.0, 0.0, 1.0]];
        let mobile: Vec<[f64; 3]> = reference.iter().map(|p| [-p[1] + 5.0, p[0] - 3.0, p[2] + 1.0]).collect();
        let sp = kabsch(&mobile, &reference).unwrap();
        assert!(sp.rmsd < 1e-6, "rmsd {}", sp.rmsd);
        // Applying the transform maps mobile back onto reference.
        for (m, r) in mobile.iter().zip(reference.iter()) {
            let a = apply(&sp, *m);
            for d in 0..3 {
                assert!((a[d] - r[d]).abs() < 1e-6, "coord {d}: {} vs {}", a[d], r[d]);
            }
        }
    }

    #[test]
    fn noncongruent_sets_positive_rmsd() {
        let reference = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mobile = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let sp = kabsch(&mobile, &reference).unwrap();
        assert!(sp.rmsd > 0.1, "rmsd {}", sp.rmsd);
    }
}
