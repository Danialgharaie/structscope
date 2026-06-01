//! DSSP-style secondary structure assignment (Kabsch & Sander 1983).
//! Emits a raw per-residue state string (H/G/I helix, E strand, C coil) per chain.
use std::collections::HashSet;
use structscope_core::{Residue, Structure};

/// Per-chain secondary structure: `ss` has one character per residue, in order.
pub struct ChainSecondaryStructure {
    pub chain_id: String,
    pub ss: String,
}

fn atom(res: &Residue, name: &str) -> Option<[f64; 3]> {
    res.atoms.iter().find(|a| a.name == name).map(|a| [a.x, a.y, a.z])
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let (dx, dy, dz) = (a[0] - b[0], a[1] - b[1], a[2] - b[2]);
    (dx * dx + dy * dy + dz * dz).sqrt()
}

struct R {
    chain: usize,
    ca: Option<[f64; 3]>,
    c: Option<[f64; 3]>,
    o: Option<[f64; 3]>,
    n: Option<[f64; 3]>,
    h: Option<[f64; 3]>,
    is_pro: bool,
}

fn rank(c: u8) -> u8 {
    match c {
        b'H' => 5,
        b'E' => 4,
        b'G' => 3,
        b'I' => 2,
        _ => 1,
    }
}

fn set(ss: &mut [u8], i: usize, c: u8) {
    if rank(c) > rank(ss[i]) {
        ss[i] = c;
    }
}

/// Assign secondary structure for every chain in the structure.
pub fn secondary_structure(structure: &Structure) -> Vec<ChainSecondaryStructure> {
    // Flatten residues across chains, preserving order; one record per residue.
    let mut rs: Vec<R> = Vec::new();
    for (ci, chain) in structure.chains.iter().enumerate() {
        for res in &chain.residues {
            rs.push(R {
                chain: ci,
                ca: atom(res, "CA"),
                c: atom(res, "C"),
                o: atom(res, "O"),
                n: atom(res, "N"),
                h: None,
                is_pro: res.name.trim() == "PRO",
            });
        }
    }
    let len = rs.len();

    // Place amide H 1.0 A from N along the previous residue's O->C direction.
    for k in 1..len {
        if rs[k].chain != rs[k - 1].chain {
            continue;
        }
        if let (Some(n), Some(c), Some(o)) = (rs[k].n, rs[k - 1].c, rs[k - 1].o) {
            let d = [c[0] - o[0], c[1] - o[1], c[2] - o[2]];
            let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            if l > 1e-6 {
                rs[k].h = Some([n[0] + d[0] / l, n[1] + d[1] / l, n[2] + d[2] / l]);
            }
        }
    }

    // Kabsch-Sander H-bonds: (i,j) means CO(i) donates to NH(j). E<-0.5 kcal/mol.
    let mut hb: HashSet<(usize, usize)> = HashSet::new();
    for i in 0..len {
        let (c, o) = match (rs[i].c, rs[i].o) {
            (Some(c), Some(o)) => (c, o),
            _ => continue,
        };
        for j in 0..len {
            if i == j || rs[j].is_pro {
                continue;
            }
            let (n, h) = match (rs[j].n, rs[j].h) {
                (Some(n), Some(h)) => (n, h),
                _ => continue,
            };
            if let (Some(cai), Some(caj)) = (rs[i].ca, rs[j].ca) {
                if dist(cai, caj) > 9.0 {
                    continue;
                }
            }
            let (r_on, r_ch, r_oh, r_cn) = (dist(o, n), dist(c, h), dist(o, h), dist(c, n));
            if r_on < 1e-3 || r_ch < 1e-3 || r_oh < 1e-3 || r_cn < 1e-3 {
                continue;
            }
            let e = 27.888 * (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn);
            if e < -0.5 {
                hb.insert((i, j));
            }
        }
    }

    let bonded = |a: i64, b: i64| -> bool {
        a >= 0 && b >= 0 && (a as usize) < len && (b as usize) < len && hb.contains(&(a as usize, b as usize))
    };
    let cont = |a: i64, b: i64| -> bool {
        a >= 0 && b >= 0 && (a as usize) < len && (b as usize) < len && rs[a as usize].chain == rs[b as usize].chain
    };

    let mut ss = vec![b'C'; len];

    // Helices: two consecutive n-turns (bonds at i-1 and i) flag residues i..i+n-1.
    for n in [4usize, 3, 5] {
        let label = match n {
            4 => b'H',
            3 => b'G',
            _ => b'I',
        };
        for i in 1..len {
            let (i64, n64) = (i as i64, n as i64);
            if i + n < len
                && cont(i64 - 1, i64 + n64)
                && bonded(i64, i64 + n64)
                && bonded(i64 - 1, i64 - 1 + n64)
            {
                for k in i..(i + n) {
                    set(&mut ss, k, label);
                }
            }
        }
    }

    // Strands: residues participating in a parallel or antiparallel bridge -> E.
    for i in 0..len {
        let ii = i as i64;
        for j in (i + 1)..len {
            if rs[i].chain == rs[j].chain && j <= i + 2 {
                continue;
            }
            let jj = j as i64;
            let parallel = (cont(ii - 1, ii) && cont(ii, ii + 1) && bonded(ii - 1, jj) && bonded(jj, ii + 1))
                || (cont(jj - 1, jj) && cont(jj, jj + 1) && bonded(jj - 1, ii) && bonded(ii, jj + 1));
            let antiparallel = (bonded(ii, jj) && bonded(jj, ii))
                || (cont(ii - 1, ii)
                    && cont(ii, ii + 1)
                    && cont(jj - 1, jj)
                    && cont(jj, jj + 1)
                    && bonded(ii - 1, jj + 1)
                    && bonded(jj - 1, ii + 1));
            if parallel || antiparallel {
                set(&mut ss, i, b'E');
                set(&mut ss, j, b'E');
            }
        }
    }

    // Re-split the flat state array back into per-chain strings.
    let mut out = Vec::new();
    let mut idx = 0;
    for chain in &structure.chains {
        let mut s = String::with_capacity(chain.residues.len());
        for _ in &chain.residues {
            s.push(ss[idx] as char);
            idx += 1;
        }
        out.push(ChainSecondaryStructure { chain_id: chain.id.clone(), ss: s });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use structscope_core::{parse_str, InputFormat, ParseOptions};

    #[test]
    fn extended_tripeptide_is_all_coil() {
        // Three residues laid out in an extended line cannot form backbone H-bonds.
        let pdb = "\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 0.00           C
ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00 0.00           O
ATOM      5  N   ALA A   2       3.332   1.540   0.000  1.00 0.00           N
ATOM      6  CA  ALA A   2       3.970   2.840   0.000  1.00 0.00           C
ATOM      7  C   ALA A   2       5.480   2.700   0.000  1.00 0.00           C
ATOM      8  O   ALA A   2       6.000   1.590   0.000  1.00 0.00           O
ATOM      9  N   ALA A   3       6.150   3.830   0.000  1.00 0.00           N
ATOM     10  CA  ALA A   3       7.600   3.870   0.000  1.00 0.00           C
ATOM     11  C   ALA A   3       8.130   5.290   0.000  1.00 0.00           C
ATOM     12  O   ALA A   3       7.380   6.270   0.000  1.00 0.00           O
";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let ss = secondary_structure(&s);
        assert_eq!(ss[0].ss, "CCC");
    }
}
