use crate::dihedral::backbone_dihedrals;
use crate::ramachandran::{class_label, classify, RamachandranClass};
use crate::sasa::vdw_radius;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use structscope_core::{Atom, Residue, Structure};

const BACKBONE_ATOMS: [&str; 4] = ["N", "CA", "C", "O"];

const QUALITY_ALLOWLIST: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL", "MSE", "SEP", "TPO", "PTR", "HYP", "CSO", "CME", "KCX", "LLY", "MHO",
];

#[derive(Debug, Clone)]
pub struct QualityParams {
    pub clash_overlap: f64,
}

#[derive(Debug, Clone, Default)]
pub struct QualitySummary {
    pub quality_residue_count: usize,
    pub ramachandran_evaluated_count: usize,
    pub ramachandran_favored_count: usize,
    pub ramachandran_allowed_count: usize,
    pub ramachandran_outlier_count: usize,
    pub clash_pair_count: usize,
    pub missing_backbone_residue_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityRecord {
    pub structure_id: String,
    pub residue_id: String,
    pub chain_id: String,
    pub residue_name: String,
    pub seq_number: i32,
    pub phi: Option<f64>,
    pub psi: Option<f64>,
    pub ramachandran: Option<String>,
    pub missing_backbone: Vec<String>,
    pub clash_count: usize,
    pub clash_atom_ids: Vec<String>,
}

struct EvaluatedResidue<'a> {
    chain_label: &'a str,
    residue: &'a Residue,
    phi: Option<f64>,
    psi: Option<f64>,
    ramachandran: Option<RamachandranClass>,
    missing_backbone: Vec<String>,
}

struct HeavyAtom {
    id: String,
    residue_id: String,
    chain_label: String,
    residue_index: usize,
    x: f64,
    y: f64,
    z: f64,
    vdw: f64,
}

pub fn is_quality_residue(name: &str) -> bool {
    let n = name.trim().to_ascii_uppercase();
    QUALITY_ALLOWLIST.iter().any(|allowed| *allowed == n)
}

pub fn missing_backbone_atoms(residue: &Residue) -> Vec<String> {
    BACKBONE_ATOMS
        .iter()
        .filter(|name| !residue.atoms.iter().any(|a| a.name.trim() == **name))
        .map(|s| s.to_string())
        .collect()
}

fn is_heavy_atom(atom: &Atom) -> bool {
    if let Some(elem) = &atom.element {
        return !elem.trim().eq_ignore_ascii_case("H");
    }
    !atom.name.trim().starts_with('H')
}

fn collect_evaluated(structure: &Structure) -> Vec<EvaluatedResidue<'_>> {
    let dihedrals = backbone_dihedrals(structure);
    let mut dihedral_map: HashMap<(String, i32), (Option<f64>, Option<f64>)> = HashMap::new();
    for d in dihedrals {
        dihedral_map.insert((d.chain_id.clone(), d.seq_number), (d.phi, d.psi));
    }

    let mut out = Vec::new();
    for chain in &structure.chains {
        for (_residue_index, residue) in chain.residues.iter().enumerate() {
            if !is_quality_residue(&residue.name) {
                continue;
            }
            let missing = missing_backbone_atoms(residue);
            let (phi, psi) = dihedral_map
                .get(&(chain.id.clone(), residue.seq_number))
                .copied()
                .unwrap_or((None, None));
            let ramachandran = match (phi, psi) {
                (Some(p), Some(s)) if missing.is_empty() => Some(classify(p, s, &residue.name)),
                _ => None,
            };
            out.push(EvaluatedResidue {
                chain_label: &chain.label,
                residue,
                phi,
                psi,
                ramachandran,
                missing_backbone: missing,
            });
        }
    }
    out
}

fn collect_heavy_atoms(structure: &Structure) -> Vec<HeavyAtom> {
    let mut out = Vec::new();
    for chain in &structure.chains {
        for (residue_index, residue) in chain.residues.iter().enumerate() {
            if !is_quality_residue(&residue.name) {
                continue;
            }
            for atom in residue.atoms.iter().filter(|a| is_heavy_atom(a)) {
                out.push(HeavyAtom {
                    id: atom.id.clone(),
                    residue_id: residue.id.clone(),
                    chain_label: chain.label.clone(),
                    residue_index,
                    x: atom.x,
                    y: atom.y,
                    z: atom.z,
                    vdw: vdw_radius(atom.element.as_deref().unwrap_or("")),
                });
            }
        }
    }
    out
}

fn sequential_pair(chain_a: &str, idx_a: usize, chain_b: &str, idx_b: usize) -> bool {
    chain_a == chain_b && idx_a.abs_diff(idx_b) == 1
}

fn clash_pairs(atoms: &[HeavyAtom], overlap: f64) -> Vec<(usize, usize)> {
    if atoms.len() < 2 {
        return Vec::new();
    }
    let max_vdw = atoms.iter().map(|a| a.vdw).fold(0.0_f64, f64::max);
    let cell = (2.0 * max_vdw + overlap).max(1.0);
    let key = |a: &HeavyAtom| {
        (
            (a.x / cell).floor() as i64,
            (a.y / cell).floor() as i64,
            (a.z / cell).floor() as i64,
        )
    };
    let mut cells: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    for (i, atom) in atoms.iter().enumerate() {
        cells.entry(key(atom)).or_default().push(i);
    }

    let mut pairs = Vec::new();
    for (i, a) in atoms.iter().enumerate() {
        let (cx, cy, cz) = key(a);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(bucket) = cells.get(&(cx + dx, cy + dy, cz + dz)) else {
                        continue;
                    };
                    for &j in bucket {
                        if j <= i {
                            continue;
                        }
                        let b = &atoms[j];
                        if a.residue_id == b.residue_id {
                            continue;
                        }
                        if sequential_pair(&a.chain_label, a.residue_index, &b.chain_label, b.residue_index) {
                            continue;
                        }
                        let dist = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt();
                        let cutoff = a.vdw + b.vdw - overlap;
                        if dist < cutoff {
                            pairs.push((i, j));
                        }
                    }
                }
            }
        }
    }
    pairs
}

fn clash_counts_by_residue(atoms: &[HeavyAtom], pairs: &[(usize, usize)]) -> HashMap<String, (usize, Vec<String>)> {
    let mut map: HashMap<String, (usize, Vec<String>)> = HashMap::new();
    for &(i, j) in pairs {
        let (a, b) = (&atoms[i], &atoms[j]);
        for (res_id, partner) in [(&a.residue_id, b.id.clone()), (&b.residue_id, a.id.clone())] {
            let entry = map.entry(res_id.clone()).or_insert((0, Vec::new()));
            entry.0 += 1;
            if !entry.1.contains(&partner) {
                entry.1.push(partner);
            }
        }
    }
    map
}

fn is_problem(record: &QualityRecord) -> bool {
    record.ramachandran.as_deref() == Some("outlier")
        || !record.missing_backbone.is_empty()
        || record.clash_count > 0
}

fn build_records(
    structure: &Structure,
    evaluated: &[EvaluatedResidue<'_>],
    clash_map: &HashMap<String, (usize, Vec<String>)>,
) -> Vec<QualityRecord> {
    evaluated
        .iter()
        .map(|r| {
            let (clash_count, clash_atom_ids) = clash_map
                .get(&r.residue.id)
                .cloned()
                .unwrap_or((0, Vec::new()));
            QualityRecord {
                structure_id: structure.id.clone(),
                residue_id: r.residue.id.clone(),
                chain_id: r.chain_label.to_string(),
                residue_name: r.residue.name.clone(),
                seq_number: r.residue.seq_number,
                phi: r.phi,
                psi: r.psi,
                ramachandran: r.ramachandran.map(class_label).map(str::to_string),
                missing_backbone: r.missing_backbone.clone(),
                clash_count,
                clash_atom_ids,
            }
        })
        .collect()
}

pub fn quality_summary(structure: &Structure, params: &QualityParams) -> QualitySummary {
    let evaluated = collect_evaluated(structure);
    let atoms = collect_heavy_atoms(structure);
    let pairs = clash_pairs(&atoms, params.clash_overlap);

    let mut summary = QualitySummary {
        quality_residue_count: evaluated.len(),
        clash_pair_count: pairs.len(),
        missing_backbone_residue_count: evaluated.iter().filter(|r| !r.missing_backbone.is_empty()).count(),
        ..Default::default()
    };

    for r in &evaluated {
        if let Some(class) = r.ramachandran {
            summary.ramachandran_evaluated_count += 1;
            match class {
                RamachandranClass::Favored => summary.ramachandran_favored_count += 1,
                RamachandranClass::Allowed => summary.ramachandran_allowed_count += 1,
                RamachandranClass::Outlier => summary.ramachandran_outlier_count += 1,
            }
        }
    }

    summary
}

pub fn per_residue_quality(
    structure: &Structure,
    params: &QualityParams,
    all_residues: bool,
) -> Vec<QualityRecord> {
    let evaluated = collect_evaluated(structure);
    let atoms = collect_heavy_atoms(structure);
    let pairs = clash_pairs(&atoms, params.clash_overlap);
    let clash_map = clash_counts_by_residue(&atoms, &pairs);
    let records = build_records(structure, &evaluated, &clash_map);
    if all_residues {
        records
    } else {
        records.into_iter().filter(|r| is_problem(r)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use structscope_core::{parse_str, InputFormat, ParseOptions};

    #[test]
    fn detects_missing_backbone_o() {
        let pdb = "\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 0.00           C
";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let summary = quality_summary(&s, &QualityParams { clash_overlap: 0.4 });
        assert_eq!(summary.missing_backbone_residue_count, 1);
    }

    #[test]
    fn hoh_not_evaluated() {
        let pdb = "HETATM    1  O   HOH A   2      10.000   0.000   0.000  1.00 0.00           O\n";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let summary = quality_summary(&s, &QualityParams { clash_overlap: 0.4 });
        assert_eq!(summary.quality_residue_count, 0);
    }

    #[test]
    fn detects_vdw_clash_between_distant_residues() {
        let pdb = "\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 0.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00 0.00           C
ATOM      3  C   ALA A   1       1.500   0.000   0.000  1.00 0.00           C
ATOM      4  O   ALA A   1       2.000   1.000   0.000  1.00 0.00           O
ATOM      5  N   ALA A   2      10.000   0.000   0.000  1.00 0.00           N
ATOM      6  CA  ALA A   2      11.458   0.000   0.000  1.00 0.00           C
ATOM      7  C   ALA A   2      12.009   1.420   0.000  1.00 0.00           C
ATOM      8  O   ALA A   2      13.000   1.000   0.000  1.00 0.00           O
ATOM      9  N   ALA A   3       2.000   0.000   0.000  1.00 0.00           N
ATOM     10  CA  ALA A   3       2.000   0.000   0.000  1.00 0.00           C
ATOM     11  C   ALA A   3       3.500   0.000   0.000  1.00 0.00           C
ATOM     12  O   ALA A   3       4.000   1.000   0.000  1.00 0.00           O
";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let summary = quality_summary(&s, &QualityParams { clash_overlap: 0.4 });
        assert!(summary.clash_pair_count >= 1);
    }

    #[test]
    fn sequential_neighbors_not_clashed() {
        let pdb = "\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 0.00           C
ATOM      4  O   ALA A   1       2.500   2.200   0.000  1.00 0.00           O
ATOM      5  N   ALA A   2       3.332   1.540   0.000  1.00 0.00           N
ATOM      6  CA  ALA A   2       3.970   2.840   0.000  1.00 0.00           C
ATOM      7  C   ALA A   2       5.480   2.700   0.000  1.00 0.00           C
ATOM      8  O   ALA A   2       6.000   3.600   0.000  1.00 0.00           O
";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let summary = quality_summary(&s, &QualityParams { clash_overlap: 0.4 });
        assert_eq!(summary.clash_pair_count, 0);
    }

    #[test]
    fn problems_only_omits_clean_residue() {
        let pdb = "\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 0.00           C
ATOM      4  O   ALA A   1       2.500   2.200   0.000  1.00 0.00           O
ATOM      5  N   ALA A   2       3.332   1.540   0.000  1.00 0.00           N
ATOM      6  CA  ALA A   2       3.970   2.840   0.000  1.00 0.00           C
ATOM      7  C   ALA A   2       5.480   2.700   0.000  1.00 0.00           C
ATOM      8  O   ALA A   2       6.000   3.600   0.000  1.00 0.00           O
";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let params = QualityParams { clash_overlap: 0.4 };
        assert!(per_residue_quality(&s, &params, false).is_empty());
        assert_eq!(per_residue_quality(&s, &params, true).len(), 2);
    }

    #[test]
    fn all_residues_emits_every_evaluated() {
        let pdb = "\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 0.00           C
";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let records = per_residue_quality(&s, &QualityParams { clash_overlap: 0.4 }, true);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].missing_backbone, vec!["O".to_string()]);
    }
}
