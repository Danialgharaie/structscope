use crate::interactions::{is_hydrogen_bond, is_hydrophobic_carbon, is_salt_bridge};
use crate::sasa::atom_sasa;
use serde::{Deserialize, Serialize};
use serde_json::{json, Map};
use structscope_core::{Atom, LigandFilter, Residue, Structure};

#[derive(Debug, Clone, Default)]
pub struct ProteinLigandSummary {
    pub ligand_count: usize,
    pub ligand_sasa_total: f64,
    pub ligand_sasa_mean: f64,
    pub binding_site_residue_count: usize,
    pub protein_ligand_hbond_count: usize,
    pub protein_ligand_salt_bridge_count: usize,
    pub protein_ligand_hydrophobic_count: usize,
    pub protein_ligand_contact_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LigandFeature {
    pub structure_id: String,
    pub ligand_id: String,
    pub residue_name: String,
    pub chain_id: String,
    pub seq_number: i32,
    pub sasa: f64,
    pub binding_site_residues: Vec<String>,
    pub interactions: Map<String, serde_json::Value>,
}

struct TaggedAtom {
    res_id: String,
    res_name: String,
    atom_name: String,
    x: f64,
    y: f64,
    z: f64,
    is_heavy: bool,
    is_ligand_carbon: bool,
    is_ligand: bool,
    is_protein: bool,
}

struct LigandContext {
    ligand_residues: Vec<(String, String, String, i32, f64)>,
    tagged_atoms: Vec<TaggedAtom>,
}

fn is_heavy_atom(atom: &Atom) -> bool {
    if let Some(elem) = &atom.element {
        return !elem.trim().eq_ignore_ascii_case("H");
    }
    !atom.name.trim().starts_with('H')
}

fn is_ligand_carbon(atom: &Atom) -> bool {
    if let Some(elem) = &atom.element {
        return elem.trim().eq_ignore_ascii_case("C");
    }
    atom.name.trim().starts_with('C')
}

fn is_protein_residue(residue: &Residue) -> bool {
    !residue.is_hetero
}

fn dist(a: &TaggedAtom, b: &TaggedAtom) -> f64 {
    ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt()
}

fn build_context(structure: &Structure, filter: &LigandFilter) -> LigandContext {
    let sasa = atom_sasa(structure);
    let mut atom_i = 0;
    let mut ligand_residues = Vec::new();
    let mut tagged_atoms = Vec::new();

    for chain in &structure.chains {
        for residue in &chain.residues {
            let is_lig = filter.is_ligand(residue);
            let is_prot = is_protein_residue(residue);
            let mut res_sasa = 0.0;

            for atom in &residue.atoms {
                if is_lig {
                    res_sasa += sasa[atom_i];
                }
                tagged_atoms.push(TaggedAtom {
                    res_id: residue.id.clone(),
                    res_name: residue.name.trim().to_string(),
                    atom_name: atom.name.trim().to_string(),
                    x: atom.x,
                    y: atom.y,
                    z: atom.z,
                    is_heavy: is_heavy_atom(atom),
                    is_ligand_carbon: is_ligand_carbon(atom),
                    is_ligand: is_lig,
                    is_protein: is_prot,
                });
                atom_i += 1;
            }

            if is_lig {
                ligand_residues.push((
                    residue.id.clone(),
                    chain.label.clone(),
                    residue.name.clone(),
                    residue.seq_number,
                    res_sasa,
                ));
            }
        }
    }

    LigandContext {
        ligand_residues,
        tagged_atoms,
    }
}

fn binding_site_residues(atoms: &[TaggedAtom], binding_distance: f64) -> Vec<String> {
    let ligand_heavy: Vec<&TaggedAtom> = atoms.iter().filter(|a| a.is_ligand && a.is_heavy).collect();
    let protein_heavy: Vec<&TaggedAtom> = atoms.iter().filter(|a| a.is_protein && a.is_heavy).collect();

    let mut sites = Vec::new();
    for protein in &protein_heavy {
        if sites.iter().any(|id: &String| id == &protein.res_id) {
            continue;
        }
        for ligand in &ligand_heavy {
            if dist(protein, ligand) <= binding_distance {
                sites.push(protein.res_id.clone());
                break;
            }
        }
    }
    sites.sort();
    sites
}

fn binding_site_for_ligand(atoms: &[TaggedAtom], ligand_res_id: &str, binding_distance: f64) -> Vec<String> {
    let ligand_heavy: Vec<&TaggedAtom> = atoms
        .iter()
        .filter(|a| a.is_ligand && a.is_heavy && a.res_id == ligand_res_id)
        .collect();
    let protein_heavy: Vec<&TaggedAtom> = atoms.iter().filter(|a| a.is_protein && a.is_heavy).collect();

    let mut sites = Vec::new();
    for protein in &protein_heavy {
        if sites.iter().any(|id: &String| id == &protein.res_id) {
            continue;
        }
        for ligand in &ligand_heavy {
            if dist(protein, ligand) <= binding_distance {
                sites.push(protein.res_id.clone());
                break;
            }
        }
    }
    sites.sort();
    sites
}

struct InteractionCounts {
    hbond: usize,
    salt_bridge: usize,
    hydrophobic: usize,
    contact: usize,
}

fn count_interactions(atoms: &[TaggedAtom], ligand_res_id: Option<&str>) -> InteractionCounts {
    let mut counts = InteractionCounts {
        hbond: 0,
        salt_bridge: 0,
        hydrophobic: 0,
        contact: 0,
    };

    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            let (a, b) = (&atoms[i], &atoms[j]);
            if !(a.is_ligand ^ b.is_ligand) {
                continue;
            }
            if !a.is_heavy || !b.is_heavy {
                continue;
            }
            if let Some(res_id) = ligand_res_id {
                let touches_ligand = a.res_id == res_id || b.res_id == res_id;
                if !touches_ligand {
                    continue;
                }
            }

            let (lig, prot) = if a.is_ligand { (a, b) } else { (b, a) };
            let d = dist(a, b);

            if d < 4.5 {
                counts.contact += 1;
            }
            if is_hydrogen_bond(&a.atom_name, &b.atom_name, d) {
                counts.hbond += 1;
            }
            if is_salt_bridge(&a.res_name, &a.atom_name, &b.res_name, &b.atom_name, d) {
                counts.salt_bridge += 1;
            }
            if d < 4.5 && lig.is_ligand_carbon && is_hydrophobic_carbon(&prot.res_name, &prot.atom_name) {
                counts.hydrophobic += 1;
            }
        }
    }

    counts
}

pub fn protein_ligand_summary(structure: &Structure, filter: &LigandFilter, binding_distance: f64) -> ProteinLigandSummary {
    let ctx = build_context(structure, filter);
    let ligand_count = ctx.ligand_residues.len();

    if ligand_count == 0 {
        return ProteinLigandSummary::default();
    }

    let ligand_sasa_total: f64 = ctx.ligand_residues.iter().map(|(_, _, _, _, sasa)| sasa).sum();
    let ligand_sasa_mean = ligand_sasa_total / ligand_count as f64;
    let binding_site = binding_site_residues(&ctx.tagged_atoms, binding_distance);
    let interactions = count_interactions(&ctx.tagged_atoms, None);

    ProteinLigandSummary {
        ligand_count,
        ligand_sasa_total,
        ligand_sasa_mean,
        binding_site_residue_count: binding_site.len(),
        protein_ligand_hbond_count: interactions.hbond,
        protein_ligand_salt_bridge_count: interactions.salt_bridge,
        protein_ligand_hydrophobic_count: interactions.hydrophobic,
        protein_ligand_contact_count: interactions.contact,
    }
}

pub fn per_ligand_features(structure: &Structure, filter: &LigandFilter, binding_distance: f64) -> Vec<LigandFeature> {
    let ctx = build_context(structure, filter);

    ctx.ligand_residues
        .into_iter()
        .map(|(ligand_id, chain_id, residue_name, seq_number, sasa)| {
            let binding_site_residues = binding_site_for_ligand(&ctx.tagged_atoms, &ligand_id, binding_distance);
            let counts = count_interactions(&ctx.tagged_atoms, Some(&ligand_id));
            let mut interactions = Map::new();
            interactions.insert("hydrogen_bond".to_string(), json!(counts.hbond));
            interactions.insert("salt_bridge".to_string(), json!(counts.salt_bridge));
            interactions.insert("hydrophobic".to_string(), json!(counts.hydrophobic));
            interactions.insert("contact".to_string(), json!(counts.contact));

            LigandFeature {
                structure_id: structure.id.clone(),
                ligand_id,
                residue_name,
                chain_id,
                seq_number,
                sasa,
                binding_site_residues,
                interactions,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use structscope_core::{parse_str, InputFormat, LigandFilter, ParseOptions};

    #[test]
    fn detects_protein_ligand_hbond() {
        let pdb = "\
ATOM      1  OD1 ASP A   1       0.000   0.000   0.000  1.00 0.00           O
HETATM    2  N1  HEM A 501       2.800   0.000   0.000  1.00 0.00           N
";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let filter = LigandFilter::default();
        let summary = protein_ligand_summary(&s, &filter, 5.0);
        assert_eq!(summary.protein_ligand_hbond_count, 1);
    }

    #[test]
    fn binding_site_counts_nearby_protein_residue() {
        let pdb = "\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 0.00           C
ATOM      2  CA  GLY A   2      10.000   0.000   0.000  1.00 0.00           C
HETATM    3  C1  HEM A 501       3.500   0.000   0.000  1.00 0.00           C
";
        let s = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let filter = LigandFilter::default();
        let summary = protein_ligand_summary(&s, &filter, 5.0);
        assert_eq!(summary.binding_site_residue_count, 1);
    }
}
