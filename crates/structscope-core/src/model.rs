use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Structure {
    pub id: StructureId,
    pub metadata: StructureMetadata,
    pub chains: Vec<Chain>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StructureMetadata {
    pub source_format: String,
    pub source_path: Option<String>,
    pub title: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Chain {
    pub id: ChainId,
    pub label: String,
    pub residues: Vec<Residue>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Residue {
    pub id: ResidueId,
    pub name: String,
    pub seq_number: i32,
    pub insertion_code: Option<String>,
    pub atoms: Vec<Atom>,
    pub is_hetero: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    pub id: AtomId,
    pub serial: Option<i32>,
    pub name: String,
    pub element: Option<String>,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub occupancy: Option<f64>,
    pub temp_factor: Option<f64>,
}

pub type StructureId = String;
pub type ChainId = String;
pub type ResidueId = String;
pub type AtomId = String;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParseSummary {
    pub structure_id: StructureId,
    pub chain_count: usize,
    pub residue_count: usize,
    pub atom_count: usize,
    pub heteroatom_count: usize,
    pub ligand_count: usize,
}

impl Structure {
    pub fn summary(&self) -> ParseSummary {
        let chain_count = self.chains.len();
        let residue_count = self.chains.iter().map(|chain| chain.residues.len()).sum();
        let atom_count = self
            .chains
            .iter()
            .flat_map(|chain| &chain.residues)
            .map(|residue| residue.atoms.len())
            .sum();
        let ligand_count = self
            .chains
            .iter()
            .flat_map(|chain| &chain.residues)
            .filter(|residue| residue.is_hetero)
            .count();
        let heteroatom_count = self
            .chains
            .iter()
            .flat_map(|chain| &chain.residues)
            .filter(|residue| residue.is_hetero)
            .flat_map(|residue| &residue.atoms)
            .count();

        ParseSummary {
            structure_id: self.id.clone(),
            chain_count,
            residue_count,
            atom_count,
            heteroatom_count,
            ligand_count,
        }
    }
}
