use crate::model::{Atom, Chain, Residue, Structure, StructureMetadata};
use anyhow::Context;
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Clone, Copy)]
pub struct ParseOptions {
    pub keep_hetero: bool,
}

impl Default for ParseOptions {
    fn default() -> Self {
        Self { keep_hetero: true }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum InputFormat {
    Pdb,
    Mmcif,
}

#[derive(Debug, Error)]
pub enum ParseError {
    #[error("unsupported input format for {0}")]
    UnsupportedFormat(String),
    #[error("failed to read input: {0}")]
    Io(#[from] std::io::Error),
    #[error("failed to parse structure: {0}")]
    Invalid(String),
}

pub fn parse_file(path: &Path, options: ParseOptions) -> Result<Structure, ParseError> {
    let contents = fs::read_to_string(path)?;
    let format = detect_format(path)?;
    let source_path = Some(path.display().to_string());
    parse_str(&contents, format, source_path, options)
}

pub fn parse_str(
    contents: &str,
    format: InputFormat,
    source_path: Option<String>,
    options: ParseOptions,
) -> Result<Structure, ParseError> {
    match format {
        InputFormat::Pdb => parse_pdb(contents, source_path, options),
        InputFormat::Mmcif => parse_mmcif(contents, source_path, options),
    }
}

fn detect_format(path: &Path) -> Result<InputFormat, ParseError> {
    match path.extension().and_then(|ext| ext.to_str()) {
        Some("pdb") => Ok(InputFormat::Pdb),
        Some("cif") | Some("mmcif") => Ok(InputFormat::Mmcif),
        _ => Err(ParseError::UnsupportedFormat(path.display().to_string())),
    }
}

fn parse_pdb(
    contents: &str,
    source_path: Option<String>,
    options: ParseOptions,
) -> Result<Structure, ParseError> {
    let mut title = None;
    let mut chains: BTreeMap<String, BTreeMap<(i32, String, bool), Vec<Atom>>> = BTreeMap::new();
    let mut residue_names: BTreeMap<(String, i32, String, bool), String> = BTreeMap::new();

    for line in contents.lines() {
        if line.starts_with("TITLE") && title.is_none() {
            title = Some(slice(line, 10, line.len()).trim().to_string());
        }

        if !(line.starts_with("ATOM") || line.starts_with("HETATM")) {
            continue;
        }

        let is_hetero = line.starts_with("HETATM");
        if is_hetero && !options.keep_hetero {
            continue;
        }

        let serial = slice(line, 6, 11).trim().parse::<i32>().ok();
        let atom_name = slice(line, 12, 16).trim().to_string();
        let residue_name = slice(line, 17, 20).trim().to_string();
        let chain_label = slice(line, 21, 22).trim().to_string();
        let residue_seq = slice(line, 22, 26).trim().parse::<i32>().unwrap_or_default();
        let insertion_code = slice(line, 26, 27).trim().to_string();
        let x = slice(line, 30, 38).trim().parse::<f64>().unwrap_or_default();
        let y = slice(line, 38, 46).trim().parse::<f64>().unwrap_or_default();
        let z = slice(line, 46, 54).trim().parse::<f64>().unwrap_or_default();
        let occupancy = slice(line, 54, 60).trim().parse::<f64>().ok();
        let element = Some(slice(line, 76, 78).trim().to_string()).filter(|s| !s.is_empty());

        let chain_key = if chain_label.is_empty() {
            "_".to_string()
        } else {
            chain_label
        };
        let insertion = insertion_code.clone();
        let atom = Atom {
            id: String::new(),
            serial,
            name: atom_name,
            element,
            x,
            y,
            z,
            occupancy,
        };
        chains
            .entry(chain_key.clone())
            .or_default()
            .entry((residue_seq, insertion.clone(), is_hetero))
            .or_default()
            .push(atom);
        residue_names.insert((chain_key, residue_seq, insertion, is_hetero), residue_name);
    }

    build_structure(
        title,
        "pdb",
        source_path,
        chains,
        residue_names,
    )
}

fn parse_mmcif(
    contents: &str,
    source_path: Option<String>,
    options: ParseOptions,
) -> Result<Structure, ParseError> {
    let mut title = None;
    let mut atom_headers: Vec<String> = Vec::new();
    let mut in_atom_loop = false;
    let mut reading_headers = false;
    let mut chains: BTreeMap<String, BTreeMap<(i32, String, bool), Vec<Atom>>> = BTreeMap::new();
    let mut residue_names: BTreeMap<(String, i32, String, bool), String> = BTreeMap::new();

    for raw in contents.lines() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        if line.starts_with("_struct.title") && title.is_none() {
            title = Some(line.trim_start_matches("_struct.title").trim().trim_matches('\'').trim_matches('"').to_string());
            continue;
        }

        if line == "loop_" {
            atom_headers.clear();
            in_atom_loop = false;
            reading_headers = true;
            continue;
        }

        if reading_headers && line.starts_with("_atom_site.") {
            in_atom_loop = true;
            atom_headers.push(line.to_string());
            continue;
        }

        if reading_headers && line.starts_with('_') {
            continue;
        }

        if in_atom_loop {
            if line.starts_with('_') || line == "loop_" {
                break;
            }

            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < atom_headers.len() || atom_headers.is_empty() {
                continue;
            }

            let get = |name: &str| -> Option<&str> {
                atom_headers
                    .iter()
                    .position(|header| header == name)
                    .and_then(|idx| fields.get(idx).copied())
            };

            let group = get("_atom_site.group_PDB").unwrap_or("ATOM");
            let is_hetero = group == "HETATM";
            if is_hetero && !options.keep_hetero {
                continue;
            }

            let chain_label = get("_atom_site.auth_asym_id")
                .or_else(|| get("_atom_site.label_asym_id"))
                .unwrap_or("_");
            let residue_name = get("_atom_site.auth_comp_id")
                .or_else(|| get("_atom_site.label_comp_id"))
                .unwrap_or("UNK");
            let residue_seq = get("_atom_site.auth_seq_id")
                .or_else(|| get("_atom_site.label_seq_id"))
                .unwrap_or("0")
                .parse::<i32>()
                .unwrap_or_default();
            let insertion = get("_atom_site.pdbx_PDB_ins_code")
                .filter(|value| *value != "." && *value != "?")
                .unwrap_or("")
                .to_string();
            let atom_name = get("_atom_site.label_atom_id")
                .or_else(|| get("_atom_site.auth_atom_id"))
                .unwrap_or("X");
            let serial = get("_atom_site.id").and_then(|value| value.parse::<i32>().ok());
            let x = get("_atom_site.Cartn_x")
                .context("missing Cartn_x")
                .map_err(|err| ParseError::Invalid(err.to_string()))?
                .parse::<f64>()
                .unwrap_or_default();
            let y = get("_atom_site.Cartn_y")
                .context("missing Cartn_y")
                .map_err(|err| ParseError::Invalid(err.to_string()))?
                .parse::<f64>()
                .unwrap_or_default();
            let z = get("_atom_site.Cartn_z")
                .context("missing Cartn_z")
                .map_err(|err| ParseError::Invalid(err.to_string()))?
                .parse::<f64>()
                .unwrap_or_default();
            let occupancy = get("_atom_site.occupancy").and_then(|value| value.parse::<f64>().ok());
            let element = get("_atom_site.type_symbol")
                .map(|value| value.to_string())
                .filter(|value| !value.is_empty());

            let chain_key = chain_label.to_string();
            let atom = Atom {
                id: String::new(),
                serial,
                name: atom_name.to_string(),
                element,
                x,
                y,
                z,
                occupancy,
            };
            chains
                .entry(chain_key.clone())
                .or_default()
                .entry((residue_seq, insertion.clone(), is_hetero))
                .or_default()
                .push(atom);
            residue_names.insert(
                (chain_key, residue_seq, insertion, is_hetero),
                residue_name.to_string(),
            );
        }
    }

    build_structure(
        title,
        "mmcif",
        source_path,
        chains,
        residue_names,
    )
}

fn build_structure(
    title: Option<String>,
    source_format: &str,
    source_path: Option<String>,
    chains_map: BTreeMap<String, BTreeMap<(i32, String, bool), Vec<Atom>>>,
    residue_names: BTreeMap<(String, i32, String, bool), String>,
) -> Result<Structure, ParseError> {
    let structure_id = source_path
        .as_deref()
        .and_then(|path| Path::new(path).file_stem())
        .and_then(|stem| stem.to_str())
        .unwrap_or("structure")
        .to_string();

    let mut chains = Vec::new();
    for (chain_label, residues_map) in chains_map {
        let mut residues = Vec::new();
        for ((seq_number, insertion_code, is_hetero), atoms) in residues_map {
            let residue_id = format!(
                "{}:{}:{}:{}",
                structure_id,
                chain_label,
                seq_number,
                if insertion_code.is_empty() {
                    "_"
                } else {
                    &insertion_code
                }
            );
            let residue_name = residue_names
                .get(&(chain_label.clone(), seq_number, insertion_code.clone(), is_hetero))
                .cloned()
                .unwrap_or_else(|| "UNK".to_string());
            let atoms = atoms
                .into_iter()
                .enumerate()
                .map(|(idx, mut atom)| {
                    atom.id = format!("{residue_id}:{idx}:{}", atom.name);
                    atom
                })
                .collect();
            residues.push(Residue {
                id: residue_id,
                name: residue_name,
                seq_number,
                insertion_code: if insertion_code.is_empty() {
                    None
                } else {
                    Some(insertion_code)
                },
                atoms,
                is_hetero,
            });
        }
        let chain_id = format!("{structure_id}:{chain_label}");
        chains.push(Chain {
            id: chain_id,
            label: chain_label,
            residues,
        });
    }

    if chains.is_empty() {
        return Err(ParseError::Invalid("no atoms found".to_string()));
    }

    Ok(Structure {
        id: structure_id,
        metadata: StructureMetadata {
            source_format: source_format.to_string(),
            source_path,
            title,
        },
        chains,
    })
}

fn slice(line: &str, start: usize, end: usize) -> &str {
    line.get(start..line.len().min(end)).unwrap_or("")
}

#[cfg(test)]
mod tests {
    use super::*;

    const PDB_SAMPLE: &str = "\
ATOM      1  N   GLY A   1      11.104  13.207   8.292  1.00 20.00           N
ATOM      2  CA  GLY A   1      12.000  12.500   8.000  1.00 20.00           C
ATOM      3  C   GLY A   2      13.100  12.800   8.900  1.00 20.00           C
HETATM    4  O   HOH A 101      15.000  10.000   6.000  1.00 20.00           O
";

    const MMCIF_SAMPLE: &str = "\
data_demo
_struct.title 'Demo structure'
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N N GLY A 1 11.104 13.207 8.292
ATOM 2 C CA GLY A 1 12.000 12.500 8.000
ATOM 3 C C GLY A 2 13.100 12.800 8.900
";

    #[test]
    fn parses_pdb_summary() {
        let structure = parse_str(PDB_SAMPLE, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let summary = structure.summary();
        assert_eq!(summary.chain_count, 1);
        assert_eq!(summary.residue_count, 3);
        assert_eq!(summary.atom_count, 4);
        assert_eq!(summary.ligand_count, 1);
    }

    #[test]
    fn parses_mmcif_summary() {
        let structure = parse_str(MMCIF_SAMPLE, InputFormat::Mmcif, None, ParseOptions::default()).unwrap();
        let summary = structure.summary();
        assert_eq!(summary.chain_count, 1);
        assert_eq!(summary.residue_count, 2);
        assert_eq!(summary.atom_count, 3);
    }
}
