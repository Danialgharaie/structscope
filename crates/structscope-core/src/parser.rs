use crate::model::{Atom, Chain, Residue, Structure, StructureMetadata};
use pdbrust::{
    parse_gzip_structure_file, parse_mmcif_string, parse_pdb_string, parse_structure_file, Atom as PdAtom,
    PdbStructure,
};
use std::collections::{BTreeMap, HashMap};
use std::path::Path;
use thiserror::Error;

use crate::bcif;

/// Format-neutral atom record that both the pdbrust and BinaryCIF paths feed into.
#[derive(Debug, Clone)]
struct RawAtom {
    serial: i32,
    name: String,
    alt_loc: Option<char>,
    residue_name: String,
    chain_id: String,
    residue_seq: i32,
    x: f64,
    y: f64,
    z: f64,
    occupancy: f64,
    temp_factor: f64,
    element: String,
    ins_code: Option<char>,
    is_hetatm: bool,
}

impl From<PdAtom> for RawAtom {
    fn from(a: PdAtom) -> Self {
        Self {
            serial: a.serial,
            name: a.name,
            alt_loc: a.alt_loc,
            residue_name: a.residue_name,
            chain_id: a.chain_id,
            residue_seq: a.residue_seq,
            x: a.x,
            y: a.y,
            z: a.z,
            occupancy: a.occupancy,
            temp_factor: a.temp_factor,
            element: a.element,
            ins_code: a.ins_code,
            is_hetatm: a.is_hetatm,
        }
    }
}

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
    BinaryCif,
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
    let format = detect_format(path)?;
    let source_path = Some(path.display().to_string());
    if matches!(format, InputFormat::BinaryCif) {
        let raw = std::fs::read(path)?;
        let bytes = if is_gzip_path(path) { gunzip(&raw)? } else { raw };
        return parse_bcif_bytes(&bytes, source_path, options);
    }
    let parsed = if is_gzip_path(path) {
        parse_gzip_structure_file(path).map_err(|err| ParseError::Invalid(err.to_string()))?
    } else {
        parse_structure_file(path).map_err(|err| ParseError::Invalid(err.to_string()))?
    };
    convert_structure(parsed, format, source_path, options)
}

pub fn parse_str(
    contents: &str,
    format: InputFormat,
    source_path: Option<String>,
    options: ParseOptions,
) -> Result<Structure, ParseError> {
    let parsed = match format {
        InputFormat::Pdb => parse_pdb_string(contents).map_err(|err| ParseError::Invalid(err.to_string()))?,
        InputFormat::Mmcif => parse_mmcif_string(contents).map_err(|err| ParseError::Invalid(err.to_string()))?,
        InputFormat::BinaryCif => return parse_bcif_bytes(contents.as_bytes(), source_path, options),
    };
    convert_structure(parsed, format, source_path, options)
}

fn gunzip(bytes: &[u8]) -> Result<Vec<u8>, ParseError> {
    use std::io::Read;
    let mut out = Vec::new();
    flate2::read::GzDecoder::new(bytes)
        .read_to_end(&mut out)
        .map_err(|err| ParseError::Invalid(format!("gzip: {err}")))?;
    Ok(out)
}

/// Decode a BinaryCIF buffer's atom_site category into a Structure.
fn parse_bcif_bytes(bytes: &[u8], source_path: Option<String>, options: ParseOptions) -> Result<Structure, ParseError> {
    let categories = bcif::decode_first_block(bytes).map_err(|err| ParseError::Invalid(err.to_string()))?;
    let atom_site = categories
        .iter()
        .find(|c| c.name == "atom_site")
        .ok_or_else(|| ParseError::Invalid("no atom_site category in BinaryCIF".to_string()))?;
    let title = categories
        .iter()
        .find(|c| c.name == "struct")
        .and_then(|c| c.column("title"))
        .and_then(|col| col.as_str(0));

    let atoms = bcif_atoms(atom_site);
    build_structure(atoms, InputFormat::BinaryCif, source_path, title, options)
}

/// Map the atom_site columns onto format-neutral RawAtoms, preferring author-assigned fields.
fn bcif_atoms(cat: &bcif::Category) -> Vec<RawAtom> {
    let col = |name: &str| cat.column(name);
    let pick = |a: &str, b: &str| col(a).or_else(|| col(b));

    let group = col("group_PDB");
    let chain = pick("label_asym_id", "auth_asym_id");
    let seq = pick("auth_seq_id", "label_seq_id");
    let comp = pick("auth_comp_id", "label_comp_id");
    let atom = pick("auth_atom_id", "label_atom_id");
    let alt = col("label_alt_id");
    let ins = col("pdbx_PDB_ins_code");
    let (x, y, z) = (col("Cartn_x"), col("Cartn_y"), col("Cartn_z"));
    let occ = col("occupancy");
    let b_fac = col("B_iso_or_equiv");
    let sym = col("type_symbol");
    let id = col("id");

    (0..cat.row_count)
        .map(|i| {
            let one_char = |c: Option<&bcif::ColumnData>| {
                c.and_then(|d| d.as_str(i)).and_then(|s| s.trim().chars().next())
            };
            RawAtom {
                serial: id.and_then(|d| d.as_i64(i)).unwrap_or((i + 1) as i64) as i32,
                name: atom.and_then(|d| d.as_str(i)).unwrap_or_default(),
                alt_loc: one_char(alt),
                residue_name: comp.and_then(|d| d.as_str(i)).unwrap_or_default(),
                chain_id: chain.and_then(|d| d.as_str(i)).unwrap_or_default(),
                residue_seq: seq.and_then(|d| d.as_i64(i)).unwrap_or(0) as i32,
                x: x.and_then(|d| d.as_f64(i)).unwrap_or(0.0),
                y: y.and_then(|d| d.as_f64(i)).unwrap_or(0.0),
                z: z.and_then(|d| d.as_f64(i)).unwrap_or(0.0),
                occupancy: occ.and_then(|d| d.as_f64(i)).unwrap_or(1.0),
                temp_factor: b_fac.and_then(|d| d.as_f64(i)).unwrap_or(0.0),
                element: sym.and_then(|d| d.as_str(i)).unwrap_or_default(),
                ins_code: one_char(ins),
                is_hetatm: group.and_then(|d| d.as_str(i)).map(|g| g.trim() == "HETATM").unwrap_or(false),
            }
        })
        .collect()
}

fn detect_format(path: &Path) -> Result<InputFormat, ParseError> {
    let name = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| ParseError::UnsupportedFormat(path.display().to_string()))?
        .to_ascii_lowercase();

    if name.ends_with(".pdb") || name.ends_with(".ent") || name.ends_with(".pdb.gz") || name.ends_with(".ent.gz") {
        Ok(InputFormat::Pdb)
    } else if name.ends_with(".cif")
        || name.ends_with(".mmcif")
        || name.ends_with(".cif.gz")
        || name.ends_with(".mmcif.gz")
    {
        Ok(InputFormat::Mmcif)
    } else if name.ends_with(".bcif") || name.ends_with(".bcif.gz") {
        Ok(InputFormat::BinaryCif)
    } else {
        Err(ParseError::UnsupportedFormat(path.display().to_string()))
    }
}

fn convert_structure(
    parsed: PdbStructure,
    format: InputFormat,
    source_path: Option<String>,
    options: ParseOptions,
) -> Result<Structure, ParseError> {
    let title = parsed.title.or(parsed.header);
    let atoms = parsed.atoms.into_iter().map(RawAtom::from).collect();
    build_structure(atoms, format, source_path, title, options)
}

fn build_structure(
    atoms: Vec<RawAtom>,
    format: InputFormat,
    source_path: Option<String>,
    title: Option<String>,
    options: ParseOptions,
) -> Result<Structure, ParseError> {
    let source_format = match format {
        InputFormat::Pdb => "pdb",
        InputFormat::Mmcif => "mmcif",
        InputFormat::BinaryCif => "bcif",
    };
    let structure_id = source_path
        .as_deref()
        .map(structure_id_from_path)
        .unwrap_or_else(|| "structure".to_string());

    let filtered_atoms = retain_selected_altlocs(atoms, options.keep_hetero);
    let mut chains_map: BTreeMap<String, BTreeMap<(i32, String, bool), Vec<RawAtom>>> = BTreeMap::new();

    for atom in filtered_atoms {
        let chain_key = if atom.chain_id.trim().is_empty() {
            "_".to_string()
        } else {
            atom.chain_id.trim().to_string()
        };
        let insertion = atom.ins_code.map(|code| code.to_string()).unwrap_or_default();
        chains_map
            .entry(chain_key)
            .or_default()
            .entry((atom.residue_seq, insertion, atom.is_hetatm))
            .or_default()
            .push(atom);
    }

    let mut chains = Vec::new();
    for (chain_label, residues_map) in chains_map {
        let chain_id = format!("{structure_id}:{chain_label}");
        let mut residues = Vec::new();
        for ((seq_number, insertion_code, is_hetero), mut atoms) in residues_map {
            atoms.sort_by(|left, right| {
                left.serial
                    .cmp(&right.serial)
                    .then_with(|| left.name.cmp(&right.name))
                    .then_with(|| left.alt_loc.cmp(&right.alt_loc))
            });

            let residue_name = atoms
                .first()
                .map(|atom| atom.residue_name.trim().to_string())
                .unwrap_or_else(|| "UNK".to_string());
            let residue_id = format!(
                "{}:{}:{}:{}",
                structure_id,
                chain_label,
                seq_number,
                if insertion_code.is_empty() { "_" } else { &insertion_code }
            );

            let atoms = atoms
                .into_iter()
                .enumerate()
                .map(|(idx, atom)| Atom {
                    id: format!("{residue_id}:{idx}:{}", atom.name.trim()),
                    serial: Some(atom.serial),
                    name: atom.name.trim().to_string(),
                    element: normalize_string(Some(atom.element)),
                    x: atom.x,
                    y: atom.y,
                    z: atom.z,
                    occupancy: Some(atom.occupancy),
                    temp_factor: Some(atom.temp_factor),
                })
                .collect();

            residues.push(Residue {
                id: residue_id,
                name: residue_name,
                seq_number,
                insertion_code: normalize_string(Some(insertion_code)),
                atoms,
                is_hetero,
            });
        }

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

fn retain_selected_altlocs(atoms: Vec<RawAtom>, keep_hetero: bool) -> Vec<RawAtom> {
    let mut residue_altloc_scores: HashMap<(String, i32, Option<char>, bool), HashMap<char, (f64, usize)>> =
        HashMap::new();

    for atom in &atoms {
        if atom.is_hetatm && !keep_hetero {
            continue;
        }
        if let Some(alt_loc) = atom.alt_loc {
            let residue_key = (
                atom.chain_id.clone(),
                atom.residue_seq,
                atom.ins_code,
                atom.is_hetatm,
            );
            let entry = residue_altloc_scores.entry(residue_key).or_default();
            let score = entry.entry(alt_loc).or_insert((0.0, 0));
            score.0 += atom.occupancy;
            score.1 += 1;
        }
    }

    let selected_altlocs: HashMap<(String, i32, Option<char>, bool), char> = residue_altloc_scores
        .into_iter()
        .filter_map(|(residue_key, choices)| {
            choices
                .into_iter()
                .max_by(|left, right| {
                    left.1
                        .0
                        .partial_cmp(&right.1 .0)
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .then_with(|| left.1 .1.cmp(&right.1 .1))
                        .then_with(|| left.0.cmp(&right.0))
                })
                .map(|(alt_loc, _)| (residue_key, alt_loc))
        })
        .collect();

    atoms
        .into_iter()
        .filter(|atom| !atom.is_hetatm || keep_hetero)
        .filter(|atom| {
            let residue_key = (
                atom.chain_id.clone(),
                atom.residue_seq,
                atom.ins_code,
                atom.is_hetatm,
            );
            match atom.alt_loc {
                None => true,
                Some(alt_loc) => selected_altlocs
                    .get(&residue_key)
                    .map(|selected| *selected == alt_loc)
                    .unwrap_or(true),
            }
        })
        .collect()
}

fn structure_id_from_path(path: &str) -> String {
    let mut name = Path::new(path)
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("structure")
        .to_string();

    for suffix in [".pdb.gz", ".cif.gz", ".mmcif.gz", ".bcif.gz", ".pdb", ".ent", ".cif", ".mmcif", ".bcif"] {
        if name.to_ascii_lowercase().ends_with(suffix) {
            let new_len = name.len().saturating_sub(suffix.len());
            name.truncate(new_len);
            break;
        }
    }

    if name.is_empty() {
        "structure".to_string()
    } else {
        name
    }
}

fn normalize_string<T: Into<String>>(value: Option<T>) -> Option<String> {
    value
        .map(Into::into)
        .map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty())
}

fn is_gzip_path(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .map(|name| name.to_ascii_lowercase().ends_with(".gz"))
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::fs;
    use std::io::Write;
    use std::time::{SystemTime, UNIX_EPOCH};

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

    #[test]
    fn keeps_highest_occupancy_altloc_per_residue() {
        let pdb = "\
ATOM      1  CA AGLY A   1      11.104  13.207   8.292  0.40 20.00           C
ATOM      2  CA BGLY A   1      21.000  23.000  18.000  0.60 20.00           C
ATOM      3  N   GLY A   1      12.000  12.500   8.000  1.00 20.00           N
";
        let structure = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
        let residue = &structure.chains[0].residues[0];
        assert_eq!(residue.atoms.len(), 2);
        let ca = residue.atoms.iter().find(|atom| atom.name == "CA").unwrap();
        assert_eq!(ca.x, 21.0);
    }

    #[test]
    fn parses_gzip_pdb_file() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!("structscope-core-{unique}"));
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("demo.pdb.gz");

        let file = fs::File::create(&path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(PDB_SAMPLE.as_bytes()).unwrap();
        encoder.finish().unwrap();

        let structure = parse_file(&path, ParseOptions::default()).unwrap();
        assert_eq!(structure.summary().atom_count, 4);

        fs::remove_file(&path).unwrap();
        fs::remove_dir_all(&dir).unwrap();
    }
}
