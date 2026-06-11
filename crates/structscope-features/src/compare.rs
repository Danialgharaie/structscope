use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value, json};
use structscope_core::{superposition_rmsd, RmsdParams, Structure};

use crate::FeatureRecord;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReferenceMode {
    FirstInput,
    ExplicitPath(PathBuf),
    AutoQuality,
    ByField { maximize: bool, field: String },
}

pub fn pick_reference_index(
    records: &[FeatureRecord],
    paths: &[PathBuf],
    mode: &ReferenceMode,
) -> Result<usize> {
    if records.is_empty() {
        bail!("no feature records");
    }
    if paths.len() != records.len() {
        bail!(
            "paths length ({}) must match records length ({})",
            paths.len(),
            records.len()
        );
    }

    match mode {
        ReferenceMode::FirstInput => Ok(0),
        ReferenceMode::ExplicitPath(path) => find_explicit_index(records, paths, path),
        ReferenceMode::AutoQuality => pick_auto_quality(records),
        ReferenceMode::ByField { maximize, field } => pick_by_field(records, *maximize, field),
    }
}

fn find_explicit_index(
    records: &[FeatureRecord],
    paths: &[PathBuf],
    target: &Path,
) -> Result<usize> {
    if let Some(index_str) = target.to_str().and_then(|s| s.strip_prefix('#')) {
        if let Ok(index) = index_str.parse::<usize>() {
            if index < records.len() {
                return Ok(index);
            }
            bail!("reference index {index} out of range (0..{})", records.len());
        }
    }

    let target_str = target.to_string_lossy();
    for (index, (record, path)) in records.iter().zip(paths.iter()).enumerate() {
        if paths_match(path, target)
            || record
                .source_path
                .as_deref()
                .is_some_and(|source| paths_match(Path::new(source), target))
            || record.structure_id == target_str
            || path.file_name().is_some_and(|name| name == target.as_os_str())
        {
            return Ok(index);
        }
    }

    bail!("reference path not found in input set: {}", target.display())
}

fn paths_match(left: &Path, right: &Path) -> bool {
    left == right
        || left
            .file_name()
            .zip(right.file_name())
            .is_some_and(|(left, right)| left == right)
        || canonicalize_if_exists(left).is_some_and(|left| {
            canonicalize_if_exists(right).is_some_and(|right| left == right)
        })
}

fn canonicalize_if_exists(path: &Path) -> Option<PathBuf> {
    path.canonicalize().ok()
}

fn pick_auto_quality(records: &[FeatureRecord]) -> Result<usize> {
    let mut best_index = 0usize;
    let mut best_key = quality_key(&records[0], 0)?;

    for (index, record) in records.iter().enumerate().skip(1) {
        let key = quality_key(record, index)?;
        if key < best_key {
            best_key = key;
            best_index = index;
        }
    }

    Ok(best_index)
}

fn quality_key(record: &FeatureRecord, tie_breaker: usize) -> Result<(u64, u64, u64, usize)> {
    Ok((
        feature_u64(record, "ramachandran_outlier_count")?,
        feature_u64(record, "clash_pair_count")?,
        feature_u64(record, "missing_backbone_residue_count")?,
        tie_breaker,
    ))
}

fn pick_by_field(records: &[FeatureRecord], maximize: bool, field: &str) -> Result<usize> {
    let mut best_index = 0usize;
    let mut best_value = feature_f64(&records[0], field)?;

    for (index, record) in records.iter().enumerate().skip(1) {
        let value = feature_f64(record, field)?;
        let is_better = if maximize {
            value > best_value
        } else {
            value < best_value
        };
        if is_better {
            best_index = index;
            best_value = value;
        }
    }

    Ok(best_index)
}

fn feature_u64(record: &FeatureRecord, field: &str) -> Result<u64> {
    let value = record
        .features
        .get(field)
        .ok_or_else(|| anyhow!("missing feature field: {field}"))?;
    value
        .as_u64()
        .or_else(|| value.as_f64().map(|v| v as u64))
        .ok_or_else(|| anyhow!("feature {field} is not numeric"))
}

fn feature_f64(record: &FeatureRecord, field: &str) -> Result<f64> {
    let value = record
        .features
        .get(field)
        .ok_or_else(|| anyhow!("missing feature field: {field}"))?;
    value_as_f64(value).ok_or_else(|| anyhow!("feature {field} is not numeric"))
}

fn value_as_f64(value: &Value) -> Option<f64> {
    match value {
        Value::Number(number) => number.as_f64(),
        _ => None,
    }
}

pub fn numeric_feature_keys(record: &FeatureRecord) -> Vec<String> {
    let mut keys: Vec<String> = record
        .features
        .iter()
        .filter(|(_, value)| value.is_number())
        .map(|(key, _)| key.clone())
        .collect();
    keys.sort();
    keys
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeltaRecord {
    pub structure_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_path: Option<String>,
    #[serde(flatten)]
    pub fields: Map<String, Value>,
}

pub fn feature_deltas(
    records: &[FeatureRecord],
    ref_idx: usize,
    fields: Option<&[String]>,
) -> Vec<DeltaRecord> {
    let ref_record = &records[ref_idx];
    let field_names: Vec<String> = match fields {
        Some(names) => names.to_vec(),
        None => {
            let mut keys: Vec<String> = records
                .iter()
                .flat_map(numeric_feature_keys)
                .collect();
            keys.sort();
            keys.dedup();
            keys
        }
    };

    records
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != ref_idx)
        .map(|(_, record)| {
            let mut delta_fields = Map::new();
            for field in &field_names {
                let Some(value) = record.features.get(field).and_then(value_as_f64) else {
                    continue;
                };
                let Some(ref_value) = ref_record.features.get(field).and_then(value_as_f64) else {
                    continue;
                };
                delta_fields.insert(field.clone(), json!(value));
                delta_fields.insert(format!("reference_{field}"), json!(ref_value));
                delta_fields.insert(format!("delta_{field}"), json!(value - ref_value));
            }
            DeltaRecord {
                structure_id: record.structure_id.clone(),
                source_path: record.source_path.clone(),
                fields: delta_fields,
            }
        })
        .collect()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RmsdMatrix {
    pub labels: Vec<String>,
    pub rmsd: Vec<Vec<Option<f64>>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParseFailure {
    pub path: String,
    pub error: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CompareResult {
    pub reference_id: String,
    pub reference_path: String,
    pub reference_mode: String,
    pub structure_count: usize,
    pub structures: Vec<String>,
    pub matrix: RmsdMatrix,
    pub deltas: Vec<DeltaRecord>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub failed: Option<Vec<ParseFailure>>,
}

pub fn reference_mode_label(mode: &ReferenceMode) -> String {
    match mode {
        ReferenceMode::FirstInput => "first_input".to_string(),
        ReferenceMode::ExplicitPath(_) => "explicit".to_string(),
        ReferenceMode::AutoQuality => "auto_quality".to_string(),
        ReferenceMode::ByField { .. } => "reference_by".to_string(),
    }
}

pub fn compare_set(
    structures: &[Structure],
    records: &[FeatureRecord],
    paths: &[PathBuf],
    rmsd_params: &RmsdParams,
    reference_mode: &ReferenceMode,
    delta_fields: Option<&[String]>,
    failed: Option<Vec<ParseFailure>>,
) -> Result<CompareResult> {
    if structures.len() != records.len() || structures.len() != paths.len() {
        bail!(
            "structures ({}), records ({}), and paths ({}) must have the same length",
            structures.len(),
            records.len(),
            paths.len()
        );
    }

    let matrix = rmsd_matrix(structures, rmsd_params);
    let ref_idx = pick_reference_index(records, paths, reference_mode)?;
    let deltas = feature_deltas(records, ref_idx, delta_fields);
    let reference_path = paths[ref_idx].to_string_lossy().into_owned();

    Ok(CompareResult {
        reference_id: records[ref_idx].structure_id.clone(),
        reference_path,
        reference_mode: reference_mode_label(reference_mode),
        structure_count: structures.len(),
        structures: matrix.labels.clone(),
        matrix,
        deltas,
        failed,
    })
}

pub fn write_compare_json(out_dir: &Path, result: &CompareResult) -> Result<()> {
    fs::create_dir_all(out_dir)
        .with_context(|| format!("failed to create output directory {}", out_dir.display()))?;

    let matrix_path = out_dir.join("matrix.json");
    let matrix_json = serde_json::to_string_pretty(&result.matrix)?;
    fs::write(&matrix_path, matrix_json)
        .with_context(|| format!("failed to write {}", matrix_path.display()))?;

    let deltas_path = out_dir.join("deltas.jsonl");
    let mut writer = BufWriter::new(
        File::create(&deltas_path)
            .with_context(|| format!("failed to create {}", deltas_path.display()))?,
    );
    for delta in &result.deltas {
        writeln!(writer, "{}", serde_json::to_string(delta)?)?;
    }
    writer.flush()?;

    Ok(())
}

pub fn write_compare_csv(out_dir: &Path, result: &CompareResult) -> Result<()> {
    fs::create_dir_all(out_dir)
        .with_context(|| format!("failed to create output directory {}", out_dir.display()))?;

    write_matrix_csv(&out_dir.join("matrix.csv"), &result.matrix)?;
    write_deltas_csv(&out_dir.join("deltas.csv"), &result.deltas)?;

    Ok(())
}

fn write_matrix_csv(path: &Path, matrix: &RmsdMatrix) -> Result<()> {
    let mut lines = vec![format!(",{}", matrix.labels.join(","))];
    for (row_index, label) in matrix.labels.iter().enumerate() {
        let cells: Vec<String> = matrix.rmsd[row_index]
            .iter()
            .map(|value| value.map(|v| v.to_string()).unwrap_or_default())
            .collect();
        lines.push(format!("{label},{}", cells.join(",")));
    }
    fs::write(path, lines.join("\n")).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn write_deltas_csv(path: &Path, deltas: &[DeltaRecord]) -> Result<()> {
    let mut field_names: Vec<String> = deltas
        .iter()
        .flat_map(|delta| delta.fields.keys().cloned())
        .collect();
    field_names.sort();
    field_names.dedup();

    let mut lines = vec![format!(
        "structure_id,source_path,{}",
        field_names.join(",")
    )];
    for delta in deltas {
        let source_path = delta.source_path.as_deref().unwrap_or("");
        let values: Vec<String> = field_names
            .iter()
            .map(|field| csv_cell(delta.fields.get(field)))
            .collect();
        lines.push(format!(
            "{},{},{}",
            csv_cell_str(&delta.structure_id),
            csv_cell_str(source_path),
            values.join(",")
        ));
    }
    fs::write(path, lines.join("\n")).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn csv_cell_str(value: &str) -> String {
    if value.contains(',') || value.contains('"') || value.contains('\n') {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

fn csv_cell(value: Option<&Value>) -> String {
    match value {
        Some(Value::String(text)) => csv_cell_str(text),
        Some(other) => other.to_string(),
        None => String::new(),
    }
}

pub fn rmsd_matrix(structures: &[Structure], params: &RmsdParams) -> RmsdMatrix {
    let n = structures.len();
    let labels: Vec<String> = structures.iter().map(|s| s.id.clone()).collect();
    let mut rmsd = vec![vec![None; n]; n];

    for i in 0..n {
        rmsd[i][i] = Some(0.0);
    }

    for i in 0..n {
        for j in (i + 1)..n {
            let value = superposition_rmsd(&structures[i], &structures[j], params.clone())
                .ok()
                .map(|result| result.rmsd);
            rmsd[i][j] = value;
            rmsd[j][i] = value;
        }
    }

    RmsdMatrix { labels, rmsd }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::{json, Map};
    use structscope_core::{Atom, Chain, Residue, Structure, StructureMetadata};

    fn feature_record(id: &str, features: Map<String, Value>) -> FeatureRecord {
        FeatureRecord {
            structure_id: id.to_string(),
            source_path: None,
            features,
        }
    }

    fn quality_record(id: &str, outliers: u64, clashes: u64, missing_backbone: u64) -> FeatureRecord {
        let mut features = Map::new();
        features.insert("ramachandran_outlier_count".to_string(), json!(outliers));
        features.insert("clash_pair_count".to_string(), json!(clashes));
        features.insert(
            "missing_backbone_residue_count".to_string(),
            json!(missing_backbone),
        );
        feature_record(id, features)
    }

    fn sasa_record(id: &str, sasa_total: f64) -> FeatureRecord {
        let mut features = Map::new();
        features.insert("sasa_total".to_string(), json!(sasa_total));
        feature_record(id, features)
    }

    #[test]
    fn auto_quality_picks_lowest_outlier_count() {
        let records = vec![
            quality_record("high_outliers", 5, 0, 0),
            quality_record("best", 1, 10, 5),
            quality_record("medium_outliers", 3, 0, 0),
        ];
        let paths = vec![
            PathBuf::from("a.pdb"),
            PathBuf::from("b.pdb"),
            PathBuf::from("c.pdb"),
        ];

        let index = pick_reference_index(&records, &paths, &ReferenceMode::AutoQuality).unwrap();

        assert_eq!(index, 1);
    }

    #[test]
    fn numeric_feature_keys_returns_sorted_numeric_fields() {
        let mut features = Map::new();
        features.insert("sasa_total".to_string(), json!(100.0));
        features.insert("centroid".to_string(), json!([1.0, 2.0, 3.0]));
        features.insert("residue_count".to_string(), json!(42));
        let record = feature_record("test", features);

        assert_eq!(
            numeric_feature_keys(&record),
            vec!["residue_count".to_string(), "sasa_total".to_string()]
        );
    }

    #[test]
    fn feature_deltas_compute_sasa_total_difference() {
        let records = vec![
            sasa_record("ref", 1180.0),
            sasa_record("mob", 1200.0),
        ];

        let deltas = feature_deltas(&records, 0, Some(&["sasa_total".to_string()]));

        assert_eq!(deltas.len(), 1);
        assert_eq!(deltas[0].structure_id, "mob");
        assert_eq!(deltas[0].fields["sasa_total"].as_f64(), Some(1200.0));
        assert_eq!(deltas[0].fields["reference_sasa_total"].as_f64(), Some(1180.0));
        assert_eq!(deltas[0].fields["delta_sasa_total"].as_f64(), Some(20.0));
    }

    #[test]
    fn min_sasa_total_picks_correct_index() {
        let records = vec![
            sasa_record("large", 1500.0),
            sasa_record("smallest", 900.0),
            sasa_record("middle", 1200.0),
        ];
        let paths = vec![
            PathBuf::from("large.pdb"),
            PathBuf::from("smallest.pdb"),
            PathBuf::from("middle.pdb"),
        ];
        let mode = ReferenceMode::ByField {
            maximize: false,
            field: "sasa_total".to_string(),
        };

        let index = pick_reference_index(&records, &paths, &mode).unwrap();

        assert_eq!(index, 1);
    }

    fn ca_atom(x: f64, y: f64, z: f64) -> Atom {
        Atom {
            id: "CA".to_string(),
            serial: Some(1),
            name: "CA".to_string(),
            element: Some("C".to_string()),
            x,
            y,
            z,
            occupancy: None,
            temp_factor: None,
        }
    }

    fn residue(name: &str, ca: [f64; 3]) -> Residue {
        Residue {
            id: name.to_string(),
            name: name.to_string(),
            seq_number: 1,
            insertion_code: None,
            atoms: vec![ca_atom(ca[0], ca[1], ca[2])],
            is_hetero: false,
        }
    }

    fn structure(id: &str, residues: &[(&str, [f64; 3])]) -> Structure {
        Structure {
            id: id.to_string(),
            metadata: StructureMetadata {
                source_format: "test".to_string(),
                source_path: None,
                title: None,
            },
            chains: vec![Chain {
                id: "A".to_string(),
                label: "A".to_string(),
                residues: residues
                    .iter()
                    .map(|(name, ca)| residue(name, *ca))
                    .collect(),
            }],
        }
    }

    fn default_params() -> RmsdParams {
        RmsdParams {
            atoms: "ca".to_string(),
            align: false,
            local: false,
        }
    }

    #[test]
    fn identical_structures_have_zero_off_diagonal_rmsd() {
        let coords = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let residues = [("ALA", coords[0]), ("GLY", coords[1]), ("VAL", coords[2])];
        let structures = [
            structure("s1", &residues),
            structure("s2", &residues),
            structure("s3", &residues),
        ];

        let matrix = rmsd_matrix(&structures, &default_params());

        assert_eq!(matrix.labels, vec!["s1", "s2", "s3"]);
        assert_eq!(matrix.rmsd.len(), 3);
        for i in 0..3 {
            assert_eq!(matrix.rmsd[i][i], Some(0.0));
            for j in 0..3 {
                if i != j {
                    let value = matrix.rmsd[i][j].expect("expected rmsd value");
                    assert!(value < 1e-9, "rmsd[{i}][{j}] = {value}");
                    assert_eq!(matrix.rmsd[j][i], matrix.rmsd[i][j]);
                }
            }
        }
    }

    #[test]
    fn offset_structures_have_nonzero_rmsd() {
        let reference = structure(
            "ref",
            &[
                ("ALA", [0.0, 0.0, 0.0]),
                ("GLY", [1.0, 0.0, 0.0]),
                ("VAL", [0.0, 1.0, 0.0]),
            ],
        );
        let mobile = structure(
            "mob",
            &[
                ("ALA", [0.0, 0.0, 0.0]),
                ("GLY", [2.0, 0.0, 0.0]),
                ("VAL", [0.0, 1.0, 0.0]),
            ],
        );

        let matrix = rmsd_matrix(&[reference, mobile], &default_params());

        assert_eq!(matrix.labels, vec!["ref", "mob"]);
        let off_diagonal = matrix.rmsd[0][1].expect("expected rmsd value");
        assert!(off_diagonal > 0.0, "rmsd = {off_diagonal}");
        assert_eq!(matrix.rmsd[1][0], Some(off_diagonal));
    }

    const MINI_PDB: &str = "\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  GLY A   2       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  CA  VAL A   3       0.000   1.000   0.000  1.00  0.00           C
";

    fn mini_pdb(offset_x: f64) -> String {
        MINI_PDB.replace("0.000   0.000   0.000", &format!("{offset_x:.3}   0.000   0.000"))
    }

    #[test]
    fn compare_set_assembles_matrix_and_deltas() {
        use structscope_core::{parse_str, InputFormat, ParseOptions};

        let structures = [
            parse_str(
                &mini_pdb(0.0),
                InputFormat::Pdb,
                Some("s1.pdb".to_string()),
                ParseOptions::default(),
            )
            .unwrap(),
            parse_str(
                &mini_pdb(0.5),
                InputFormat::Pdb,
                Some("s2.pdb".to_string()),
                ParseOptions::default(),
            )
            .unwrap(),
            parse_str(
                &mini_pdb(1.0),
                InputFormat::Pdb,
                Some("s3.pdb".to_string()),
                ParseOptions::default(),
            )
            .unwrap(),
        ];
        let records = vec![
            sasa_record("s1", 1000.0),
            sasa_record("s2", 1100.0),
            sasa_record("s3", 1200.0),
        ];
        let paths = vec![
            PathBuf::from("s1.pdb"),
            PathBuf::from("s2.pdb"),
            PathBuf::from("s3.pdb"),
        ];

        let result = compare_set(
            &structures,
            &records,
            &paths,
            &default_params(),
            &ReferenceMode::FirstInput,
            Some(&["sasa_total".to_string()]),
            None,
        )
        .unwrap();

        assert_eq!(result.reference_id, "s1");
        assert_eq!(result.reference_path, "s1.pdb");
        assert_eq!(result.reference_mode, "first_input");
        assert_eq!(result.structure_count, 3);
        assert_eq!(result.structures, vec!["s1", "s2", "s3"]);
        assert_eq!(result.matrix.labels.len(), 3);
        assert_eq!(result.deltas.len(), 2);
        assert_eq!(result.deltas[0].structure_id, "s2");
        assert_eq!(result.deltas[0].fields["delta_sasa_total"].as_f64(), Some(100.0));
        assert_eq!(result.deltas[1].fields["delta_sasa_total"].as_f64(), Some(200.0));
        assert!(result.failed.is_none());

        for i in 0..3 {
            assert_eq!(result.matrix.rmsd[i][i], Some(0.0));
        }
    }

    #[test]
    fn compare_set_writes_json_and_csv_outputs() {
        use structscope_core::{parse_str, InputFormat, ParseOptions};

        let structures = [
            parse_str(
                &mini_pdb(0.0),
                InputFormat::Pdb,
                Some("s1.pdb".to_string()),
                ParseOptions::default(),
            )
            .unwrap(),
            parse_str(
                &mini_pdb(1.0),
                InputFormat::Pdb,
                Some("s2.pdb".to_string()),
                ParseOptions::default(),
            )
            .unwrap(),
        ];
        let records = vec![sasa_record("s1", 900.0), sasa_record("s2", 950.0)];
        let paths = vec![PathBuf::from("s1.pdb"), PathBuf::from("s2.pdb")];
        let result = compare_set(
            &structures,
            &records,
            &paths,
            &default_params(),
            &ReferenceMode::FirstInput,
            Some(&["sasa_total".to_string()]),
            None,
        )
        .unwrap();

        let out_dir = std::env::temp_dir().join(format!(
            "structscope-compare-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&out_dir);
        write_compare_json(&out_dir, &result).unwrap();
        write_compare_csv(&out_dir, &result).unwrap();

        assert!(out_dir.join("matrix.json").is_file());
        assert!(out_dir.join("deltas.jsonl").is_file());
        assert!(out_dir.join("matrix.csv").is_file());
        assert!(out_dir.join("deltas.csv").is_file());

        let matrix_json = fs::read_to_string(out_dir.join("matrix.json")).unwrap();
        assert!(matrix_json.contains("\"labels\""));
        let deltas_jsonl = fs::read_to_string(out_dir.join("deltas.jsonl")).unwrap();
        assert!(deltas_jsonl.contains("delta_sasa_total"));

        let _ = fs::remove_dir_all(&out_dir);
    }
}
