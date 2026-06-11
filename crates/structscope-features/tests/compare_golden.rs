use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use serde_json::{json, Map};
use structscope_core::{parse_str, InputFormat, ParseOptions, RmsdParams};
use structscope_features::compare::{
    compare_set, feature_deltas, rmsd_matrix, CompareResult, DeltaRecord, ReferenceMode,
    RmsdMatrix,
};
use structscope_features::FeatureRecord;

const MINI_PDB: &str = "\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  GLY A   2       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  CA  VAL A   3       0.000   1.000   0.000  1.00  0.00           C
";

const PATHS: [&str; 3] = ["s1.pdb", "s2.pdb", "s3.pdb"];
const OFFSETS: [f64; 3] = [0.0, 0.5, 1.0];
const SASA_TOTALS: [f64; 3] = [1000.0, 1100.0, 1200.0];
const RMSD_ABS_TOLERANCE: f64 = 1e-6;

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/compare")
}

fn mini_pdb(offset_x: f64) -> String {
    MINI_PDB.replace("0.000   0.000   0.000", &format!("{offset_x:.3}   0.000   0.000"))
}

fn default_rmsd_params() -> RmsdParams {
    RmsdParams {
        atoms: "ca".to_string(),
        align: false,
        local: false,
    }
}

#[derive(Debug, Serialize, Deserialize)]
struct CompareGoldenRmsdParams {
    atoms: String,
    align: bool,
    local: bool,
}

#[derive(Debug, Serialize, Deserialize)]
struct CompareGoldenParams {
    rmsd: CompareGoldenRmsdParams,
    reference_mode: String,
    delta_fields: Vec<String>,
}

#[derive(Debug, Serialize, Deserialize)]
struct CompareGoldenFeature {
    structure_id: String,
    sasa_total: f64,
}

#[derive(Debug, Serialize, Deserialize, PartialEq)]
struct CompareGoldenMatrix {
    labels: Vec<String>,
    rmsd: Vec<Vec<Option<f64>>>,
}

#[derive(Debug, Serialize, Deserialize, PartialEq)]
struct CompareGoldenDelta {
    structure_id: String,
    source_path: Option<String>,
    sasa_total: f64,
    reference_sasa_total: f64,
    delta_sasa_total: f64,
}

#[derive(Debug, Serialize, Deserialize, PartialEq)]
struct CompareGoldenExpected {
    reference_id: String,
    reference_path: String,
    reference_mode: String,
    structure_count: usize,
    structures: Vec<String>,
    matrix: CompareGoldenMatrix,
    deltas_count: usize,
    deltas: Vec<CompareGoldenDelta>,
}

#[derive(Debug, Serialize, Deserialize)]
struct CompareGoldenFixture {
    pdbs: Vec<String>,
    paths: Vec<String>,
    params: CompareGoldenParams,
    features: Vec<CompareGoldenFeature>,
    expected: CompareGoldenExpected,
}

fn rmsd_params_from_fixture(params: &CompareGoldenRmsdParams) -> RmsdParams {
    RmsdParams {
        atoms: params.atoms.clone(),
        align: params.align,
        local: params.local,
    }
}

fn reference_mode_from_fixture(mode: &str) -> ReferenceMode {
    match mode {
        "first_input" => ReferenceMode::FirstInput,
        other => panic!("unsupported reference_mode in fixture: {other}"),
    }
}

fn feature_records_from_fixture(features: &[CompareGoldenFeature]) -> Vec<FeatureRecord> {
    features
        .iter()
        .map(|feature| {
            let mut map = Map::new();
            map.insert("sasa_total".to_string(), json!(feature.sasa_total));
            FeatureRecord {
                structure_id: feature.structure_id.clone(),
                source_path: None,
                features: map,
            }
        })
        .collect()
}

fn matrix_from(actual: &RmsdMatrix) -> CompareGoldenMatrix {
    CompareGoldenMatrix {
        labels: actual.labels.clone(),
        rmsd: actual.rmsd.clone(),
    }
}

fn delta_from(actual: &DeltaRecord) -> CompareGoldenDelta {
    CompareGoldenDelta {
        structure_id: actual.structure_id.clone(),
        source_path: actual.source_path.clone(),
        sasa_total: actual.fields["sasa_total"].as_f64().expect("sasa_total"),
        reference_sasa_total: actual.fields["reference_sasa_total"]
            .as_f64()
            .expect("reference_sasa_total"),
        delta_sasa_total: actual.fields["delta_sasa_total"]
            .as_f64()
            .expect("delta_sasa_total"),
    }
}

fn expected_from(result: &CompareResult) -> CompareGoldenExpected {
    CompareGoldenExpected {
        reference_id: result.reference_id.clone(),
        reference_path: result.reference_path.clone(),
        reference_mode: result.reference_mode.clone(),
        structure_count: result.structure_count,
        structures: result.structures.clone(),
        matrix: matrix_from(&result.matrix),
        deltas_count: result.deltas.len(),
        deltas: result.deltas.iter().map(delta_from).collect(),
    }
}

fn assert_optional_f64(actual: Option<f64>, expected: Option<f64>, label: &str) {
    match (actual, expected) {
        (Some(a), Some(e)) => {
            let diff = (a - e).abs();
            assert!(
                diff <= RMSD_ABS_TOLERANCE,
                "{label}: expected {e}, got {a} (diff {diff})"
            );
        }
        (None, None) => {}
        (a, e) => panic!("{label}: expected {e:?}, got {a:?}"),
    }
}

fn assert_matrix(actual: &RmsdMatrix, expected: &CompareGoldenMatrix) {
    assert_eq!(actual.labels, expected.labels);
    assert_eq!(actual.rmsd.len(), expected.rmsd.len());
    for (row_index, (actual_row, expected_row)) in actual.rmsd.iter().zip(&expected.rmsd).enumerate() {
        assert_eq!(actual_row.len(), expected_row.len(), "row {row_index} width");
        for (col_index, (actual_value, expected_value)) in actual_row.iter().zip(expected_row).enumerate() {
            assert_optional_f64(*actual_value, *expected_value, &format!("rmsd[{row_index}][{col_index}]"));
        }
    }
}

fn assert_delta(actual: &DeltaRecord, expected: &CompareGoldenDelta) {
    assert_eq!(actual.structure_id, expected.structure_id);
    assert_eq!(actual.source_path, expected.source_path);
    let actual_delta = delta_from(actual);
    assert_eq!(actual_delta, *expected);
}

fn compute_and_assert(fixture: &CompareGoldenFixture) {
    let structures: Vec<_> = fixture
        .pdbs
        .iter()
        .zip(&fixture.paths)
        .map(|(pdb, path)| {
            parse_str(
                pdb,
                InputFormat::Pdb,
                Some(path.clone()),
                ParseOptions::default(),
            )
            .unwrap()
        })
        .collect();
    let records = feature_records_from_fixture(&fixture.features);
    let paths: Vec<PathBuf> = fixture.paths.iter().map(PathBuf::from).collect();
    let rmsd_params = rmsd_params_from_fixture(&fixture.params.rmsd);
    let reference_mode = reference_mode_from_fixture(&fixture.params.reference_mode);
    let delta_fields: Vec<String> = fixture.params.delta_fields.clone();

    let result = compare_set(
        &structures,
        &records,
        &paths,
        &rmsd_params,
        &reference_mode,
        Some(&delta_fields),
        None,
    )
    .unwrap();

    assert_eq!(result.reference_id, fixture.expected.reference_id);
    assert_eq!(result.reference_path, fixture.expected.reference_path);
    assert_eq!(result.reference_mode, fixture.expected.reference_mode);
    assert_eq!(result.structure_count, fixture.expected.structure_count);
    assert_eq!(result.structures, fixture.expected.structures);
    assert_matrix(&result.matrix, &fixture.expected.matrix);
    assert_eq!(result.deltas.len(), fixture.expected.deltas_count);
    assert_eq!(result.deltas.len(), fixture.expected.deltas.len());
    for (actual, expected) in result.deltas.iter().zip(&fixture.expected.deltas) {
        assert_delta(actual, expected);
    }

    let matrix_only = rmsd_matrix(&structures, &rmsd_params);
    assert_matrix(&matrix_only, &fixture.expected.matrix);

    let ref_idx = 0;
    let deltas_only = feature_deltas(&records, ref_idx, Some(&delta_fields));
    assert_eq!(deltas_only.len(), fixture.expected.deltas_count);
    for (actual, expected) in deltas_only.iter().zip(&fixture.expected.deltas) {
        assert_delta(actual, expected);
    }
}

fn triplet_pdbs() -> Vec<String> {
    OFFSETS.iter().map(|offset| mini_pdb(*offset)).collect()
}

fn triplet_features() -> Vec<CompareGoldenFeature> {
    PATHS
        .iter()
        .zip(SASA_TOTALS)
        .map(|(path, sasa_total)| CompareGoldenFeature {
            structure_id: Path::new(path)
                .file_stem()
                .and_then(|stem| stem.to_str())
                .unwrap_or("structure")
                .to_string(),
            sasa_total,
        })
        .collect()
}

fn write_fixture(path: &str) {
    let pdbs = triplet_pdbs();
    let paths: Vec<String> = PATHS.iter().map(|p| (*p).to_string()).collect();
    let structures: Vec<_> = pdbs
        .iter()
        .zip(&paths)
        .map(|(pdb, path)| {
            parse_str(
                pdb,
                InputFormat::Pdb,
                Some(path.clone()),
                ParseOptions::default(),
            )
            .unwrap()
        })
        .collect();
    let features = triplet_features();
    let records = feature_records_from_fixture(&features);
    let path_bufs: Vec<PathBuf> = paths.iter().map(PathBuf::from).collect();
    let rmsd_params = default_rmsd_params();
    let delta_fields = vec!["sasa_total".to_string()];

    let result = compare_set(
        &structures,
        &records,
        &path_bufs,
        &rmsd_params,
        &ReferenceMode::FirstInput,
        Some(&delta_fields),
        None,
    )
    .unwrap();

    let fixture = CompareGoldenFixture {
        pdbs,
        paths,
        params: CompareGoldenParams {
            rmsd: CompareGoldenRmsdParams {
                atoms: rmsd_params.atoms,
                align: rmsd_params.align,
                local: rmsd_params.local,
            },
            reference_mode: "first_input".to_string(),
            delta_fields,
        },
        features,
        expected: expected_from(&result),
    };

    let out_path = fixtures_dir().join(path);
    fs::create_dir_all(out_path.parent().unwrap()).unwrap();
    let json = serde_json::to_string_pretty(&fixture).unwrap();
    fs::write(&out_path, format!("{json}\n")).unwrap();
    println!("wrote {}", out_path.display());
    println!("{json}");
}

#[test]
#[ignore]
fn write_synthetic_triplet_golden_fixture() {
    write_fixture("synthetic_triplet.json");
}

#[test]
fn synthetic_triplet_matches_golden() {
    let fixture: CompareGoldenFixture =
        serde_json::from_str(include_str!("../../../tests/fixtures/compare/synthetic_triplet.json")).unwrap();
    compute_and_assert(&fixture);
}
