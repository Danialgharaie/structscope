use std::fs;
use std::path::PathBuf;

use serde::{Deserialize, Serialize};
use structscope_core::{parse_str, InputFormat, ParseOptions};
use structscope_features::quality::{per_residue_quality, quality_summary, QualityParams, QualityRecord};

/// Three-residue clash panel (res 1 vs 3) plus ALA 4 missing O.
const PANEL: &str = "\
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
ATOM     13  N   ALA A   4       5.000   0.000   0.000  1.00 0.00           N
ATOM     14  CA  ALA A   4       6.458   0.000   0.000  1.00 0.00           C
ATOM     15  C   ALA A   4       7.009   1.420   0.000  1.00 0.00           C
";

/// Clean extended dipeptide; sequential neighbors should not clash.
const DIPEPTIDE: &str = "\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 0.00           C
ATOM      4  O   ALA A   1       2.500   2.200   0.000  1.00 0.00           O
ATOM      5  N   ALA A   2       3.332   1.540   0.000  1.00 0.00           N
ATOM      6  CA  ALA A   2       3.970   2.840   0.000  1.00 0.00           C
ATOM      7  C   ALA A   2       5.480   2.700   0.000  1.00 0.00           C
ATOM      8  O   ALA A   2       6.000   3.600   0.000  1.00 0.00           O
";

const COUNT_TOLERANCE: usize = 0;
const PHI_PSI_ABS_TOLERANCE: f64 = 1e-6;

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/quality")
}

#[derive(Debug, Serialize, Deserialize)]
struct QualityGoldenParams {
    clash_overlap: f64,
}

#[derive(Debug, Serialize, Deserialize, PartialEq)]
struct QualityGoldenSummary {
    quality_residue_count: usize,
    ramachandran_evaluated_count: usize,
    ramachandran_favored_count: usize,
    ramachandran_allowed_count: usize,
    ramachandran_outlier_count: usize,
    clash_pair_count: usize,
    missing_backbone_residue_count: usize,
}

#[derive(Debug, Serialize, Deserialize, PartialEq)]
struct QualityGoldenRecord {
    chain_id: String,
    residue_name: String,
    seq_number: i32,
    phi: Option<f64>,
    psi: Option<f64>,
    ramachandran: Option<String>,
    missing_backbone: Vec<String>,
    clash_count: usize,
}

#[derive(Debug, Serialize, Deserialize)]
struct QualityGoldenFixture {
    pdb: String,
    params: QualityGoldenParams,
    expected: QualityGoldenExpected,
}

#[derive(Debug, Serialize, Deserialize)]
struct QualityGoldenExpected {
    summary: QualityGoldenSummary,
    problems_only_count: usize,
    problems: Vec<QualityGoldenRecord>,
    all_residues_count: usize,
}

fn params_from_fixture(fixture: &QualityGoldenFixture) -> QualityParams {
    QualityParams {
        clash_overlap: fixture.params.clash_overlap,
    }
}

fn summary_from(actual: &structscope_features::quality::QualitySummary) -> QualityGoldenSummary {
    QualityGoldenSummary {
        quality_residue_count: actual.quality_residue_count,
        ramachandran_evaluated_count: actual.ramachandran_evaluated_count,
        ramachandran_favored_count: actual.ramachandran_favored_count,
        ramachandran_allowed_count: actual.ramachandran_allowed_count,
        ramachandran_outlier_count: actual.ramachandran_outlier_count,
        clash_pair_count: actual.clash_pair_count,
        missing_backbone_residue_count: actual.missing_backbone_residue_count,
    }
}

fn record_from(actual: &QualityRecord) -> QualityGoldenRecord {
    QualityGoldenRecord {
        chain_id: actual.chain_id.clone(),
        residue_name: actual.residue_name.clone(),
        seq_number: actual.seq_number,
        phi: actual.phi,
        psi: actual.psi,
        ramachandran: actual.ramachandran.clone(),
        missing_backbone: actual.missing_backbone.clone(),
        clash_count: actual.clash_count,
    }
}

fn assert_usize_eq(actual: usize, expected: usize, label: &str) {
    assert_eq!(
        actual, expected,
        "{label}: expected {expected}, got {actual} (tol {COUNT_TOLERANCE})"
    );
}

fn assert_optional_f64(actual: Option<f64>, expected: Option<f64>, label: &str) {
    match (actual, expected) {
        (Some(a), Some(e)) => {
            let diff = (a - e).abs();
            assert!(
                diff <= PHI_PSI_ABS_TOLERANCE,
                "{label}: expected {e}, got {a} (diff {diff})"
            );
        }
        (None, None) => {}
        (a, e) => panic!("{label}: expected {e:?}, got {a:?}"),
    }
}

fn assert_record(actual: &QualityRecord, expected: &QualityGoldenRecord) {
    assert_eq!(actual.chain_id, expected.chain_id);
    assert_eq!(actual.residue_name, expected.residue_name);
    assert_eq!(actual.seq_number, expected.seq_number);
    assert_optional_f64(actual.phi, expected.phi, "phi");
    assert_optional_f64(actual.psi, expected.psi, "psi");
    assert_eq!(actual.ramachandran, expected.ramachandran);
    assert_eq!(actual.missing_backbone, expected.missing_backbone);
    assert_eq!(actual.clash_count, expected.clash_count);
}

fn compute_and_assert(fixture: &QualityGoldenFixture) {
    let structure = parse_str(&fixture.pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
    let params = params_from_fixture(fixture);

    let summary = quality_summary(&structure, &params);
    let actual_summary = summary_from(&summary);
    assert_eq!(actual_summary, fixture.expected.summary);

    let problems = per_residue_quality(&structure, &params, false);
    assert_usize_eq(problems.len(), fixture.expected.problems_only_count, "problems_only_count");
    assert_eq!(problems.len(), fixture.expected.problems.len());
    for (actual, expected) in problems.iter().zip(&fixture.expected.problems) {
        assert_record(actual, expected);
    }

    let all = per_residue_quality(&structure, &params, true);
    assert_usize_eq(all.len(), fixture.expected.all_residues_count, "all_residues_count");
}

fn write_fixture(path: &str, pdb: &str, clash_overlap: f64) {
    let structure = parse_str(pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
    let params = QualityParams { clash_overlap };
    let summary = quality_summary(&structure, &params);
    let problems = per_residue_quality(&structure, &params, false);
    let all = per_residue_quality(&structure, &params, true);

    let fixture = QualityGoldenFixture {
        pdb: pdb.to_string(),
        params: QualityGoldenParams { clash_overlap },
        expected: QualityGoldenExpected {
            summary: summary_from(&summary),
            problems_only_count: problems.len(),
            problems: problems.iter().map(record_from).collect(),
            all_residues_count: all.len(),
        },
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
fn write_synthetic_panel_golden_fixture() {
    write_fixture("synthetic_panel.json", PANEL, 0.4);
}

#[test]
#[ignore]
fn write_synthetic_dipeptide_golden_fixture() {
    write_fixture("synthetic_dipeptide.json", DIPEPTIDE, 0.4);
}

#[test]
fn synthetic_panel_matches_golden() {
    let fixture: QualityGoldenFixture =
        serde_json::from_str(include_str!("../../../tests/fixtures/quality/synthetic_panel.json")).unwrap();
    compute_and_assert(&fixture);
}

#[test]
fn synthetic_dipeptide_matches_golden() {
    let fixture: QualityGoldenFixture =
        serde_json::from_str(include_str!("../../../tests/fixtures/quality/synthetic_dipeptide.json")).unwrap();
    compute_and_assert(&fixture);
}

#[test]
#[ignore = "requires offline MolProbity reference values in molprobity_reference.json"]
fn molprobity_validated_structure() {
    let fixture_path = fixtures_dir().join("molprobity_reference.json");
    assert!(
        fixture_path.exists(),
        "missing {}; run offline MolProbity validation and add fixture",
        fixture_path.display()
    );
    let fixture_json = fs::read_to_string(&fixture_path).unwrap();
    let fixture: QualityGoldenFixture = serde_json::from_str(&fixture_json).unwrap();
    compute_and_assert(&fixture);
}
