use std::fs;
use std::path::PathBuf;

use serde::{Deserialize, Serialize};
use structscope_core::{parse_str, InputFormat, ParseOptions};
use structscope_features::interface::{per_interface_features, InterfaceParams};

const DIMER: &str = "\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 0.00           C
ATOM      2  CA  ALA B   1       3.500   0.000   0.000  1.00 0.00           C
";

const SYNTHETIC_TOLERANCE: f64 = 1e-6;
const CCP4_BSA_REL_TOLERANCE: f64 = 0.02;
const CCP4_AREA_REL_TOLERANCE: f64 = 0.02;
const CCP4_SC_ABS_TOLERANCE: f64 = 0.05;

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/interface")
}

#[derive(Debug, Serialize, Deserialize)]
struct InterfaceGoldenParams {
    contact_distance: f64,
    area_distance: f64,
    sc_distance: f64,
}

#[derive(Debug, Serialize, Deserialize)]
struct InterfaceGoldenExpected {
    chain_a: String,
    chain_b: String,
    bsa: f64,
    interface_area: f64,
    shape_complementarity: f64,
}

#[derive(Debug, Serialize, Deserialize)]
struct InterfaceGoldenFixture {
    pdb: String,
    params: InterfaceGoldenParams,
    expected: InterfaceGoldenExpected,
}

fn default_params() -> InterfaceParams {
    InterfaceParams {
        contact_distance: 8.0,
        area_distance: 5.0,
        sc_distance: 5.0,
    }
}

fn params_from_fixture(fixture: &InterfaceGoldenFixture) -> InterfaceParams {
    InterfaceParams {
        contact_distance: fixture.params.contact_distance,
        area_distance: fixture.params.area_distance,
        sc_distance: fixture.params.sc_distance,
    }
}

fn assert_abs_within(actual: f64, expected: f64, tolerance: f64, label: &str) {
    let diff = (actual - expected).abs();
    assert!(
        diff <= tolerance,
        "{label}: expected {expected}, got {actual} (diff {diff}, tol {tolerance})"
    );
}

fn assert_rel_within(actual: f64, expected: f64, rel_tolerance: f64, label: &str) {
    let denom = expected.abs().max(1e-12);
    let rel_diff = (actual - expected).abs() / denom;
    assert!(
        rel_diff <= rel_tolerance,
        "{label}: expected {expected}, got {actual} (rel diff {rel_diff}, tol {rel_tolerance})"
    );
}

fn compute_and_assert(fixture: &InterfaceGoldenFixture, bsa_tol: f64, area_tol: f64, sc_tol: f64, use_rel: bool) {
    let structure = parse_str(&fixture.pdb, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
    let params = params_from_fixture(fixture);
    let records = per_interface_features(&structure, &params);
    assert_eq!(records.len(), 1);
    let record = &records[0];
    assert_eq!(record.chain_a, fixture.expected.chain_a);
    assert_eq!(record.chain_b, fixture.expected.chain_b);

    if use_rel {
        assert_rel_within(record.bsa, fixture.expected.bsa, bsa_tol, "bsa");
        assert_rel_within(record.interface_area, fixture.expected.interface_area, area_tol, "interface_area");
    } else {
        assert_abs_within(record.bsa, fixture.expected.bsa, bsa_tol, "bsa");
        assert_abs_within(record.interface_area, fixture.expected.interface_area, area_tol, "interface_area");
    }
    assert_abs_within(
        record.shape_complementarity,
        fixture.expected.shape_complementarity,
        sc_tol,
        "shape_complementarity",
    );
}

#[test]
#[ignore]
fn write_synthetic_golden_fixture() {
    let structure = parse_str(DIMER, InputFormat::Pdb, None, ParseOptions::default()).unwrap();
    let params = default_params();
    let records = per_interface_features(&structure, &params);
    assert_eq!(records.len(), 1);
    let record = &records[0];

    let fixture = InterfaceGoldenFixture {
        pdb: DIMER.to_string(),
        params: InterfaceGoldenParams {
            contact_distance: params.contact_distance,
            area_distance: params.area_distance,
            sc_distance: params.sc_distance,
        },
        expected: InterfaceGoldenExpected {
            chain_a: record.chain_a.clone(),
            chain_b: record.chain_b.clone(),
            bsa: record.bsa,
            interface_area: record.interface_area,
            shape_complementarity: record.shape_complementarity,
        },
    };

    let path = fixtures_dir().join("synthetic_dimer.json");
    fs::create_dir_all(path.parent().unwrap()).unwrap();
    let json = serde_json::to_string_pretty(&fixture).unwrap();
    fs::write(&path, format!("{json}\n")).unwrap();
    println!("wrote {}", path.display());
    println!("{json}");
}

#[test]
fn synthetic_dimer_matches_golden() {
    let fixture: InterfaceGoldenFixture =
        serde_json::from_str(include_str!("../../../tests/fixtures/interface/synthetic_dimer.json")).unwrap();
    compute_and_assert(&fixture, SYNTHETIC_TOLERANCE, SYNTHETIC_TOLERANCE, SYNTHETIC_TOLERANCE, false);
}

#[test]
#[ignore = "requires offline CCP4/PISA values in ccp4_dimer.json"]
fn ccp4_validated_dimer() {
    let fixture_path = fixtures_dir().join("ccp4_dimer.json");
    assert!(
        fixture_path.exists(),
        "missing {}; run offline CCP4/PISA validation and add fixture",
        fixture_path.display()
    );
    let fixture_json = fs::read_to_string(&fixture_path).unwrap();
    let fixture: InterfaceGoldenFixture = serde_json::from_str(&fixture_json).unwrap();
    compute_and_assert(
        &fixture,
        CCP4_BSA_REL_TOLERANCE,
        CCP4_AREA_REL_TOLERANCE,
        CCP4_SC_ABS_TOLERANCE,
        true,
    );
}
