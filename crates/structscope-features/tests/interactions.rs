use std::path::PathBuf;
use structscope_core::{parse_file, ParseOptions};
use structscope_features::interactions::interactions;

fn corpus_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/corpus/rcsb")
}

/// 1bju (beta-trypsin) has six disulfide bridges and many polar contacts.
#[test]
fn detects_disulfides_and_contacts_on_real_structure() {
    let s = parse_file(&corpus_root().join("mmcif/1bju.cif.gz"), ParseOptions::default()).unwrap();
    let found = interactions(&s);
    let count = |k: &str| found.iter().filter(|i| i.kind == k).count();
    assert_eq!(count("disulfide"), 6, "beta-trypsin has six disulfides");
    assert!(count("hydrogen_bond") > 100, "expected many polar contacts");
    // Every reported distance must be within the loosest cutoff used (6.0A for cation-pi).
    assert!(found.iter().all(|i| i.distance <= 6.0));
    assert!(count("cation_pi") > 0, "expected some cation-pi interactions");
    assert!(count("hydrophobic") > 0, "expected many hydrophobic contacts");
}
