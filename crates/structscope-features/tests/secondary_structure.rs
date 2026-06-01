use std::path::PathBuf;
use structscope_core::{parse_file, ParseOptions};
use structscope_features::ss::secondary_structure;

fn corpus_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/corpus/rcsb")
}

/// 1nkd (NK-lysin) is an all-alpha saposin fold: helices present, no beta strands.
#[test]
fn all_alpha_protein_has_helix_no_strand() {
    let s = parse_file(&corpus_root().join("mmcif/1nkd.cif.gz"), ParseOptions::default()).unwrap();
    let ss: String = secondary_structure(&s).iter().map(|c| c.ss.clone()).collect();
    let helix = ss.chars().filter(|c| matches!(c, 'H' | 'G' | 'I')).count();
    let strand = ss.chars().filter(|&c| c == 'E').count();
    assert!(helix > 10, "expected substantial helix content, got {helix}");
    assert_eq!(strand, 0, "all-alpha fold should have no strands, got {strand}");
}

/// 2aco is a mixed alpha/beta structure: both helix and strand must be present.
#[test]
fn mixed_protein_has_helix_and_strand() {
    let s = parse_file(&corpus_root().join("mmcif/2aco.cif.gz"), ParseOptions::default()).unwrap();
    let ss: String = secondary_structure(&s).iter().map(|c| c.ss.clone()).collect();
    assert!(ss.chars().any(|c| matches!(c, 'H' | 'G' | 'I')), "expected helices");
    assert!(ss.contains('E'), "expected strands");
}
