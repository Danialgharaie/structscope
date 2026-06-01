use std::fs;
use std::path::{Path, PathBuf};
use structscope_core::{parse_file, ParseOptions};

const TARGET_PDB_FILES: usize = 100;
const TARGET_MMCIF_FILES: usize = 100;

#[test]
fn parses_checked_in_rcsb_corpus() {
    let corpus_root = corpus_root();
    let mmcif_files = collect_fixture_files(&corpus_root.join("mmcif"), "cif.gz");
    let pdb_files = collect_fixture_files(&corpus_root.join("pdb"), "pdb.gz");

    assert_eq!(
        mmcif_files.len(),
        TARGET_MMCIF_FILES,
        "unexpected mmCIF fixture count under {}",
        corpus_root.join("mmcif").display()
    );
    assert_eq!(
        pdb_files.len(),
        TARGET_PDB_FILES,
        "unexpected PDB fixture count under {}",
        corpus_root.join("pdb").display()
    );

    for path in &mmcif_files {
        assert_parses(path, "mmCIF");
    }
    for path in &pdb_files {
        assert_parses(path, "PDB");
    }
}

fn corpus_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/corpus/rcsb")
}

fn collect_fixture_files(dir: &Path, suffix: &str) -> Vec<PathBuf> {
    let mut files = fs::read_dir(dir)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", dir.display()))
        .filter_map(|entry| entry.ok().map(|entry| entry.path()))
        .filter(|path| {
            path.file_name()
                .and_then(|name| name.to_str())
                .map(|name| name.ends_with(suffix))
                .unwrap_or(false)
        })
        .collect::<Vec<_>>();
    files.sort();
    files
}

fn assert_parses(path: &Path, label: &str) {
    let metadata = fs::metadata(path)
        .unwrap_or_else(|err| panic!("failed to stat {}: {err}", path.display()));
    assert!(
        metadata.len() > 0,
        "{label} fixture is empty: {}",
        path.display()
    );

    let structure = parse_file(path, ParseOptions::default())
        .unwrap_or_else(|err| panic!("failed to parse {label} fixture {}: {err}", path.display()));
    let summary = structure.summary();
    assert!(summary.chain_count > 0, "{label} fixture has no chains: {}", path.display());
    assert!(
        summary.residue_count > 0,
        "{label} fixture has no residues: {}",
        path.display()
    );
    assert!(summary.atom_count > 0, "{label} fixture has no atoms: {}", path.display());
}

#[test]
fn bcif_matches_mmcif_summary() {
    let root = corpus_root();
    for id in ["100d", "1nkd"] {
        let bcif = parse_file(&root.join(format!("bcif/{id}.bcif")), ParseOptions::default())
            .unwrap_or_else(|err| panic!("failed to parse {id}.bcif: {err}"))
            .summary();
        let mmcif = parse_file(&root.join(format!("mmcif/{id}.cif.gz")), ParseOptions::default())
            .unwrap_or_else(|err| panic!("failed to parse {id}.cif.gz: {err}"))
            .summary();
        assert_eq!(bcif.chain_count, mmcif.chain_count, "{id} chains");
        assert_eq!(bcif.residue_count, mmcif.residue_count, "{id} residues");
        assert_eq!(bcif.atom_count, mmcif.atom_count, "{id} atoms");
        assert_eq!(bcif.heteroatom_count, mmcif.heteroatom_count, "{id} heteroatoms");
    }
}
