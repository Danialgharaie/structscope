use std::path::PathBuf;
use structscope_core::{parse_file, ParseOptions};
use structscope_features::dihedral::backbone_dihedrals;

fn corpus_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/corpus/rcsb")
}

/// On a real protein, all defined torsions must lie in [-180, 180], and an
/// all-alpha fold (1nkd) should have a clearly negative mean phi (helical region).
#[test]
fn dihedrals_valid_range_and_helical_phi() {
    let s = parse_file(&corpus_root().join("mmcif/1nkd.cif.gz"), ParseOptions::default()).unwrap();
    let d = backbone_dihedrals(&s);
    assert!(d.len() > 100);

    // phi is a circular quantity, so count residues whose phi lies in the helical
    // half-plane (-180,0) rather than averaging. An all-alpha fold is dominated by it.
    let (mut neg_phi, mut total_phi) = (0u32, 0u32);
    for r in &d {
        for a in [r.phi, r.psi, r.omega].into_iter().flatten() {
            assert!((-180.0..=180.0).contains(&a), "angle out of range: {a}");
        }
        if let Some(phi) = r.phi {
            total_phi += 1;
            if (-180.0..0.0).contains(&phi) {
                neg_phi += 1;
            }
        }
    }
    assert!(
        neg_phi * 2 > total_phi,
        "all-alpha fold: most phi should be negative ({neg_phi}/{total_phi})"
    );
}