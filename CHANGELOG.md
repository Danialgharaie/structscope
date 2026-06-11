# Changelog

All notable changes to this project are documented here. This project follows
[Keep a Changelog](https://keepachangelog.com/) conventions.

## [Unreleased]

### Added

- Parallel execution support for `structscope featurize` using Rayon, allowing high-performance concurrent parsing and feature extraction controlled via a new `--jobs` / `-j` CLI flag.
- B-factor (temperature factor) support: extended core `Atom` model, parsed B-factors from PDB, mmCIF, and BinaryCIF, and computed structure-level (`bfactor_mean`, `bfactor_std`, `bfactor_min`, `bfactor_max`) and residue-level statistics.
- Advanced geometric interaction detectors: cation-pi, parallel and perpendicular pi-pi stacking, and hydrophobic carbon-carbon contacts.
- Enhanced contact graphs: integrated chemical/geometric interactions as prioritized edges overlaying standard distance contacts.
- New contact graph formats: added custom exporters for GML and node-link JSON formats, with automatic file extension matching.
- Thread-safe background event and provenance logging architecture using a message-passing channel (`mpsc`) to safely write events to SQLite and JSONL from worker threads.
- BinaryCIF (`.bcif` / `.bcif.gz`) parsing: a hand-written MessagePack-based
  decoder covering all seven column encodings (ByteArray, IntegerPacking,
  RunLength, Delta, FixedPoint, IntervalQuantization, StringArray) and value
  masks. structscope now ingests PDB, mmCIF, and BinaryCIF.
- Structural primitives computed directly from coordinates:
  - Solvent accessible surface area (Shrake-Rupley), per-atom and total.
  - Relative solvent accessibility (RSA), per residue, normalised by
    residue-type maxima (Tien et al. 2013).
  - DSSP-style secondary structure (Kabsch-Sander hydrogen bonds).
  - Backbone dihedrals (phi/psi/omega).
  - Optimal superposition and RMSD (quaternion/Kabsch).
  - Typed interactions: disulfides, salt bridges, hydrogen bonds.
- Sequence alignment primitive (Needleman-Wunsch) for residue correspondence.
- CLI commands:
  - `rmsd <reference> <mobile>` — optimal-superposition RMSD over matched
    atoms (`--atoms ca|backbone|all`), with `--align` for sequence-based
    correspondence between structures of different lengths.
  - `residues <input>` — per-residue features (SASA, RSA, secondary
    structure, dihedrals) as JSONL.
  - `ligands <input>` — per-ligand features (SASA, binding-site residues,
    protein–ligand interaction counts) as JSONL.
- Protein–ligand features: configurable ligand filter (default excludes water
  and common ions), structure-level protein–ligand interaction counts,
  binding-site residue count, and ligand SASA.
- Structure-level features: `sasa_total`, `helix/strand/coil_residue_count`,
  `disulfide/salt_bridge/hydrogen_bond_count`, `buried/exposed_residue_count`,
  `ligand_sasa_total`, `ligand_sasa_mean`, `binding_site_residue_count`,
  `protein_ligand_hbond_count`, `protein_ligand_salt_bridge_count`,
  `protein_ligand_hydrophobic_count`, `protein_ligand_contact_count`.

### Changed

- **Breaking:** `ligand_count` in `featurize` now uses the filtered ligand
  definition (hetero residues minus the default denylist and CLI overrides),
  not the raw hetero residue count.

### Notes

- All primitives emit raw quantities; downstream interpretation is left to the
  user. Per-residue and per-atom detail is exposed as library functions.
