# Changelog

All notable changes to this project are documented here. This project follows
[Keep a Changelog](https://keepachangelog.com/) conventions.

## [Unreleased]

## [0.4.0] - 2026-06-11

### Added

- Multi-structure compare: pairwise RMSD matrix (CA atoms with sequence
  alignment by default) and numeric feature deltas against a chosen reference.
- CLI command:
  - `compare <input>` — compare two or more structures from a file or
    directory; prints JSON to stdout or writes `matrix.json` + `deltas.jsonl`
    (or `matrix.csv` + `deltas.csv` with `--format csv`) to `--out`.
- Reference selection (first match wins): `--reference`, `--reference-by
  min:field|max:field`, `--auto-reference` (lowest Ramachandran outliers,
  clashes, missing backbone), else first input.
- `--delta-fields` to restrict feature delta columns; `--atoms`, `--align`,
  and `--local` for RMSD correspondence (same semantics as `rmsd`).
- Ligand, interface, and quality flags shared with `featurize`.

## [0.3.1] - 2026-06-11

### Added

- Golden fixtures and regression tests for structure quality metrics
  (`tests/fixtures/quality/`, `quality_golden.rs`).

## [0.3.0] - 2026-06-11

### Added

- Structure quality metrics: MolProbity-style Ramachandran classification
  (favored / allowed / outlier, with Gly/Pro regions), steric clash detection
  (heavy atoms, configurable VdW overlap), and missing backbone atom checks
  (N, CA, C, O) over canonical and common variant residues.
- CLI command:
  - `quality <input>` — per-residue quality records as JSONL (problems only
    by default; `--all-residues` for full output).
- Structure-level quality aggregates in `featurize`: `quality_residue_count`,
  Ramachandran counts, `clash_pair_count`, and `missing_backbone_residue_count`.
- `--clash-overlap` flag on `featurize` and `quality` (default `0.4` Å).

## [0.2.0] - 2026-06-11

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
- Protein–protein interface metrics: buried surface area (BSA), interface patch
  area, and Lawrence–Colman shape complementarity per contacting chain pair.
- CLI command:
  - `interfaces <input>` — per chain-pair interface features (BSA, area, shape
    complementarity, contact and residue counts) as JSONL.
- Structure-level interface aggregates and largest-interface fields in
  `featurize`: `interface_pair_count`, total and max BSA/area/SC, and
  largest-interface chain IDs.
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
