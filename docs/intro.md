# structscope

<p align="center">
  <img src="../assets/structscope-logo.svg" alt="structscope" width="360">
</p>

> **Work in progress** — APIs, CLI flags, and output schemas may change between
> releases. Pin a version tag for reproducible workflows.

`structscope` is a Rust-native structural bioinformatics toolkit for canonical
protein structure parsing, graph-native representations, reproducible feature
extraction, and analytical outputs.

It parses PDB, mmCIF, and BinaryCIF (with gzip support) into a canonical model,
builds residue/atom/interface graphs, and computes raw structural primitives:
solvent accessible surface area, relative accessibility, DSSP-style secondary
structure, backbone dihedrals, optimal superposition/RMSD, and typed
interactions. All primitives emit raw quantities; downstream interpretation is
left to the user.

See [CLI Usage](cli.md) to get started, [Architecture](architecture.md) for the
crate layout, and the [Changelog](changelog.md) for recent additions.
