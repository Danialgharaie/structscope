# structscope

`structscope` is a Rust-native structural bioinformatics toolkit for canonical protein structure parsing, graph-native representations, reproducible feature extraction, and analytical outputs.

This repository currently contains a bootstrap implementation with:

- workspace scaffolding for all planned crates
- crate-backed PDB, mmCIF, and BinaryCIF parsing with gzip input support
- canonical structure normalization
- residue, atom, and interface graph construction (GraphML export)
- structural primitives: solvent accessible surface area (Shrake-Rupley), DSSP-style secondary structure, backbone dihedrals, optimal superposition/RMSD (Kabsch), and typed interactions (disulfides, salt bridges, hydrogen bonds)
- basic and graph-derived feature extraction
- JSONL and Parquet feature export
- DuckDB-backed SQL querying over feature Parquet (build with `--features duckdb`)
- optional SQLite/JSONL provenance
- CLI entrypoints for parse, featurize, graph, query, rmsd, and provenance

Querying is gated behind a Cargo feature because it bundles DuckDB:

```
cargo build -p structscope-cli --features duckdb
structscope query <features.parquet|out-dir> --sql "SELECT * FROM features"
```

Feature records are exposed to SQL as a `features` table.

Current limitations:

- the eBPF guard crate is scaffolded only
