# structscope

`structscope` is a Rust-native structural bioinformatics toolkit for canonical protein structure parsing, graph-native representations, reproducible feature extraction, and analytical outputs.

This repository currently contains a bootstrap implementation with:

- workspace scaffolding for all planned crates
- plain-text PDB and mmCIF parsing
- canonical structure normalization
- residue-graph construction
- basic and graph-derived feature extraction
- JSONL feature export
- optional SQLite/JSONL provenance
- CLI entrypoints for parse, featurize, graph, query, and provenance

Current limitations:

- transparent gzip and BinaryCIF are not implemented yet
- Parquet writing is not implemented in this bootstrap slice
- DuckDB-backed querying is a stub until the `duckdb` integration is added
- the eBPF guard crate is scaffolded only
