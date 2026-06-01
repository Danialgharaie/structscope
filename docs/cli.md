# CLI Usage

`structscope <command> [options]`. All commands accept PDB, mmCIF, and
BinaryCIF inputs, including gzip-compressed variants (`.pdb.gz`, `.cif.gz`,
`.bcif.gz`).

## parse

Summarise a structure (chains, residues, atoms, heteroatoms, ligands).

```
structscope parse 1nkd.cif.gz
structscope parse 1nkd.bcif --format json
```

## featurize

Compute structure-level features and write them to an output directory
(JSONL + Parquet). Accepts a single file or a directory to batch-process.

```
structscope featurize 1nkd.cif.gz --out ./out
structscope featurize ./structures --out ./out --provenance
```

Emitted features include counts (atoms, residues, chains, ligands), graph
metrics (contacts, density, clustering), geometry (radius of gyration, SASA),
secondary-structure composition, typed-interaction counts, and
buried/exposed residue counts.

## residues

Emit one JSON record per residue (SASA, RSA, secondary structure,
phi/psi/omega) as JSONL, to stdout or a file.

```
structscope residues 1nkd.cif.gz
structscope residues 1nkd.cif.gz --out residues.jsonl
```

## rmsd

Optimal-superposition RMSD between two structures.

```
# Equal-length structures, matched by atom order:
structscope rmsd ref.pdb mobile.pdb --atoms ca       # or: backbone, all

# Different-length but related structures (sequence-aligned CA atoms):
structscope rmsd ref.pdb mobile.pdb --align

# Partial or domain-level overlap (local Smith-Waterman alignment):
structscope rmsd ref.pdb fragment.pdb --local
```

Without `--align`, the two selections must have equal atom counts; the error
message hints at `--align` when they differ.

## graph

Export a residue, atom, or interface contact graph as GraphML.

```
structscope graph 1nkd.cif.gz --graph-type residue --out graph.graphml
structscope graph complex.cif.gz --graph-type interface
```

## query

Run SQL over a feature Parquet file or featurize output directory. Requires a
build with the `duckdb` feature.

```
cargo build -p structscope-cli --features duckdb
structscope query ./out --sql "SELECT structure_id, sasa_total FROM features"
```

## provenance

Inspect a provenance SQLite database produced by `featurize --provenance`.

```
structscope provenance ./out/run.sqlite
```
