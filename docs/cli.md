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
structscope featurize ./structures --out ./out --provenance -j 4
structscope featurize ./structures --out ./out --ligand-exclude SO4,PO4
structscope featurize dimer.pdb --out ./out --interface-distance 8.0
```

Optional ligand flags (also on `ligands`):

- `--ligand-exclude RES[,RES...]` — add residue names to the default denylist
- `--ligand-include RES[,RES...]` — allowlist mode; only these hetero residues count
- `--binding-distance <Å>` — binding-site cutoff (default `5.0`)

Optional interface flags (also on `interfaces`):

- `--interface-distance <Å>` — chain-pair contact cutoff (default `8.0`)
- `--interface-area-distance <Å>` — interface patch area cutoff (default `5.0`)
- `--interface-sc-distance <Å>` — shape complementarity surface cutoff (default `5.0`)

Emitted features include counts (atoms, residues, chains, ligands), graph
metrics (contacts, density, clustering), geometry (radius of gyration, SASA),
secondary-structure composition, typed-interaction counts, protein–ligand
interaction counts, binding-site residue count, ligand SASA, buried/exposed
residue counts, and protein–protein interface summaries (pair count, total
and max BSA/area/SC, largest-interface chain IDs).

## ligands

Emit one JSON record per filtered ligand (SASA, binding-site residues,
interaction counts).

```
structscope ligands 1nkd.cif.gz
structscope ligands complex.cif.gz --out ligands.jsonl
structscope ligands complex.cif.gz --ligand-include HEM,NAG --binding-distance 4.0
```

## interfaces

Emit one JSON record per contacting chain pair (BSA, interface patch area,
shape complementarity, contact and residue counts).

```
structscope interfaces 1nkd.cif.gz
structscope interfaces dimer.pdb --out interfaces.jsonl
structscope interfaces dimer.pdb --interface-distance 8.0 --interface-area-distance 5.0
```

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

Export a residue, atom, or interface contact graph. Supported formats are `graphml` (default), `gml`, and `json`. Chemical and geometric interactions (disulfides, salt bridges, hydrogen bonds, cation-pi, pi-pi, and hydrophobic contacts) are automatically resolved and embedded as prioritized edges in residue and interface graphs.

```
structscope graph 1nkd.cif.gz --graph-type residue --format gml
structscope graph complex.cif.gz --graph-type interface --format json
structscope graph 1nkd.cif.gz --graph-type residue --out graph.graphml
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
