# Smoke Checklist

- `structscope parse <file>` prints a deterministic summary.
- `structscope featurize <dir> --out <dir>` writes `features.jsonl` and `manifest.json`.
- `structscope graph <file> --out graph.graphml` writes a residue graph.
- `structscope provenance run.sqlite` lists recorded runs when provenance is enabled.
- `structscope query ...` currently reports the missing DuckDB integration explicitly.
