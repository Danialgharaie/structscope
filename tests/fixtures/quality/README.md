# Quality golden fixtures

Golden expected values for structscope structure quality metrics (Slice C).

## Synthetic fixtures

| File | Scenario |
|------|----------|
| `synthetic_panel.json` | Clash between distant residues (1 vs 3), missing backbone O on residue 4 |
| `synthetic_dipeptide.json` | Clean extended dipeptide; no problems, sequential neighbors not clashed |

Regenerate with:

```bash
cargo test -p structscope-features write_synthetic_panel_golden_fixture -- --ignored --nocapture
cargo test -p structscope-features write_synthetic_dipeptide_golden_fixture -- --ignored --nocapture
```

Run regression tests:

```bash
cargo test -p structscope-features quality_golden -- --nocapture
```

## MolProbity validated fixtures (offline)

`molprobity_reference.json` (future) — reference Ramachandran and clash counts from offline MolProbity on a real PDB.

Regeneration steps:

1. Run `structscope quality <pdb> --all-residues` with default `--clash-overlap 0.4`.
2. Compare Ramachandran outlier counts and clash totals to MolProbity report.
3. Update the JSON expected summary and problem records.

Suggested tolerances when validating against MolProbity reference:

- Summary counts: exact match for synthetic fixtures; ±1–2% on clash pairs for real structures if VdW radii differ slightly
- Ramachandran class per residue: exact match (same polygon definitions)
