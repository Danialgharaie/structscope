# Interface golden fixtures

Golden expected values for structscope protein–protein interface metrics (Slice B).

## Synthetic fixtures

`synthetic_dimer.json` — minimal two-chain CA dimer used in unit tests. Regenerate with:

```bash
cargo test -p structscope-features write_synthetic_golden_fixture -- --ignored --nocapture
```

## CCP4/PISA validated fixtures (offline)

`ccp4_dimer.json` (future) — reference values from offline CCP4 `sc` / PISA BSA comparison on a real dimer PDB.

Regeneration steps:

1. Run `structscope interfaces` on the fixture PDB with default distance params (8/5/5).
2. Compare BSA and interface area to PISA; compare shape complementarity to CCP4 `sc`.
3. Update the JSON expected values.

Tolerances when validating against CCP4/PISA reference:

- BSA: ±2%
- Interface area: ±2%
- Shape complementarity (SC): ±0.05
