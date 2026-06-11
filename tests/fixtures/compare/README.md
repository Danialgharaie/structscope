# Compare golden fixtures

Golden expected values for structscope multi-structure comparison (Slice D).

## Synthetic fixtures

| File | Scenario |
|------|----------|
| `synthetic_triplet.json` | Three mini CA-only structures with x-offsets 0.0 / 0.5 / 1.0 Å; first-input reference; `sasa_total` feature deltas |

Regenerate with:

```bash
cargo test -p structscope-features write_synthetic_triplet_golden_fixture -- --ignored --nocapture
```

Run regression tests:

```bash
cargo test -p structscope-features compare_golden -- --nocapture
```

Suggested tolerances:

- RMSD matrix cells: absolute ±1e-6 for synthetic fixtures
- Feature deltas: exact match for synthetic numeric fields
