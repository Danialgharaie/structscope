# Architecture

The first implementation slice is CLI-first and keeps crate boundaries aligned with the long-term design:

- `structscope-core`: data model and parsers
- `structscope-graphs`: graph builders and export
- `structscope-features`: scientific features
- `structscope-store`: persisted outputs and query adapter boundary
- `structscope-events`: structured events
- `structscope-provenance`: optional lineage capture
- `structscope-agent`: optional eBPF integration boundary
- `structscope-cli`: orchestration

The scientific path is independent of provenance and eBPF so batch processing can remain portable.
