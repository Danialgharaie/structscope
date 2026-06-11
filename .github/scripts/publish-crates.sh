#!/usr/bin/env bash
set -euo pipefail

# Publish workspace crates to crates.io in dependency order.
# Requires CARGO_REGISTRY_TOKEN (or login via `cargo login`).

crates=(
  structscope-core
  structscope-events
  structscope-agent
  structscope-graphs
  structscope-features
  structscope-store
  structscope-provenance
  structscope-cli
)

for pkg in "${crates[@]}"; do
  echo "Publishing ${pkg}..."
  cargo publish -p "${pkg}"
  if [[ "${pkg}" != "structscope-cli" ]]; then
    echo "Waiting for crates.io index to update..."
    sleep 45
  fi
done

echo "All crates published."
