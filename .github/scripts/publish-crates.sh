#!/usr/bin/env bash
set -euo pipefail

# Publish workspace crates to crates.io in dependency order.
# Requires CARGO_REGISTRY_TOKEN (set via `cargo login` or GitHub Actions secret).

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

VERSION="$(grep '^version' Cargo.toml | head -1 | sed -E 's/.*"([^"]+)".*/\1/')"
PUBLISH_SLEEP="${PUBLISH_SLEEP:-90}"

publish_crate() {
  local pkg="$1"
  local attempt
  for attempt in 1 2 3 4 5; do
    set +e
    output="$(cargo publish -p "${pkg}" 2>&1)"
    status=$?
    set -e
    echo "${output}"
    if [[ ${status} -eq 0 ]]; then
      return 0
    fi
    if echo "${output}" | grep -q "already exists on crates.io index"; then
      echo "${pkg} ${VERSION} already published; skipping."
      return 0
    fi
    if echo "${output}" | grep -q "429 Too Many Requests"; then
      wait_hint="$(echo "${output}" | sed -n 's/.*try again after \(.*\) and see.*/\1/p' | head -1)"
      if [[ -n "${wait_hint}" ]]; then
        wait_until="$(date -d "${wait_hint}" +%s 2>/dev/null || true)"
        now="$(date +%s)"
        if [[ -n "${wait_until}" && ${wait_until} -gt ${now} ]]; then
          sleep_for=$((wait_until - now + 5))
          echo "Rate limited; sleeping ${sleep_for}s until ${wait_hint}..."
          sleep "${sleep_for}"
          continue
        fi
      fi
      backoff=$((attempt * 60))
      echo "Rate limited; backing off ${backoff}s (attempt ${attempt}/5)..."
      sleep "${backoff}"
      continue
    fi
    return "${status}"
  done
  echo "Failed to publish ${pkg} after retries."
  return 1
}

for pkg in "${crates[@]}"; do
  echo "Publishing ${pkg} ${VERSION}..."
  publish_crate "${pkg}"
  if [[ "${pkg}" != "structscope-cli" ]]; then
    echo "Waiting ${PUBLISH_SLEEP}s for crates.io index..."
    sleep "${PUBLISH_SLEEP}"
  fi
done

echo "All crates published (or already present)."
