#!/usr/bin/env bash
# Download the 13 large raw result files (>= 10 MB; kept as GitHub release assets, not in git) into results/
# and verify them against results/RAW_RESULTS.sha256. Files that are already present and correct are skipped.
#   bash scripts/fetch_raw_results.sh
# Needs curl, xz and sha256sum. The JSON/Markdown summaries in results/ do not need these files; they are only
# required to re-run analyses that pool the raw per-run data (e.g. perlin-mole rescore, perlin_mole.export).
set -euo pipefail
cd "$(dirname "$0")/.."
REPO="Shuvam-Banerji-Seal/perlin-noise-molecular-search"
TAG="${RAW_RESULTS_TAG:-v1.0.0}"
BASE="https://github.com/$REPO/releases/download/$TAG"
while read -r sum path; do
  if [ -f "$path" ] && echo "$sum  $path" | sha256sum --check --status -; then
    echo "ok        $path"
    continue
  fi
  asset="$(echo "$path" | sed 's#/#__#g').xz"
  mkdir -p "$(dirname "$path")"
  echo "download  $path"
  curl -fL --retry 3 --progress-bar -o "$path.xz" "$BASE/$asset"
  xz -d -f "$path.xz"
  echo "$sum  $path" | sha256sum --check -
done < results/RAW_RESULTS.sha256
