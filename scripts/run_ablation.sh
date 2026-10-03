#!/usr/bin/env bash
# Corrected ground-truth grids, then the ablation benchmark (all cores).
set -euo pipefail
cd "$(dirname "$0")/.."
rm -f results/validation/systematic.json
uv run python scripts/systematic_reference.py ala1 decane ala3 imatinib tetraglyme
ABL=reset,gauss-iid,perlin,gauss-corr,perlin-gaussmarg,gauss-perlinmarg,perlin-oct1,gauss-corr-oct1,perlin-oct4,gauss-corr-oct4,perlin-norot,gauss-corr-norot,perlin-shuffled,gauss-corr-shuffled,perlin-drift,gauss-ar1
uv run perlin-mole bench --molecules ala1,decane,tetraglyme,ala3,imatinib --methods $ABL \
  --seeds 20 --steps 500 --corr-lengths 4 --sigma 60 --out results/ablation
