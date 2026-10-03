#!/usr/bin/env bash
# One command to install, test and try the package (a few minutes on a laptop):
#   1. install Python 3.12 and the locked dependencies with uv
#   2. run the unit tests
#   3. run a 10-second basin-hopping benchmark (Perlin vs matched Gaussian vs i.i.d. vs ETKDG on n-decane)
#   4. regenerate the figures from the shipped results (figures that need the large raw files are skipped
#      until scripts/fetch_raw_results.sh has downloaded them)
# Needs uv: https://docs.astral.sh/uv/getting-started/installation/
set -euo pipefail
cd "$(dirname "$0")/.."
uv sync --extra plot
uv run pytest -q
uv run perlin-mole bench --molecules decane --seeds 4 --steps 100 --corr-lengths 4 --out runs/demo
uv run --extra plot python scripts/make_figures.py
echo
echo "demo benchmark: runs/demo/summary.md and runs/demo/coverage.png; figures: figures/"
