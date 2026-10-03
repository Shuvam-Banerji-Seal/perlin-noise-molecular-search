#!/usr/bin/env bash
# "Hunt" for a regime where Perlin helps: larger/floppier molecules, two kick sizes,
# coarse-to-fine (anneal) and time-coherent (drift) variants, each with its matched Gaussian.
set -euo pipefail
cd "$(dirname "$0")/.."
METHODS=reset,gauss-iid,gauss-corr,perlin,perlin-anneal,gauss-anneal,perlin-drift,gauss-ar1
for sg in ${SIGMAS:-60 120}; do
  uv run python - "$sg" "$METHODS" <<'PY'
import sys
from perlin_mole import benchmark
from perlin_mole.benchmark import run_benchmark
HUNT = {
    "eicosane": "CCCCCCCCCCCCCCCCCCCC",
    "hexaglyme": "COCCOCCOCCOCCOCCOCCOC",
    "ala6": "CC(=O)N" + "[C@@H](C)C(=O)N" * 6 + "C",
    "methyl-oleate": "CCCCCCCC/C=C\\CCCCCCCC(=O)OC",
}
benchmark.KNOWN.update(HUNT)
sg, methods = float(sys.argv[1]), sys.argv[2].split(",")
run_benchmark(list(HUNT), methods, 10, 500, sg, (8.0,), 1.0, 3.0, None, f"results/hunt/sigma{int(sg)}")
PY
done
