#!/usr/bin/env bash
# Regenerate every result and figure of the study from scratch. Uses all cores.
# On a 24-core node + one A100: roughly 8-10 h in total (MMFF relaxations dominate; the tetraglyme grid is
# 531k relaxations; the LJ38 confirmation/replication is ~50 CPU-h; DFT ~1-2 GPU-h).
set -euo pipefail
cd "$(dirname "$0")/.."
uv sync --extra plot --extra gpu --extra dft

# formulation A: cost and gradient correlation of a full-dimensional Perlin PES
uv run perlin-mole landscape --max-dim 14 | tee results/landscape.txt

# S1: flexible organics at three kick sizes (raw.pkl kept for pooling)
MOLS=decane,tetraglyme,ala3,imatinib
uv run perlin-mole bench --molecules $MOLS --seeds 10 --steps 500 --corr-lengths 2,4,8 --sigma 60 --out results/sigma60
for sg in 30 120; do
  uv run perlin-mole bench --molecules $MOLS --methods gauss-iid,gauss-corr,perlin --seeds 10 --steps 500 \
    --corr-lengths 4,8 --sigma $sg --out results/sigma$sg
done

# S2-S4: ground truth (alanine dipeptide, exhaustive grids) and the 16-variant ablation
uv run python scripts/validate_alanine.py
./scripts/run_ablation.sh                      # grids + ablation benchmark
uv run python -m perlin_mole.ablation results/ablation

# S5: predicted vs real (experimental structures are in data/experimental; re-fetch with fetch_experimental.py)
uv run python scripts/run_experimental.py --jobs "$(nproc)"

# S6/S7: isomers and enantiomers
uv run python scripts/run_isomers.py conformational enantiomers diastereomers constitutional --seeds 10 --steps 300

# S8: Lennard-Jones clusters
uv run python scripts/run_clusters.py --jobs "$(nproc)" --seeds 30

# S9-S11: regime hunt, temperature sweep, Cartesian kicks, 2-D Perlin bias
./scripts/run_hunt.sh
uv run python scripts/run_hunt_kt.py
uv run python scripts/run_cartesian.py --seeds 10 --steps 400
uv run python scripts/run_bias2d.py --jobs "$(nproc)"

# pre-registered confirmations and replication (see PREREGISTRATION.md in each results folder)
uv run python scripts/run_lj38_confirm.py "$(nproc)"                                  # LJ38, seeds 30-129
LJ38_OUT=clusters_replication LJ38_SEEDS=130-430 LJ38_METHODS=perlin@2,gauss-field@2,iid-perlinmarg \
  uv run python scripts/run_lj38_confirm.py "$(nproc)"                                # replication + mechanism arm
uv run python scripts/run_lowT_confirm.py                                             # low-temperature hint
uv run python scripts/run_isomers.py enantiomers --guard --seeds 10 --steps 300      # stereo-guarded enantiomers

# pooled discoveries (SDF export, common-reference coverage, recall vs grid)
uv run python -m perlin_mole.export results/sigma30 results/sigma60 results/sigma120 --out results/discoveries

# second and third judges on the GPU: ANI-2x (TorchANI) and DFT (GPU4PySCF, B3LYP-D3(BJ))
export TORCHANI_NO_WARN_EXTENSIONS=1
uv run --extra gpu --extra plot python scripts/run_ml_judge.py
uv run --extra gpu python scripts/run_ani_recall.py
uv run --extra dft --extra gpu python scripts/run_dft.py d1 d2
uv run --extra dft --extra gpu python scripts/run_dft.py d3   # def2-TZVP global-minimum arbitration

# figures and animations from the saved results
uv run python scripts/rescore_ground_truth.py   # corrected ground truth; low-T sensitivity analysis
uv run --extra plot python scripts/make_figures.py
uv run --extra plot python scripts/make_animations.py   # README animations (docs/animations/)
