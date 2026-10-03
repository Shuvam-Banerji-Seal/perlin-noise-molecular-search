"""Literal Cartesian Perlin kicks on molecules, scored against grid ground truth. All cores.

Usage: uv run python scripts/run_cartesian.py [--seeds 10 --steps 400 --s 50]
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # one BLAS thread per worker process; set before numpy is imported

import argparse
import pickle
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from perlin_mole import ablation
from perlin_mole.benchmark import KNOWN, MOLECULES, Task
from perlin_mole.cartesian import run_cartesian_task

ROOT = Path(__file__).resolve().parents[1]
METHODS = ["cart-iid", "cart-perlin@2", "cart-gauss-field@2", "cart-perlin@4", "cart-gauss-field@4"]

ablation.COMPARISONS[:] = [
    ("cart-perlin@2", "cart-gauss-field@2", "Cartesian Perlin vs matched field (L=2 A)"),
    ("cart-perlin@4", "cart-gauss-field@4", "Cartesian Perlin vs matched field (L=4 A)"),
    ("cart-perlin@2", "cart-iid", "Cartesian Perlin (L=2 A) vs i.i.d."),
    ("cart-gauss-field@2", "cart-iid", "Cartesian correlated field vs i.i.d."),
]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--steps", type=int, default=400)
    ap.add_argument("--s", type=float, default=50.0, help="displacement std in 0.01 A")
    ap.add_argument("--molecules", default="ala1,decane,ala3,imatinib")
    a = ap.parse_args()
    out = ROOT / "results" / "cartesian"
    out.mkdir(parents=True, exist_ok=True)
    smi = {**MOLECULES, **KNOWN}
    tasks = [
        Task(m, smi[m], meth, sd, a.steps, a.s, 0.0, 1.0)
        for m in a.molecules.split(",")
        for meth in METHODS
        for sd in range(a.seeds)
    ]
    tasks.sort(key=lambda t: -len(t.smiles))
    with ProcessPoolExecutor() as pool:
        runs = list(pool.map(run_cartesian_task, tasks, chunksize=1))
    pickle.dump(runs, open(out / "raw.pkl", "wb"))
    for m in a.molecules.split(","):
        for meth in METHODS:
            rr = [r for r in runs if r["task"]["molecule"] == m and r["task"]["method"] == meth]
            print(
                m,
                meth,
                "integrity",
                round(1 - np.mean([r["broken"] for r in rr]) / a.steps, 3),
                "acc",
                round(np.mean([r["acceptance"] for r in rr]), 3),
            )
    rep = ablation.analyze(out, ROOT / "results" / "validation", out)
    print(ablation.format_markdown(rep))
