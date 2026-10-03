"""Pre-registered confirmatory low-temperature test (see results/lowT_confirm/PREREGISTRATION.md)."""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json
import pickle
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

from perlin_mole import ablation
from perlin_mole.benchmark import MOLECULES, Task, run_task

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "lowT_confirm"
MOLS = ["decane", "tetraglyme", "ala3"]
SEEDS = range(10, 40)

if __name__ == "__main__":
    tasks = [
        Task(m, MOLECULES[m], meth, s, 500, 60.0, 4.0, 0.3)
        for m in MOLS
        for meth in ("perlin", "gauss-corr")
        for s in SEEDS
    ]
    with ProcessPoolExecutor() as pool:
        runs = list(pool.map(run_task, tasks, chunksize=1))
    pickle.dump(runs, open(OUT / "raw.pkl", "wb"))
    ablation.COMPARISONS[:] = [("perlin(L=4)", "gauss-corr(L=4)", "low-T confirmation")]
    rep = ablation.analyze(OUT, ROOT / "results" / "validation", OUT)
    rec = {(r["molecule"], r["seed"], r["label"]): r["recall_final"] for r in rep["rows"]}
    keys = sorted({(m, s) for m, s, _ in rec})
    d = np.array([rec[(m, s, "perlin(L=4)")] - rec[(m, s, "gauss-corr(L=4)")] for m, s in keys])
    res = {
        "n_pairs": len(d),
        "mean_diff": float(d.mean()),
        "wins": int((d > 0).sum()),
        "losses": int((d < 0).sum()),
        "wilcoxon_p_two_sided": float(wilcoxon(d).pvalue),
        "per_molecule_mean_diff": {
            m: float(np.mean([x for (mm, _), x in zip(keys, d, strict=True) if mm == m])) for m in MOLS
        },
        "bootstrap_ci95": next(t["ci95"] for t in rep["tests"] if t["metric"] == "recall_final"),
    }
    (OUT / "result.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
