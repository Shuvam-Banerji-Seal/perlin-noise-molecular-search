"""Unrelaxed rigid-rotor torsion MC: do Perlin's higher-order statistics matter when nothing is relaxed?

Usage: uv run python scripts/run_rigid_mc.py [--steps 20000 --seeds 6]   (all cores)
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # one BLAS thread per worker process; set before numpy is imported

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

from perlin_mole.benchmark import MOLECULES, Task
from perlin_mole.rigid_mc import run_rigid_task

ROOT = Path(__file__).resolve().parents[1]
METHODS = ["gauss-iid", "gauss-corr", "perlin", "perlin-gaussmarg", "gauss-perlinmarg"]
SIGMAS = [5.0, 10.0, 20.0]
MOLS = ["decane", "tetraglyme", "ala3"]
PAIRS = [
    ("perlin", "gauss-corr"),
    ("perlin-gaussmarg", "perlin"),
    ("gauss-perlinmarg", "gauss-corr"),
    ("gauss-corr", "gauss-iid"),
]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=20000)
    ap.add_argument("--seeds", type=int, default=6)
    ap.add_argument("--jobs", type=int, default=None)
    a = ap.parse_args()
    tasks = [
        Task(m, MOLECULES[m], meth, sd, a.steps, sg, 4.0, 0.596)
        for m in MOLS
        for sg in SIGMAS
        for meth in METHODS
        for sd in range(a.seeds)
    ]
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        rows = list(pool.map(run_rigid_task, tasks, chunksize=1))
    out = ROOT / "results" / "rigid_mc"
    out.mkdir(parents=True, exist_ok=True)
    summary = {}
    for m in MOLS:
        for sg in SIGMAS:
            for meth in METHODS:
                rr = [r for r in rows if r["molecule"] == m and r["sigma_deg"] == sg and r["method"] == meth]
                summary[f"{m}|{sg:g}|{meth}"] = {
                    k: float(np.mean([r[k] for r in rr]))
                    for k in ("acceptance", "esjd", "transitions_per_1000", "mean_energy")
                }
    tests = []
    for a_, b_ in PAIRS:
        for metric in ("esjd", "transitions_per_1000", "acceptance"):
            key = lambda r: (r["molecule"], r["sigma_deg"], r["seed"])  # noqa: E731
            A = {key(r): r[metric] for r in rows if r["method"] == a_}
            B = {key(r): r[metric] for r in rows if r["method"] == b_}
            ks = sorted(set(A) & set(B))
            ratio = np.array([np.log((A[k] + 1e-9) / (B[k] + 1e-9)) for k in ks])
            d = np.array([A[k] - B[k] for k in ks])
            p = float(wilcoxon(d).pvalue) if np.any(d != 0) else 1.0
            tests.append(
                {
                    "a": a_,
                    "b": b_,
                    "metric": metric,
                    "n": len(ks),
                    "median_log_ratio": float(np.median(ratio)),
                    "wins": int((d > 0).sum()),
                    "losses": int((d < 0).sum()),
                    "p": p,
                }
            )
    ps = [t["p"] for t in tests]
    order = np.argsort(ps)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (len(ps) - rank) * ps[i]))
        tests[i]["p_holm"] = run
    json.dump({"rows": rows, "summary": summary, "tests": tests}, open(out / "rigid_mc.json", "w"), indent=1)
    lines = [
        "# Unrelaxed rigid-rotor torsion MC (300 K)",
        "",
        "| comparison | metric | median log(a/b) | a better / b better | p (Holm) |",
        "|---|---|---|---|---|",
    ]
    for t in tests:
        lines.append(
            f"| {t['a']} vs {t['b']} | {t['metric']} | {t['median_log_ratio']:+.3f} | {t['wins']} / {t['losses']} "
            f"| {t['p']:.3g} ({t['p_holm']:.2g}) |"
        )
    lines += [
        "",
        "| molecule | sigma | method | acceptance | ESJD (rad^2) | transitions/1000 | <E> |",
        "|---|---|---|---|---|---|---|",
    ]
    for k, v in summary.items():
        m, sg, meth = k.split("|")
        lines.append(
            f"| {m} | {sg} | {meth} | {v['acceptance']:.3f} | {v['esjd']:.4f} | {v['transitions_per_1000']:.2f} "
            f"| {v['mean_energy']:.2f} |"
        )
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:20]))
