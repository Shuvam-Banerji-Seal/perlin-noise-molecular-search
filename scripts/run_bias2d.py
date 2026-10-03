"""Formulation A1 on alanine dipeptide: does a Perlin bias accelerate barrier crossing and keep free energies right?

Usage: uv run --extra plot python scripts/run_bias2d.py [--steps 400000] [--jobs 24]
Writes results/bias2d/{surface.npy,bias2d.json,summary.md} and figures/bias2d.{pdf,png}.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # one BLAS thread per worker process; set before numpy is imported

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from perlin_mole.bias2d import KT_300, MatchedGRFBias, PerlinBias, Surface, exact_dF, relaxed_scan, run_mc

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "bias2d"
STEP = 10.0
UPDATE = 200
DRIFT = 0.25 / UPDATE  # Perlin time coordinate advances 0.25 lattice cells per bias epoch


def lag_corr(pb: PerlinBias, lag: float, n: int = 20000) -> float:
    rng = np.random.default_rng(3)
    p, q = rng.uniform(-180, 180, n), rng.uniform(-180, 180, n)
    z = rng.uniform(0, 100, n)
    a = pb.field(np.stack([(p + 180) / 360 * pb.cells, (q + 180) / 360 * pb.cells, z], -1))
    b = pb.field(np.stack([(p + 180) / 360 * pb.cells, (q + 180) / 360 * pb.cells, z + lag], -1))
    return float(np.corrcoef(a, b)[0, 1])


def make_bias(kind: str, lam: float, seed: int):
    if kind == "none" or kind == "metad":
        return None
    drift = DRIFT if kind.endswith("drift") else 0.0
    pb = PerlinBias(lam, cells=4, drift=drift, seed=seed)
    if kind.startswith("perlin"):
        return pb
    rho = lag_corr(pb, DRIFT * UPDATE) if drift else 0.0
    return MatchedGRFBias(pb, n=96, rho=rho, seed=seed)


def job(args):
    grid, kind, lam, seed, steps, walkers = args
    surf = Surface(grid, STEP)
    metad = {"width": 15.0, "height": 0.2, "gamma": 8.0} if kind == "metad" else None
    r = run_mc(
        surf, make_bias(kind, lam, seed), n_steps=steps, n_walkers=walkers, update_every=UPDATE, metad=metad, seed=seed
    )
    return {
        "kind": kind,
        "lambda": lam,
        "seed": seed,
        "transitions": r["transitions"].tolist(),
        "dF": [float(x) for x in r["dF"]],
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=400_000)
    ap.add_argument("--walkers", type=int, default=8)
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--jobs", type=int, default=24)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    sp = OUT / "surface.npy"
    if not sp.exists():
        np.save(sp, relaxed_scan(STEP, jobs=a.jobs))  # 16 start embeddings (converged; see relaxed_scan)
    grid = np.load(sp)
    surf = Surface(grid, STEP)
    ref = exact_dF(surf)
    a1 = np.arange(-180, 180, 1.0)
    barrier = {
        "phi=0": float(min(surf(np.zeros_like(a1), a1))),
        "phi=180": float(min(surf(np.full_like(a1, 180.0), a1))),
    }
    conds = [("none", 0.0), ("metad", 0.0)] + [
        (k, lam) for k in ("perlin", "grf", "perlin-drift", "grf-drift") for lam in (2.0, 4.0, 6.0)
    ]
    tasks = [(grid, k, lam, s, a.steps, a.walkers) for k, lam in conds for s in range(a.seeds)]
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        rows = list(pool.map(job, tasks))
    summary = {}
    for k, lam in conds:
        rr = [r for r in rows if r["kind"] == k and r["lambda"] == lam]
        tr = np.concatenate([r["transitions"] for r in rr]).astype(float)
        df = np.concatenate([r["dF"] for r in rr])
        fin = np.isfinite(df)
        summary[f"{k}@{lam:g}"] = {
            "kind": k,
            "lambda": lam,
            "transitions_per_1e5": float(tr.mean() / a.steps * 1e5),
            "transitions_se": float(tr.std(ddof=1) / np.sqrt(len(tr)) / a.steps * 1e5),
            "frac_walkers_crossed": float((tr > 0).mean()),
            "dF_mean": float(df[fin].mean()) if fin.any() else None,
            "dF_rmse": float(np.sqrt(np.mean((df[fin] - ref) ** 2))) if fin.any() else None,
            "frac_finite_dF": float(fin.mean()),
        }
    res = {
        "exact_dF": ref,
        "barrier_min_E": barrier,
        "kT": KT_300,
        "steps": a.steps,
        "walkers": a.walkers,
        "seeds": a.seeds,
        "summary": summary,
    }
    (OUT / "bias2d.json").write_text(json.dumps(res, indent=1))
    lines = [
        "# Perlin bias on the alanine-dipeptide (phi, psi) surface (formulation A1)",
        "",
        f"exact dF(phi>0 - phi<0) = {ref:.3f} kcal/mol; min E on phi=0: {barrier['phi=0']:.2f}, "
        f"on phi=180: {barrier['phi=180']:.2f} kcal/mol; kT = {KT_300:.3f}",
        "",
        "| condition | transitions / 1e5 steps | walkers crossed | dF mean | dF RMSE | finite dF |",
        "|---|---|---|---|---|---|",
    ]
    for k, v in summary.items():
        dm = "-" if v["dF_mean"] is None else f"{v['dF_mean']:.2f}"
        de = "-" if v["dF_rmse"] is None else f"{v['dF_rmse']:.2f}"
        lines.append(
            f"| {k} | {v['transitions_per_1e5']:.2f} ± {v['transitions_se']:.2f} | "
            f"{v['frac_walkers_crossed']:.2f} | {dm} | {de} | {v['frac_finite_dF']:.2f} |"
        )
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
