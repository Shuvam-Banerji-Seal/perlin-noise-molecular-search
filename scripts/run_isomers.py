"""Isomer / enantiomer / conformational-isomer suites (all cores).

Usage: uv run python scripts/run_isomers.py [--seeds 10 --steps 300] [suite ...]
Suites: enantiomers diastereomers constitutional conformational
Outputs: results/isomers/<suite>/ (raw.pkl, report.json, ablation.json/md), results/isomers/grids/, isomers.json/md
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from perlin_mole import ablation, benchmark
from perlin_mole.benchmark import run_benchmark
from perlin_mole.isomers import (
    CONFORMATIONAL,
    CONSTITUTIONAL,
    DIASTEREOMERS,
    ENANTIOMERS,
    atom_order_map,
    cip_labels,
    mirror_match,
    mirror_smiles,
    stereo_fidelity,
)
from perlin_mole.search import SymmetryRMSD
from perlin_mole.validate import systematic_search

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "isomers"
GRIDS = OUT / "grids"
METHODS = ["reset", "gauss-iid", "gauss-corr", "perlin"]


def grid_values(k: int) -> list[float]:
    """Finest grid that keeps the number of starts manageable."""
    for step in (15, 30, 60, 120):
        if (360 // step) ** k <= 60000:
            return [float(a) for a in range(-180, 180, step)]
    return [60.0, 180.0, 300.0]


def build_grid(name: str, smiles: str) -> dict:
    path = GRIDS / f"{name}_systematic.npz"
    from perlin_mole.chem import TorsionSystem

    k = TorsionSystem.from_smiles(smiles).n_dof
    vals = grid_values(k) if k <= 6 else [60.0, 180.0, 300.0]
    if not path.exists():
        e, heavy, full, n, n_unconv = systematic_search(smiles, vals, return_unconverged=True)
        np.savez_compressed(path, energies=e, heavy=heavy, full=full)
        meta = {"n_starts": n, "n_unconverged_discarded": n_unconv}
    else:
        meta = {}
    g = np.load(path)
    return {
        "k": k,
        "grid_step": vals[1] - vals[0],
        "n_unique": int(len(g["energies"])),
        "e_min": float(g["energies"][0]),
        "energies": g["energies"].tolist()[:12],
        **meta,
    }


def suite_molecules(suite: str) -> dict[str, str]:
    if suite == "enantiomers":
        out = {}
        for n, s in ENANTIOMERS.items():
            lab = "".join(cip_labels(s).values())
            out[f"{n}-{lab}"] = s
            out[f"{n}-{''.join(cip_labels(mirror_smiles(s)).values())}"] = mirror_smiles(s)
        return out
    return {"diastereomers": DIASTEREOMERS, "constitutional": CONSTITUTIONAL, "conformational": CONFORMATIONAL}[suite]


def run_suite(suite: str, seeds: int, steps: int, guard: bool = False) -> dict:
    mols = suite_molecules(suite)
    benchmark.KNOWN.update(mols)  # make names resolvable in workers and in the analysis
    grids = {n: build_grid(n, s) for n, s in mols.items()}
    res: dict = {"molecules": {n: {"smiles": s, "grid": grids[n]} for n, s in mols.items()}}
    if suite == "conformational":
        return res  # judge validation only: the grid is the complete answer for 1-2 torsions
    out = OUT / (suite + ("-guarded" if guard else ""))
    run_benchmark(list(mols), METHODS, seeds, steps, 60.0, (4.0,), 1.0, 3.0, None, out, stereo_guard=guard)
    abl = ablation.analyze(out, GRIDS, out)
    res["ablation_summary"] = abl["summary"]
    res["ablation_tests"] = abl["tests"]
    if suite in ("enantiomers", "diastereomers"):
        import pickle

        runs = pickle.load(open(out / "raw.pkl", "rb"))
        fid = {}
        for n, s in mols.items():
            xs = [m[1] for r in runs if r["task"]["molecule"] == n for m in r["minima"]]
            fid[n] = {"n_minima_checked": len(xs), "stereo_fidelity": stereo_fidelity(s, xs)}
        res["stereo_fidelity"] = fid
    if suite == "enantiomers":
        import pickle

        runs = pickle.load(open(out / "raw.pkl", "rb"))
        pairs = {}
        names = list(mols)
        for a, b in zip(names[::2], names[1::2], strict=True):
            rm = SymmetryRMSD(benchmark._system(mols[a]).heavy_mol())
            perm = atom_order_map(mols[a], mols[b])  # b's coordinates re-indexed into a's atom order

            def pooled(name: str, idx=None) -> list[tuple[float, np.ndarray]]:
                ms = [(m[0], m[1]) for r in runs if r["task"]["molecule"] == name for m in r["minima"]]
                e0 = min(e for e, _ in ms)
                return [(e, x if idx is None else x[idx]) for e, x in ms if e <= e0 + 3.0]

            pa, pb = pooled(a), pooled(b, perm)
            ga, gb = np.load(GRIDS / f"{a}_systematic.npz"), np.load(GRIDS / f"{b}_systematic.npz")
            pairs[f"{a} / {b}"] = {
                "grid_Emin_diff": float(ga["energies"][0] - gb["energies"][0]),
                "search_Emin_diff": float(min(e for e, _ in pa) - min(e for e, _ in pb)),
                "mirror_match_a_in_b": mirror_match(pa, pb, rm),
                "mirror_match_b_in_a": mirror_match(pb, pa, rm),
                "grid_mirror_match": mirror_match(
                    [(e, x) for e, x in zip(ga["energies"], ga["heavy"], strict=True) if e <= ga["energies"][0] + 3],
                    [
                        (e, x[perm])
                        for e, x in zip(gb["energies"], gb["heavy"], strict=True)
                        if e <= gb["energies"][0] + 3
                    ],
                    rm,
                ),
            }
            # handedness bias: per-method recall of R vs S (paired by seed)
            rows = [r for r in abl["rows"] if r["molecule"] in (a, b)]
            for meth in {r["label"] for r in rows}:
                ra = {r["seed"]: r["recall_final"] for r in rows if r["molecule"] == a and r["label"] == meth}
                rb = {r["seed"]: r["recall_final"] for r in rows if r["molecule"] == b and r["label"] == meth}
                d = np.array([ra[k] - rb[k] for k in sorted(ra)])
                pairs[f"{a} / {b}"][f"recall_diff[{meth}]"] = [float(d.mean()), float(d.std(ddof=1) / np.sqrt(len(d)))]
        res["enantiomer_pairs"] = pairs
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("suites", nargs="*", default=["conformational", "enantiomers", "diastereomers", "constitutional"])
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--steps", type=int, default=300)
    ap.add_argument("--guard", action="store_true", help="reject relaxed minima with inverted stereocentres")
    a = ap.parse_args()
    GRIDS.mkdir(parents=True, exist_ok=True)
    path = OUT / "isomers.json"
    allres = json.loads(path.read_text()) if path.exists() else {}
    for suite in a.suites:
        allres[suite + ("_guarded" if a.guard else "")] = run_suite(suite, a.seeds, a.steps, a.guard)
        path.write_text(json.dumps(allres, indent=1, default=float))
        print("done", suite, flush=True)


if __name__ == "__main__":
    main()
