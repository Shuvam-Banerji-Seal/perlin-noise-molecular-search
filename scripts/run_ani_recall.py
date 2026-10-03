"""Is the proposal comparison judge-agnostic? Re-score sigma=60 runs on the ANI-2x low-energy set.

For each molecule, the ANI-low set is the pooled minima within 1 kcal/mol of the ANI-2x minimum
(ANI energies from results/ml_judge). Each sigma=60 run is credited with the members it found
(identity as elsewhere: |dE_MMFF| < 0.1 kcal/mol and symmetry RMSD < 0.25 A). Then paired tests.
"""

import json
import pickle
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.stats import wilcoxon

from perlin_mole.search import SymmetryRMSD

ROOT = Path(__file__).resolve().parents[1]
judge = json.loads((ROOT / "results" / "ml_judge" / "ml_judge.json").read_text())["rerank"]
runs = pickle.load(open(ROOT / "results" / "sigma60" / "raw.pkl", "rb"))
out = {}
for mol, info in judge.items():
    ms = [
        m
        for m in Chem.SDMolSupplier(
            str(ROOT / "results" / "discoveries" / "sdf" / f"{mol}_low_energy.sdf"), removeHs=False
        )
    ]
    e_mmff = np.array([float(m.GetProp("E_MMFF94_kcal")) for m in ms])
    heavy = np.array([a.GetIdx() for a in ms[0].GetAtoms() if a.GetAtomicNum() > 1])
    xyz = np.stack([m.GetConformer().GetPositions()[heavy] for m in ms])
    d_ani = np.array(info["dE_ani"])
    target = np.flatnonzero(d_ani <= 1.0)
    rm = SymmetryRMSD(Chem.RemoveHs(ms[0]))
    rec = {}
    for r in runs:
        if r["task"]["molecule"] != mol:
            continue
        re_ = np.array([m[0] for m in r["minima"]])
        rx = np.stack([m[1] for m in r["minima"]])
        hit = 0
        for k in target:
            c = np.flatnonzero(np.abs(re_ - e_mmff[k]) < 0.1)
            hit += bool(c.size and (rm.many(xyz[k], rx[c]) < 0.25).any())
        rec[(r["label"], r["task"]["seed"])] = hit / len(target)
    out[mol] = {"n_ani_low": int(len(target)), "rec": rec}
    labs = sorted({k[0] for k in rec})
    print(
        mol,
        "ANI-low set size",
        len(target),
        {lab: round(np.mean([v for (l2, _), v in rec.items() if l2 == lab]), 3) for lab in labs},
    )

res = {}
for a, b in [
    ("perlin(L=4)", "gauss-corr(L=4)"),
    ("perlin(L=2)", "gauss-corr(L=2)"),
    ("perlin(L=8)", "gauss-corr(L=8)"),
    ("gauss-corr(L=4)", "gauss-iid"),
    ("etkdg", "gauss-iid"),
]:
    d = np.array([o["rec"][(a, s)] - o["rec"][(b, s)] for o in out.values() for s in range(10) if (a, s) in o["rec"]])
    p = float(wilcoxon(d).pvalue) if np.any(d != 0) else 1.0
    res[f"{a} - {b}"] = {"mean_diff": float(d.mean()), "wins": int((d > 0).sum()), "losses": int((d < 0).sum()), "p": p}
    print(f"{a:>16s} - {b:<16s} mean {d.mean():+.3f}  W/L {(d > 0).sum()}/{(d < 0).sum()}  p {p:.3g}")
json.dump(
    {"tests": res, "ani_low_sizes": {k: v["n_ani_low"] for k, v in out.items()}},
    open(ROOT / "results" / "ml_judge" / "ani_recall.json", "w"),
    indent=1,
)
