"""Ground-truth check on alanine dipeptide: exhaustive grid vs basin-hopping proposals."""

import json
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from rdkit.Chem import rdMolTransforms

from perlin_mole.benchmark import Task, run_task
from perlin_mole.search import SymmetryRMSD
from perlin_mole.validate import (
    ALA_PHI,
    ALA_PSI,
    ALANINE_DIPEPTIDE,
    _sys,
    classify_ala,
    recall,
    systematic_search,
)

SMI = ALANINE_DIPEPTIDE
s = _sys(SMI)
print("DOFs:", s.n_dof, s.torsions)

grid = list(np.arange(-180, 180, 15.0))  # 24 x 24 = 576 starts on (phi, psi)
ref_e, ref_heavy, ref_full, n_starts = systematic_search(SMI, grid)
conf = s.mol.GetConformer()
rows = []
for e, xyz in zip(ref_e, ref_full, strict=True):
    for a, p in enumerate(xyz):
        conf.SetAtomPosition(a, p.tolist())
    phi = np.rad2deg(rdMolTransforms.GetDihedralRad(conf, *ALA_PHI))
    psi = np.rad2deg(rdMolTransforms.GetDihedralRad(conf, *ALA_PSI))
    rows.append({"dE": float(e - ref_e[0]), "phi": float(phi), "psi": float(psi), "basin": classify_ala(phi, psi)})
print(f"systematic search: {n_starts} starts -> {len(ref_e)} unique MMFF94 minima")
for r in rows:
    print(f"  dE={r['dE']:6.2f}  phi={r['phi']:7.1f}  psi={r['psi']:7.1f}  {r['basin']}")

rmsd = SymmetryRMSD(s.heavy_mol())
methods = ["etkdg", "reset", "gauss-iid", "gauss-corr", "perlin"]
tasks = [Task("ala1", SMI, m, seed, 100, 60.0, 2.0, 1.0) for m in methods for seed in range(10)]
with ProcessPoolExecutor() as pool:
    runs = list(pool.map(run_task, tasks))
summary = {}
for m in methods:
    rec = [
        recall(ref_e, ref_heavy, [(e, x) for e, x, _ in r["minima"]], rmsd, 5.0)[0]
        for r in runs
        if r["task"]["method"] == m
    ]
    summary[m] = {"recall_mean": float(np.mean(rec)), "recall_min": float(np.min(rec))}
    print(
        f"{m:11s} recall of minima within 5 kcal/mol after 100 relaxations: "
        f"mean {np.mean(rec):.2f}  min {np.min(rec):.2f}  (10 seeds)"
    )
# union of everything found by every method vs the grid
allfound = [(e, x) for r in runs for e, x, _ in r["minima"]]
u, _ = recall(ref_e, ref_heavy, allfound, rmsd, 100.0)
print(f"union of all methods recovers {u:.2f} of ALL grid minima")
json.dump(
    {"n_starts": n_starts, "minima": rows, "methods": summary, "union_recall_all": u},
    open("results/validation/alanine_dipeptide.json", "w"),
    indent=1,
)
