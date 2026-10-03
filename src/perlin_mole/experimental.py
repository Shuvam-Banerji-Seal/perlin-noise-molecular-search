"""Predicted vs. experimentally observed conformations.

For molecules with an experimental 3-D structure (PDB ligand instances or
small-molecule crystal structures), every search method is scored on whether
and how fast it produces a local minimum close to the observed conformation:

* ``min_rmsd``       smallest heavy-atom RMSD of any minimum found to experiment
* ``min_rmsd_10``    the same, restricted to minima within 10 kcal/mol of the run's best
* ``rmsd_best``      RMSD of the run's lowest-energy minimum (the "prediction")
* ``success_1``/``success_05``  min_rmsd <= 1.0 / 0.5 Angstrom
* ``t_1``            first relaxation step yielding a minimum within 1.0 Angstrom (right-censored)
* ``match_rank``/``match_dE``   energy rank / dE of the best-matching minimum within the run

RMSD is symmetry-aware (minimised over graph automorphisms) and uses proper
rotations only, so mirror-image conformers of chiral molecules never match.
"""

from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers as ffh
from scipy.stats import wilcoxon

from perlin_mole.benchmark import Task, run_task
from perlin_mole.search import SymmetryRMSD

THRESHOLDS = (1.0, 0.5)
METHODS = ("etkdg", "reset", "gauss-iid", "gauss-corr", "perlin")


# --------------------------------------------------------------------------- structures
def canonical_with_coords(template: Chem.Mol, inst: Chem.Mol) -> tuple[str, Chem.Mol] | None:
    """Map an experimental instance onto the template; return (SMILES, heavy mol in SMILES order).

    ``inst`` must carry the template's bond orders (any atom order) and one conformer.
    Stereochemistry is read from the 3-D coordinates. Elements the template
    specifies must agree, and elements the template leaves unspecified are
    taken from the structure. The returned SMILES encodes the observed
    stereoisomer, so searches sample the same stereoisomer as the experiment.
    """
    m3d = Chem.Mol(inst)
    try:
        Chem.SanitizeMol(m3d)
        Chem.AssignStereochemistryFrom3D(m3d)
    except Exception:
        return None
    match = m3d.GetSubstructMatch(template, useChirality=False)
    if not match:
        return None

    # stereo elements the template does specify must agree with 3-D
    probe = Chem.Mol(m3d)
    for t_atom in template.GetAtoms():
        if t_atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            probe.GetAtomWithIdx(match[t_atom.GetIdx()]).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for t_bond in template.GetBonds():
        if t_bond.GetStereo() == Chem.BondStereo.STEREONONE:
            b = probe.GetBondBetweenAtoms(match[t_bond.GetBeginAtomIdx()], match[t_bond.GetEndAtomIdx()])
            if b is not None and b.GetBondType() == Chem.BondType.DOUBLE:
                b.SetStereo(Chem.BondStereo.STEREONONE)
    if Chem.MolToSmiles(probe) != Chem.MolToSmiles(template):
        return None

    smiles = Chem.MolToSmiles(m3d)
    mol = Chem.MolFromSmiles(smiles)
    order = m3d.GetSubstructMatch(mol, useChirality=True)
    if not order or len(order) != mol.GetNumAtoms():
        return None
    xyz = m3d.GetConformer().GetPositions()[list(order)]
    out = attach_coords(mol, xyz)
    check = Chem.Mol(out)
    Chem.AssignStereochemistryFrom3D(check)
    if Chem.MolToSmiles(check) != smiles:
        return None
    return smiles, out


def attach_coords(mol: Chem.Mol, xyz: np.ndarray) -> Chem.Mol:
    out = Chem.Mol(mol)
    out.RemoveAllConformers()
    conf = Chem.Conformer(out.GetNumAtoms())
    for i, p in enumerate(np.asarray(xyz, dtype=float)):
        conf.SetAtomPosition(i, p.tolist())
    conf.Set3D(True)
    out.AddConformer(conf, assignId=True)
    return out


def write_structure(path: Path, mol: Chem.Mol, name: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    m = Chem.Mol(mol)
    m.SetProp("_Name", name)
    m.SetProp("smiles", Chem.MolToSmiles(mol))
    w = Chem.SDWriter(str(path))
    w.write(m)
    w.close()
    return path


def load_structure(path: Path, smiles: str) -> np.ndarray:
    """Heavy-atom coordinates in ``Chem.MolFromSmiles(smiles)`` order (verified element by element)."""
    m = next(iter(Chem.SDMolSupplier(str(path), removeHs=True)))
    ref = Chem.MolFromSmiles(smiles)
    if [a.GetSymbol() for a in m.GetAtoms()] != [a.GetSymbol() for a in ref.GetAtoms()]:
        raise ValueError(f"{path}: atom order does not match SMILES")
    return m.GetConformer().GetPositions()


def relax_experimental(smiles: str, heavy_xyz: np.ndarray) -> dict:
    """MMFF94 view of the experimental geometry: H-only optimisation, then full relaxation."""
    mol = Chem.AddHs(attach_coords(Chem.MolFromSmiles(smiles), heavy_xyz), addCoords=True)
    props = ffh.MMFFGetMoleculeProperties(mol)
    ff = ffh.MMFFGetMoleculeForceField(mol, props)
    n_heavy = len(heavy_xyz)
    for i in range(n_heavy):
        ff.AddFixedPoint(i)
    ff.Minimize(maxIts=2000)
    e_frozen = ff.CalcEnergy()
    ff = ffh.MMFFGetMoleculeForceField(mol, props)
    status = ff.Minimize(maxIts=5000)
    relaxed = mol.GetConformer().GetPositions()[:n_heavy]
    rmsd = SymmetryRMSD(Chem.MolFromSmiles(smiles))(heavy_xyz, relaxed)
    return {
        "e_exp_frozen": float(e_frozen),
        "e_exp_relaxed": float(ff.CalcEnergy()),
        "rmsd_exp_to_own_minimum": float(rmsd),
        "relax_converged": status == 0,
        "relaxed_heavy": relaxed,
    }


# --------------------------------------------------------------------------- metrics
def run_metrics(run: dict, exp_xyz: np.ndarray, rmsd: SymmetryRMSD, n_steps: int) -> dict:
    e = np.array([m[0] for m in run["minima"]])
    x = np.stack([m[1] for m in run["minima"]])
    steps = np.array([m[2] for m in run["minima"]])
    r = rmsd.many(exp_xyz, x)
    best = int(np.argmin(e))
    k = int(np.argmin(r))
    within10 = e <= e.min() + 10.0
    hit = steps[r <= 1.0]
    rank = int((e < e[k]).sum()) + 1
    return {
        "min_rmsd": float(r.min()),
        "min_rmsd_10": float(r[within10].min()),
        "rmsd_best": float(r[best]),
        "success_1": bool(r.min() <= 1.0),
        "success_05": bool(r.min() <= 0.5),
        "t_1": int(hit.min()) if hit.size else None,
        "t_1_censored": int(hit.min()) if hit.size else n_steps + 1,
        "match_rank": rank,
        "match_dE": float(e[k] - e.min()),
        "best_energy": float(e.min()),
        "n_minima": int(len(e)),
        "best_match_xyz": x[k],
        "best_energy_xyz": x[best],
    }


_CTX: dict = {}


def _init(ctx: dict) -> None:
    _CTX.update(ctx)


def _work(task: Task) -> dict:
    run = run_task(task)
    exp = _CTX[task.molecule]
    m = run_metrics(run, exp["xyz"], SymmetryRMSD(Chem.MolFromSmiles(task.smiles)), task.n_steps)
    first_hit_curve = None
    return {"molecule": task.molecule, "method": task.method, "seed": task.seed, **m, "curve": first_hit_curve}


def _holm(ps: list[float]) -> list[float]:
    order = np.argsort(ps)
    m, adj, run = len(ps), [1.0] * len(ps), 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (m - rank) * ps[i]))
        adj[i] = run
    return adj


def _boot(d: np.ndarray, n: int = 4000) -> list[float]:
    rng = np.random.default_rng(0)
    means = d[rng.integers(0, len(d), size=(n, len(d)))].mean(axis=1)
    return [float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))]


def paired_tests(rows: list[dict], pairs: list[tuple[str, str]]) -> list[dict]:
    """Wilcoxon + bootstrap over (molecule, seed) pairs; Holm within each metric."""
    metrics = {"min_rmsd": "lower", "rmsd_best": "lower", "t_1_censored": "lower", "success_1": "higher"}
    out = []
    for metric, better in metrics.items():
        fam = []
        for a, b in pairs:
            va = {(r["molecule"], r["seed"]): float(r[metric]) for r in rows if r["method"] == a}
            vb = {(r["molecule"], r["seed"]): float(r[metric]) for r in rows if r["method"] == b}
            keys = sorted(set(va) & set(vb))
            d = np.array([va[k] - vb[k] for k in keys])
            nz = d[np.abs(d) > 1e-12]
            p = float(wilcoxon(d).pvalue) if nz.size >= 5 else 1.0
            sign = -1 if better == "lower" else 1
            fam.append(
                {
                    "a": a,
                    "b": b,
                    "metric": metric,
                    "n": len(d),
                    "mean_diff": float(d.mean()),
                    "ci95": _boot(d),
                    "a_better": int((sign * d > 1e-12).sum()),
                    "b_better": int((sign * d < -1e-12).sum()),
                    "p": p,
                }
            )
        for t, ph in zip(fam, _holm([t["p"] for t in fam]), strict=True):
            t["p_holm"] = ph
        out += fam
    return out


def evaluate(
    manifest: dict,
    data_dir: Path,
    methods: tuple[str, ...] = METHODS,
    seeds: int = 10,
    n_steps: int = 300,
    jobs: int | None = None,
) -> dict:
    ctx, exp_info = {}, {}
    for name, info in manifest.items():
        xyz = load_structure(data_dir / info["file"], info["smiles"])
        ctx[name] = {"xyz": xyz}
        rel = relax_experimental(info["smiles"], xyz)
        exp_info[name] = {**info, **{k: v for k, v in rel.items() if k != "relaxed_heavy"}}
    tasks = [
        Task(name, manifest[name]["smiles"], m, seed, n_steps, 60.0, 4.0, 1.0)
        for name in manifest
        for m in methods
        for seed in range(seeds)
    ]
    # heaviest molecules first keeps the pool busy
    tasks.sort(key=lambda t: -manifest[t.molecule]["n_heavy"])
    jobs = jobs or os.cpu_count() or 1
    with ProcessPoolExecutor(max_workers=jobs, initializer=_init, initargs=(ctx,)) as pool:
        rows = list(pool.map(_work, tasks, chunksize=1))

    for name in manifest:
        best = min(r["best_energy"] for r in rows if r["molecule"] == name)
        exp_info[name]["e_best_found"] = best
        exp_info[name]["e_exp_relaxed_minus_best"] = exp_info[name]["e_exp_relaxed"] - best
        exp_info[name]["e_exp_frozen_minus_best"] = exp_info[name]["e_exp_frozen"] - best

    summary = {}
    for m in methods:
        rr = [r for r in rows if r["method"] == m]
        tt = np.array([r["t_1_censored"] for r in rr])
        summary[m] = {
            "success_1": float(np.mean([r["success_1"] for r in rr])),
            "success_05": float(np.mean([r["success_05"] for r in rr])),
            "median_min_rmsd": float(np.median([r["min_rmsd"] for r in rr])),
            "mean_min_rmsd": float(np.mean([r["min_rmsd"] for r in rr])),
            "median_min_rmsd_10": float(np.median([r["min_rmsd_10"] for r in rr])),
            "median_rmsd_best": float(np.median([r["rmsd_best"] for r in rr])),
            "rmst_1": float(np.minimum(tt, n_steps).mean()),
            "survival_1": [float((tt > s).mean()) for s in range(n_steps + 1)],
            "per_molecule": {
                name: {
                    "success_1": float(np.mean([r["success_1"] for r in rr if r["molecule"] == name])),
                    "median_min_rmsd": float(np.median([r["min_rmsd"] for r in rr if r["molecule"] == name])),
                    "median_rmsd_best": float(np.median([r["rmsd_best"] for r in rr if r["molecule"] == name])),
                }
                for name in manifest
            },
        }
    pairs = [("perlin", "gauss-corr"), ("perlin", "gauss-iid"), ("perlin", "etkdg"), ("perlin", "reset")]
    tests = paired_tests(rows, pairs)
    return {
        "config": {"methods": list(methods), "seeds": seeds, "n_steps": n_steps, "sigma_deg": 60.0, "corr_length": 4.0},
        "structures": exp_info,
        "summary": summary,
        "tests": tests,
        "rows": rows,
    }


def strip_arrays(report: dict) -> dict:
    """JSON-safe copy (drops per-run coordinate arrays)."""
    rows = [{k: v for k, v in r.items() if not isinstance(v, np.ndarray) and k != "curve"} for r in report["rows"]]
    return {**report, "rows": rows}


def format_markdown(rep: dict) -> str:
    cfg = rep["config"]
    lines = [
        "# Predicted vs. experimental conformations",
        "",
        f"{len(rep['structures'])} molecules with experimental 3-D structures; {cfg['seeds']} seeds x "
        f"{cfg['n_steps']} MMFF94 relaxations per method; "
        f"sigma = {cfg['sigma_deg']:g} deg, L = {cfg['corr_length']:g}.",
        "RMSD: heavy atoms, symmetry-aware, proper rotations only.",
        "",
        "## Structures",
        "",
        "| molecule | source | ID | resolution / R | torsions | heavy atoms | E(exp, H-opt) − E* | "
        "E(exp relaxed) − E* | "
        "RMSD exp→own minimum |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for name, s in rep["structures"].items():
        q = f"{s['resolution_A']:.2f} Å" if s.get("resolution_A") else f"R = {s.get('r_factor', float('nan')):.3f}"
        lines.append(
            f"| {name} | {s['source']} | {s['id']} | {q} | {s['n_torsions']} | {s['n_heavy']} | "
            f"{s['e_exp_frozen_minus_best']:.1f} | {s['e_exp_relaxed_minus_best']:.1f} | "
            f"{s['rmsd_exp_to_own_minimum']:.2f} |"
        )
    lines += [
        "",
        "## Methods (pooled over molecules and seeds)",
        "",
        "| method | success ≤1.0 Å | success ≤0.5 Å | median min-RMSD | median min-RMSD (≤10 kcal) | "
        "median RMSD of lowest-E minimum | restricted mean steps to ≤1 Å |",
        "|---|---|---|---|---|---|---|",
    ]
    for m, s in rep["summary"].items():
        lines.append(
            f"| {m} | {s['success_1']:.2f} | {s['success_05']:.2f} | {s['median_min_rmsd']:.2f} | "
            f"{s['median_min_rmsd_10']:.2f} | {s['median_rmsd_best']:.2f} | {s['rmst_1']:.0f} |"
        )
    lines += [
        "",
        "## Paired tests over (molecule, seed)",
        "",
        "| metric | comparison | mean diff [95% CI] | a better / b better | p (Holm) |",
        "|---|---|---|---|---|",
    ]
    for t in rep["tests"]:
        lines.append(
            f"| {t['metric']} | {t['a']} − {t['b']} | "
            f"{t['mean_diff']:+.3f} [{t['ci95'][0]:+.3f}, {t['ci95'][1]:+.3f}] | "
            f"{t['a_better']} / {t['b_better']} | {t['p']:.3g} ({t['p_holm']:.2g}) |"
        )
    lines += ["", "## Per molecule: success ≤1.0 Å / median min-RMSD", ""]
    methods = list(rep["summary"])
    lines += ["| molecule | " + " | ".join(methods) + " |", "|---" * (len(methods) + 1) + "|"]
    for name in rep["structures"]:
        cells = [
            f"{rep['summary'][m]['per_molecule'][name]['success_1']:.1f} / "
            f"{rep['summary'][m]['per_molecule'][name]['median_min_rmsd']:.2f}"
            for m in methods
        ]
        lines.append(f"| {name} | " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def save(report: dict, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "experimental.json").write_text(json.dumps(strip_arrays(report), indent=1))
    (out_dir / "summary.md").write_text(format_markdown(report))


def best_alignment(exp_xyz: np.ndarray, pred_xyz: np.ndarray, perms: np.ndarray) -> tuple[np.ndarray, float]:
    """Superpose ``pred`` onto ``exp`` (best automorphism + proper rotation); return (aligned pred, RMSD)."""
    best = (None, np.inf)
    p0 = exp_xyz - exp_xyz.mean(0)
    for perm in perms:
        q = pred_xyz[perm]
        q0 = q - q.mean(0)
        u, _, vt = np.linalg.svd(q0.T @ p0)
        d = np.sign(np.linalg.det(u @ vt))
        rot = u @ np.diag([1.0, 1.0, d]) @ vt
        aligned = q0 @ rot + exp_xyz.mean(0)
        r = float(np.sqrt(((aligned - exp_xyz) ** 2).sum(1).mean()))
        if r < best[1]:
            best = (aligned, r)
    return best
