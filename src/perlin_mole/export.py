"""Pool all runs (all sigmas) + the systematic grid reference into one minima set per molecule.

Produces what the per-sweep summaries cannot:

* coverage on a **common** reference, so different sigmas are comparable;
* **recall against ground truth** (the exhaustive grid's low-energy minima);
* **discoveries**: low-energy minima found by basin hopping but *not* by the
  grid, and grid minima that every method missed;
* SDF files of all unique low-energy minima. Hydrogens are rebuilt by MMFF
  with heavy atoms frozen, and the recomputed energy is checked against the
  recorded one.
"""

from __future__ import annotations

import json
import os
import pickle
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers as ffh

from perlin_mole.benchmark import _system
from perlin_mole.search import SymmetryRMSD, deduplicate
from perlin_mole.validate import keeps_fixed_features

GRID = "grid"


def _label(r: dict) -> str:
    return f"{r['label']} @{r['task']['sigma_deg']:g}"


def _mean_sd(v: list[float]) -> list[float]:
    return [float(np.mean(v)), float(np.std(v, ddof=1)) if len(v) > 1 else 0.0]


def _with_h(smiles: str, heavy_xyz: np.ndarray) -> tuple[Chem.Mol, float]:
    """Rebuild hydrogens on frozen heavy atoms and relax them with MMFF94."""
    mol = Chem.MolFromSmiles(smiles)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, p in enumerate(heavy_xyz):
        conf.SetAtomPosition(i, p.tolist())
    mol.AddConformer(conf, assignId=True)
    mol = Chem.AddHs(mol, addCoords=True)
    props = ffh.MMFFGetMoleculeProperties(mol)
    ff = ffh.MMFFGetMoleculeForceField(mol, props)
    for i in range(len(heavy_xyz)):
        ff.AddFixedPoint(i)
    ff.Minimize(maxIts=2000)
    return mol, ff.CalcEnergy()


def analyze(
    raw_dirs: list[str | Path],
    out_dir: str | Path = "results/discoveries",
    systematic_dir: str | Path | None = "results/validation",
    window: float = 3.0,
    jobs: int | None = None,
) -> dict:
    jobs = jobs or os.cpu_count() or 1
    out = Path(out_dir)
    (out / "sdf").mkdir(parents=True, exist_ok=True)
    runs: list[dict] = []
    for d in raw_dirs:
        with open(Path(d) / "raw.pkl", "rb") as fh:
            runs += pickle.load(fh)

    report: dict = {"window": window, "sources": [str(d) for d in raw_dirs], "molecules": {}}
    for mol in dict.fromkeys(r["task"]["molecule"] for r in runs):
        mruns = [r for r in runs if r["task"]["molecule"] == mol]
        smiles = mruns[0]["task"]["smiles"]
        e = [m[0] for r in mruns for m in r["minima"]]
        x = [m[1] for r in mruns for m in r["minima"]]
        src = [k for k, r in enumerate(mruns) for _ in r["minima"]]
        grid_file = Path(systematic_dir) / f"{mol}_systematic.npz" if systematic_dir else None
        has_grid = bool(grid_file and grid_file.exists())
        if has_grid:
            g = np.load(grid_file)
            ok = [i for i in range(len(g["energies"])) if keeps_fixed_features(smiles, g["heavy"][i])]
            e += g["energies"][ok].tolist()
            x += list(g["heavy"][ok])
            src += [-1] * len(ok)
        e, x, src = np.asarray(e), np.stack(x), np.asarray(src)
        order = np.argsort(e, kind="stable")
        e, x, src = e[order], x[order], src[order]
        rep = deduplicate(e, x, SymmetryRMSD(_system(smiles).heavy_mol()), jobs=jobs)
        unique = np.flatnonzero(rep == np.arange(len(rep)))
        e_star = float(e[0])
        low = [int(i) for i in unique if e[i] <= e_star + window]

        found_by: dict[int, set[str]] = defaultdict(set)
        runs_finding: dict[int, set[int]] = defaultdict(set)
        for i, k in enumerate(src):
            r = int(rep[i])
            found_by[r].add(GRID if k < 0 else _label(mruns[k]))
            if k >= 0:
                runs_finding[r].add(int(k))

        # coverage on the common reference, per (method, sigma)
        cov: dict[str, list[float]] = defaultdict(list)
        grid_low = [i for i in low if GRID in found_by[i]]
        rec: dict[str, list[float]] = defaultdict(list)
        for k, r in enumerate(mruns):
            hits = {i for i in low if k in runs_finding[i]}
            cov[_label(r)].append(len(hits) / len(low))
            if grid_low:
                rec[_label(r)].append(len(hits & set(grid_low)) / len(grid_low))

        beyond = [i for i in low if GRID not in found_by[i]] if has_grid else []
        missed = [i for i in grid_low if len(found_by[i]) == 1]

        # SDF of every unique low-energy minimum, H rebuilt and energy re-checked
        writer = Chem.SDWriter(str(out / "sdf" / f"{mol}_low_energy.sdf"))
        e_check = []
        for rank, i in enumerate(low):
            m, e_h = _with_h(smiles, x[i])
            e_check.append(abs(e_h - e[i]))
            m.SetProp("_Name", f"{mol}_min{rank:03d}")
            m.SetProp("E_MMFF94_kcal", f"{e[i]:.4f}")
            m.SetProp("dE_kcal", f"{e[i] - e_star:.4f}")
            m.SetProp("E_recomputed_with_H", f"{e_h:.4f}")
            m.SetProp("in_grid_reference", "yes" if GRID in found_by[i] else ("no" if has_grid else "n/a"))
            m.SetProp("found_by", ";".join(sorted(found_by[i])))
            m.SetProp("n_runs_found", str(len(runs_finding[i])))
            writer.write(m)
        writer.close()

        report["molecules"][mol] = {
            "smiles": smiles,
            "e_star": e_star,
            "n_unique_total": int(len(unique)),
            "n_low": len(low),
            "has_grid": has_grid,
            "n_grid_low": len(grid_low),
            "n_beyond_grid_low": len(beyond),
            "beyond_grid_dE": [round(float(e[i] - e_star), 3) for i in beyond],
            "n_grid_low_missed_by_all_methods": len(missed),
            "global_min_found_by": sorted(found_by[low[0]]),
            "global_min_in_grid": GRID in found_by[low[0]],
            "coverage_common": {k: _mean_sd(v) for k, v in cov.items()},
            "recall_vs_grid": {k: _mean_sd(v) for k, v in rec.items()},
            "h_rebuild_max_abs_dE": float(max(e_check)) if e_check else 0.0,
        }

    (out / "discoveries.json").write_text(json.dumps(report, indent=1))
    (out / "discoveries.md").write_text(format_discoveries(report))
    return report


def format_discoveries(report: dict) -> str:
    lines = ["# Discoveries: pooled minima across all runs and the grid reference", ""]
    for mol, r in report["molecules"].items():
        more = " …" if r["n_beyond_grid_low"] > 12 else ""
        lines += [
            f"## {mol}",
            "",
            f"- unique minima (all sources): {r['n_unique_total']}; within {report['window']} kcal/mol "
            f"of E* = {r['e_star']:.3f}: **{r['n_low']}**",
        ]
        if r["has_grid"]:
            lines += [
                f"- grid reference low-energy minima: {r['n_grid_low']}; "
                f"missed by every method: **{r['n_grid_low_missed_by_all_methods']}**",
                f"- low-energy minima found by basin hopping but **not** by the grid: "
                f"**{r['n_beyond_grid_low']}** (ΔE: {r['beyond_grid_dE'][:12]}{more})",
                f"- global minimum in grid reference: {r['global_min_in_grid']}",
            ]
        lines += [
            f"- global minimum found by: {', '.join(r['global_min_found_by'])}",
            f"- H-rebuild energy check: max |ΔE| = {r['h_rebuild_max_abs_dE']:.3f} kcal/mol",
            "",
            "| method @ sigma | coverage (common ref) | recall vs grid |",
            "|---|---|---|",
        ]
        for k in sorted(r["coverage_common"], key=lambda k: -r["coverage_common"][k][0]):
            c = r["coverage_common"][k]
            g = r["recall_vs_grid"].get(k)
            gs = f"{g[0]:.3f} ± {g[1]:.3f}" if g else "-"
            lines.append(f"| {k} | {c[0]:.3f} ± {c[1]:.3f} | {gs} |")
        lines.append("")
    return "\n".join(lines)


def main() -> None:  # pragma: no cover - thin CLI wrapper
    import argparse

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("dirs", nargs="+", help="benchmark output dirs containing raw.pkl")
    p.add_argument("--out", default="results/discoveries")
    p.add_argument("--systematic", default="results/validation")
    p.add_argument("--window", type=float, default=3.0)
    p.add_argument("--jobs", type=int, default=None)
    a = p.parse_args()
    print(format_discoveries(analyze(a.dirs, a.out, a.systematic, a.window, a.jobs)))


if __name__ == "__main__":
    main()
