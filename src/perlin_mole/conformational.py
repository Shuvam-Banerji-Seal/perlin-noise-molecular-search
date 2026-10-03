"""Predicted (MMFF94, exhaustive grid) vs experimental gas-phase conformer energy differences."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from rdkit.Chem import rdMolTransforms

from perlin_mole.chem import TorsionSystem

ROOT = Path(__file__).resolve().parents[2]


def _label(phi_deg: float) -> str:
    a = abs(phi_deg)
    return "a" if a > 150 else ("g" if 30 < a < 100 else "x")


def classify(smiles: str, xyz: np.ndarray) -> str | None:
    """Conformer class from full-atom coordinates: 'a'/'g' (one torsion) or 'tt'/'tg'/'gg+'/'gg-' (two)."""
    s = TorsionSystem.from_smiles(smiles)
    conf = s.mol.GetConformer()
    for i, p in enumerate(xyz):
        conf.SetAtomPosition(i, p.tolist())
    phis = [np.rad2deg(rdMolTransforms.GetDihedralRad(conf, *t)) for t in s.torsions]
    labs = [_label(p) for p in phis]
    if "x" in labs:
        return None
    if len(phis) == 2:
        if labs == ["a", "a"]:
            return "tt"
        if "a" in labs:
            return "tg"
        return "gg+" if np.sign(phis[0]) == np.sign(phis[1]) else "gg-"
    return labs[0]


def conformer_energies(name: str, smiles: str, grid_dir: Path) -> dict[str, float]:
    """Lowest energy per conformer class (e.g. 'a', 'g', 'tt', 'tg', 'gg+', 'gg-') from the grid minima."""
    g = np.load(grid_dir / f"{name}_systematic.npz")
    s = TorsionSystem.from_smiles(smiles)
    conf = s.mol.GetConformer()
    best: dict[str, float] = {}
    for e, xyz in zip(g["energies"], g["full"], strict=True):
        for i, p in enumerate(xyz):
            conf.SetAtomPosition(i, p.tolist())
        phis = [np.rad2deg(rdMolTransforms.GetDihedralRad(conf, *t)) for t in s.torsions]
        labs = [_label(p) for p in phis]
        if "x" in labs:
            continue
        if len(phis) == 2:  # pentane: tt, tg, g+g+ (same sign), g+g- (opposite sign)
            if labs == ["a", "a"]:
                key = "tt"
            elif "a" in labs:
                key = "tg"
            else:
                key = "gg+" if np.sign(phis[0]) == np.sign(phis[1]) else "gg-"
        else:
            key = labs[0]
        best[key] = min(best.get(key, np.inf), float(e))
    return best


def compare(grid_dir: Path, exp_path: Path = ROOT / "data" / "experimental_conformer_energies.json") -> list[dict]:
    from perlin_mole.isomers import CONFORMATIONAL

    exp = json.loads(exp_path.read_text())
    rows = []
    for name, smi in CONFORMATIONAL.items():
        if name not in exp or not (grid_dir / f"{name}_systematic.npz").exists():
            continue
        b = conformer_energies(name, smi, grid_dir)
        ref = exp[name]
        if name == "n-pentane":
            pairs = [
                ("tg - tt", b.get("tg"), b.get("tt"), ref["value"], ref["err"]),
                ("g+g+ - tt", b.get("gg+"), b.get("tt"), ref["gg"]["value"], ref["gg"]["err"]),
            ]
        else:
            pairs = [("gauche - anti", b.get("g"), b.get("a"), ref["value"], ref["err"])]
        for q, e1, e0, v, err in pairs:
            if e1 is None or e0 is None:
                continue
            rows.append(
                {
                    "molecule": name,
                    "quantity": q,
                    "mmff94": e1 - e0,
                    "experiment": v,
                    "err": err,
                    "ref": ref["ref"],
                    "status": ref["status"],
                }
            )
    return rows
