"""DFT validation on the GPU (GPU4PySCF): the highest-level judge in the study.

D1  conformational isomers vs gas-phase experiment: B3LYP-D3(BJ)/def2-TZVP geometry optimisation of the
    lowest anti/gauche (tt/tg/g+g+) conformer of each molecule (start: MMFF grid minima), Delta E.
D2  re-judging every discovered low-energy minimum: B3LYP-D3(BJ)/def2-SVP single points at the MMFF94
    geometries, compared with MMFF94 and with ANI-2x single points at the same geometries.
    Basis-set check: def2-TZVP single points for Ac-(Ala)3-NHMe.

Usage: uv run --extra dft --extra gpu python scripts/run_dft.py [d1] [d2] [d3]
Checkpoints after every molecule to results/dft/dft.json.
"""

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "dft"
OUT.mkdir(parents=True, exist_ok=True)
HARTREE = 627.509474
XC = "b3lyp"


def mf_for(symbols, xyz, basis, charge=0):
    from gpu4pyscf.dft import rks
    from pyscf import gto

    mol = gto.M(atom=list(zip(symbols, map(tuple, xyz), strict=True)), basis=basis, charge=charge, verbose=0)
    mf = rks.RKS(mol, xc=XC).density_fit()
    mf.disp = "d3bj"
    mf.conv_tol = 1e-9
    return mol, mf


def single_point(symbols, xyz, basis):
    _, mf = mf_for(symbols, xyz, basis)
    e = mf.kernel()
    return float(e) * HARTREE, bool(mf.converged)


def optimise(symbols, xyz, basis):
    from pyscf.geomopt.geometric_solver import optimize

    mol, mf = mf_for(symbols, xyz, basis)
    mol_eq = optimize(mf, maxsteps=200)
    e, conv = single_point([mol_eq.atom_symbol(i) for i in range(mol_eq.natm)], mol_eq.atom_coords(unit="Ang"), basis)
    return e, mol_eq.atom_coords(unit="Ang"), conv


CKPT = os.environ.get("DFT_JSON", "dft.json")  # separate files let D1 and D2 run in parallel processes


def load():
    p = OUT / CKPT
    return json.loads(p.read_text()) if p.exists() else {"d1": [], "d2": {}}


def save(res):
    (OUT / CKPT).write_text(json.dumps(res, indent=1))


def run_d1(res):
    from perlin_mole.conformational import classify
    from perlin_mole.isomers import CONFORMATIONAL

    exp = json.loads((ROOT / "data" / "experimental_conformer_energies.json").read_text())
    done = {(r["molecule"], r["quantity"]) for r in res["d1"]}
    for name, smi in CONFORMATIONAL.items():
        g = np.load(ROOT / "results" / "isomers" / "grids" / f"{name}_systematic.npz")
        symbols = [a.GetSymbol() for a in Chem.AddHs(Chem.MolFromSmiles(smi)).GetAtoms()]
        best = {}
        for xyz, e in zip(g["full"], g["energies"], strict=True):
            k = classify(smi, xyz)
            if k is not None and (k not in best or e < best[k][1]):
                best[k] = (xyz, e)
        t0 = time.time()
        dft = {}
        for k, (xyz, _) in best.items():
            e, xopt, conv = optimise(symbols, xyz, "def2-tzvp")
            k_after = classify(smi, xopt)
            dft[k] = {"E": e, "class_after_opt": k_after, "converged": conv}
        ref = exp[name]
        pairs = (
            [
                ("tg - tt", "tg", "tt", ref["value"], ref["err"]),
                ("g+g+ - tt", "gg+", "tt", ref["gg"]["value"], ref["gg"]["err"]),
            ]
            if name == "n-pentane"
            else [("gauche - anti", "g", "a", ref["value"], ref["err"])]
        )
        for q, k1, k0, v, err in pairs:
            if (name, q) in done or k1 not in dft or k0 not in dft:
                continue
            row = {
                "molecule": name,
                "quantity": q,
                "dft": dft[k1]["E"] - dft[k0]["E"],
                "experiment": v,
                "err": err,
                "mmff94": float(best[k1][1] - best[k0][1]),
                "classes_kept": dft[k1]["class_after_opt"] == k1 and dft[k0]["class_after_opt"] == k0,
                "seconds": round(time.time() - t0, 1),
            }
            res["d1"].append(row)
            print(row, flush=True)
        save(res)


def run_d2(res):
    from perlin_mole.ml_judge import energies as ani_energies
    from perlin_mole.ml_judge import supported

    for sdf in sorted((ROOT / "results" / "discoveries" / "sdf").glob("*_low_energy.sdf")):
        mol_name = sdf.name.replace("_low_energy.sdf", "")
        if mol_name in res["d2"]:
            continue
        ms = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=False)]
        symbols = [a.GetSymbol() for a in ms[0].GetAtoms()]
        z = np.array([a.GetAtomicNum() for a in ms[0].GetAtoms()])
        X = np.stack([m.GetConformer().GetPositions() for m in ms])
        e_mmff = np.array([float(m.GetProp("E_MMFF94_kcal")) for m in ms])
        e_ani = ani_energies(z, X) if supported(z) else np.full(len(ms), np.nan)
        t0 = time.time()
        e_dft, conv = [], []
        for x in X:
            e, c = single_point(symbols, x, "def2-svp")
            e_dft.append(e)
            conv.append(c)
        e_dft = np.array(e_dft)
        rel = lambda a: a - np.nanmin(a)  # noqa: E731
        d = {
            "n": len(ms),
            "seconds": round(time.time() - t0, 1),
            "n_converged": int(sum(conv)),
            "dE_mmff": rel(e_mmff).tolist(),
            "dE_ani": rel(e_ani).tolist(),
            "dE_dft": rel(e_dft).tolist(),
            "found_by": [m.GetProp("found_by") for m in ms],
            "spearman_dft_mmff": float(spearmanr(e_dft, e_mmff).statistic),
            "spearman_dft_ani": float(spearmanr(e_dft, e_ani).statistic) if supported(z) else None,
            "mae_mmff_vs_dft": float(np.mean(np.abs(rel(e_mmff) - rel(e_dft)))),
            "mae_ani_vs_dft": float(np.mean(np.abs(rel(e_ani) - rel(e_dft)))) if supported(z) else None,
            "dft_gm_rank_under_mmff": int(np.argsort(np.argsort(e_mmff))[int(np.argmin(e_dft))]),
            "dft_gm_rank_under_ani": int(np.argsort(np.argsort(e_ani))[int(np.argmin(e_dft))])
            if supported(z)
            else None,
            "mmff_gm_rank_under_dft": int(np.argsort(np.argsort(e_dft))[int(np.argmin(e_mmff))]),
            "ani_gm_rank_under_dft": int(np.argsort(np.argsort(e_dft))[int(np.nanargmin(e_ani))])
            if supported(z)
            else None,
        }
        if mol_name == "ala3":  # basis-set check
            e_tz = np.array([single_point(symbols, x, "def2-tzvp")[0] for x in X])
            d["dE_dft_tzvp"] = rel(e_tz).tolist()
            d["spearman_svp_tzvp"] = float(spearmanr(e_dft, e_tz).statistic)
            d["mae_svp_vs_tzvp"] = float(np.mean(np.abs(rel(e_dft) - rel(e_tz))))
        res["d2"][mol_name] = d
        print(mol_name, {k: v for k, v in d.items() if not isinstance(v, list)}, flush=True)
        save(res)


def run_d3(res):
    """def2-TZVP arbitration of the global minimum: candidates = MMFF GM, ANI GM, top-3 by def2-SVP DFT.

    def2-SVP single points over-stabilise compact, dispersion-bound conformers (basis-set superposition):
    for n-decane they put a 4-gauche conformer 0.39 kcal/mol below all-trans; def2-TZVP puts it 0.40 above.
    The global-minimum question is therefore settled at def2-TZVP.
    """
    d2 = {}
    for f in ("dft.json", "dft_d2.json"):
        p = OUT / f
        if p.exists():
            d2.update(json.loads(p.read_text()).get("d2", {}))
    res.setdefault("d3", {})
    for mol_name, d in d2.items():
        if mol_name in res["d3"]:
            continue
        ms = [
            m
            for m in Chem.SDMolSupplier(
                str(ROOT / "results" / "discoveries" / "sdf" / f"{mol_name}_low_energy.sdf"), removeHs=False
            )
        ]
        symbols = [a.GetSymbol() for a in ms[0].GetAtoms()]
        svp, mm, an = (np.array(d[k], dtype=float) for k in ("dE_dft", "dE_mmff", "dE_ani"))
        cand = list(
            dict.fromkeys(
                [int(np.argmin(mm))]
                + ([int(np.nanargmin(an))] if not np.isnan(an).all() else [])
                + [int(i) for i in np.argsort(svp)[:3]]
            )
        )
        e = {i: single_point(symbols, ms[i].GetConformer().GetPositions(), "def2-tzvp")[0] for i in cand}
        e0 = min(e.values())
        gm = min(e, key=e.get)
        res["d3"][mol_name] = {
            "candidates": {
                str(i): {
                    "dE_tzvp": e[i] - e0,
                    "dE_svp": float(svp[i]),
                    "dE_mmff": float(mm[i]),
                    "dE_ani": None if np.isnan(an[i]) else float(an[i]),
                }
                for i in cand
            },
            "tzvp_gm": gm,
            "mmff_gm_is_tzvp_gm": int(np.argmin(mm)) == gm,
            "ani_gm_is_tzvp_gm": (not np.isnan(an).all()) and int(np.nanargmin(an)) == gm,
            "svp_gm_is_tzvp_gm": int(np.argmin(svp)) == gm,
            "mmff_gm_dE_tzvp": e[int(np.argmin(mm))] - e0,
            "ani_gm_dE_tzvp": None if np.isnan(an).all() else e[int(np.nanargmin(an))] - e0,
            "svp_gm_dE_tzvp": e[int(np.argmin(svp))] - e0,
        }
        print(mol_name, {k: v for k, v in res["d3"][mol_name].items() if k != "candidates"}, flush=True)
        save(res)


if __name__ == "__main__":
    parts = sys.argv[1:] or ["d1", "d2"]
    res = load()
    if "d1" in parts:
        run_d1(res)
    if "d2" in parts:
        run_d2(res)
    if "d3" in parts:
        run_d3(res)
    print("done", flush=True)
