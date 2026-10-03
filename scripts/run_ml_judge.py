"""Re-judge discovered minima and conformational isomers with ANI-2x on the GPU.

Usage: TORCHANI_NO_WARN_EXTENSIONS=1 uv run --extra gpu --extra plot python scripts/run_ml_judge.py
Writes results/ml_judge/{ml_judge.json,summary.md} and figures/ml_judge.{pdf,png}.
"""

import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.stats import spearmanr

from perlin_mole.conformational import classify
from perlin_mole.isomers import CONFORMATIONAL
from perlin_mole.ml_judge import fire_relax, supported
from perlin_mole.search import SymmetryRMSD

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "ml_judge"
OUT.mkdir(parents=True, exist_ok=True)
EXP = json.loads((ROOT / "data" / "experimental_conformer_energies.json").read_text())


# ---------------------------------------------------------------- A. conformational isomers vs experiment
conf_rows = []
for name, smi in CONFORMATIONAL.items():
    g = np.load(ROOT / "results" / "isomers" / "grids" / f"{name}_systematic.npz")
    z = np.array([a.GetAtomicNum() for a in Chem.AddHs(Chem.MolFromSmiles(smi)).GetAtoms()])
    if not supported(z) or name not in EXP:
        print("skip (elements)", name)
        continue
    xr, er, conv = fire_relax(z, g["full"])
    best_mmff, best_ani = {}, {}
    for xm, em, xa, ea in zip(g["full"], g["energies"], xr, er, strict=True):
        k = classify(smi, xa)
        if k is None:
            continue
        best_ani[k] = min(best_ani.get(k, np.inf), float(ea))
        km = classify(smi, xm)
        if km is not None:
            best_mmff[km] = min(best_mmff.get(km, np.inf), float(em))
    ref = EXP[name]
    pairs = (
        [
            ("tg - tt", "tg", "tt", ref["value"], ref["err"]),
            ("g+g+ - tt", "gg+", "tt", ref["gg"]["value"], ref["gg"]["err"]),
        ]
        if name == "n-pentane"
        else [("gauche - anti", "g", "a", ref["value"], ref["err"])]
    )
    for q, k1, k0, v, err in pairs:
        conf_rows.append(
            {
                "molecule": name,
                "quantity": q,
                "experiment": v,
                "err": err,
                "mmff94": best_mmff[k1] - best_mmff[k0],
                "ani2x": best_ani[k1] - best_ani[k0],
                "ani_converged": int(conv.sum()),
                "n": int(len(conv)),
            }
        )
        print(conf_rows[-1])

# ---------------------------------------------------------------- B. re-ranking discovered minima
rerank = {}
for sdf in sorted((ROOT / "results" / "discoveries" / "sdf").glob("*_low_energy.sdf")):
    mol_name = sdf.name.replace("_low_energy.sdf", "")
    ms = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=False)]
    z = np.array([a.GetAtomicNum() for a in ms[0].GetAtoms()])
    if not supported(z):
        continue
    e_mmff = np.array([float(m.GetProp("E_MMFF94_kcal")) for m in ms])
    X = np.stack([m.GetConformer().GetPositions() for m in ms])
    xr, e_ani, conv = fire_relax(z, X, max_steps=3000)
    heavy = np.array([a.GetIdx() for a in ms[0].GetAtoms() if a.GetAtomicNum() > 1])
    rm = SymmetryRMSD(Chem.RemoveHs(ms[0]))
    shift = np.array([rm(X[i][heavy], xr[i][heavy]) for i in range(len(ms))])
    rho = spearmanr(e_mmff, e_ani).statistic
    order_ani = np.argsort(e_ani)
    rank_of_mmff_gm = int(np.where(order_ani == int(np.argmin(e_mmff)))[0][0])
    gm_ani = int(order_ani[0])
    rerank[mol_name] = {
        "n_minima": len(ms),
        "ani_converged": int(conv.sum()),
        "spearman_mmff_vs_ani": float(rho),
        "mmff_gm_rank_under_ani": rank_of_mmff_gm,
        "ani_gm_dE_mmff": float(e_mmff[gm_ani] - e_mmff.min()),
        "ani_gm_found_by": ms[gm_ani].GetProp("found_by"),
        "ani_gm_in_grid": ms[gm_ani].GetProp("in_grid_reference"),
        "median_geometry_shift_A": float(np.median(shift)),
        "max_geometry_shift_A": float(shift.max()),
        "dE_mmff": (e_mmff - e_mmff.min()).tolist(),
        "dE_ani": (e_ani - e_ani.min()).tolist(),
    }
    print(mol_name, {k: v for k, v in rerank[mol_name].items() if not isinstance(v, list)})

res = {"conformational": conf_rows, "rerank": rerank}
(OUT / "ml_judge.json").write_text(json.dumps(res, indent=1))

lines = [
    "# ANI-2x (GPU) as a second judge",
    "",
    "## Conformational isomers: MMFF94 vs ANI-2x vs gas-phase experiment",
    "",
    "| molecule | quantity | MMFF94 | ANI-2x | experiment |",
    "|---|---|---|---|---|",
]
for r in conf_rows:
    lines.append(
        f"| {r['molecule']} | {r['quantity']} | {r['mmff94']:+.2f} | {r['ani2x']:+.2f} | "
        f"{r['experiment']:+.2f} ± {r['err']:.2f} |"
    )
if conf_rows:
    mae_m = np.mean([abs(r["mmff94"] - r["experiment"]) for r in conf_rows])
    mae_a = np.mean([abs(r["ani2x"] - r["experiment"]) for r in conf_rows])
    lines += ["", f"MAE vs experiment: MMFF94 {mae_m:.2f}, ANI-2x {mae_a:.2f} kcal/mol (same molecules)"]
lines += [
    "",
    "## Re-ranking the discovered low-energy minima",
    "",
    "| molecule | minima | Spearman(MMFF, ANI) | rank of MMFF global min under ANI "
    "| MMFF dE of ANI global min | median geometry shift (A) |",
    "|---|---|---|---|---|---|",
]
for k, v in rerank.items():
    lines.append(
        f"| {k} | {v['n_minima']} | {v['spearman_mmff_vs_ani']:.2f} | {v['mmff_gm_rank_under_ani']} | "
        f"{v['ani_gm_dE_mmff']:.2f} | {v['median_geometry_shift_A']:.2f} |"
    )
(OUT / "summary.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))

# figure
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

fig, axes = plt.subplots(
    1, 1 + len(rerank), figsize=(7.4, 2.5), gridspec_kw={"width_ratios": [1.6] + [1] * len(rerank)}
)
ax = axes[0]
y = np.arange(len(conf_rows))
ax.errorbar(
    [r["experiment"] for r in conf_rows],
    y,
    xerr=[r["err"] for r in conf_rows],
    fmt="s",
    color="k",
    ms=4,
    capsize=2,
    label="experiment",
)
ax.scatter([r["mmff94"] for r in conf_rows], y, color="#c0392b", s=18, zorder=3, label="MMFF94")
ax.scatter([r["ani2x"] for r in conf_rows], y, color="#2471a3", marker="D", s=18, zorder=3, label="ANI-2x (GPU)")
ax.set_yticks(y, [f"{r['molecule'].replace('1,2-', '')}\n{r['quantity']}" for r in conf_rows], fontsize=5.5)
ax.invert_yaxis()
ax.axvline(0, color="0.8", lw=0.6)
ax.set_xlabel(r"$\Delta E$ (kcal/mol)", fontsize=7)
ax.legend(frameon=False, fontsize=5.5, loc="lower right")
ax.set_title("(a) predicted vs real", fontsize=7.5)
for ax, (k, v) in zip(axes[1:], rerank.items(), strict=True):
    ax.scatter(v["dE_mmff"], v["dE_ani"], s=6, color="#2471a3", alpha=0.7)
    lim = max(max(v["dE_mmff"]), max(v["dE_ani"])) * 1.05
    ax.plot([0, lim], [0, lim], color="0.6", lw=0.6)
    ax.set_title(f"{k}\n" + rf"$\rho$={v['spearman_mmff_vs_ani']:.2f}", fontsize=7)
    ax.set_xlabel(r"MMFF94 $\Delta E$", fontsize=6.5)
    ax.tick_params(labelsize=6)
axes[1].set_ylabel(r"ANI-2x $\Delta E$", fontsize=6.5)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(ROOT / "figures" / f"ml_judge.{ext}", dpi=200, bbox_inches="tight")
