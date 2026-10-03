"""Predicted vs. experimental conformations: run the evaluation and draw the figures.

uv run --extra plot python scripts/run_experimental.py                 # full run (10 seeds x 300 steps)
uv run --extra plot python scripts/run_experimental.py --figures-only  # redraw from saved results
"""

from __future__ import annotations

import argparse
import json
import pickle
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem

from perlin_mole.experimental import METHODS, best_alignment, evaluate, load_structure, save
from perlin_mole.search import SymmetryRMSD

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental"
FIG = ROOT / "figures"
DEFAULT = [
    # PDB ligands (protein-bound, X-ray)
    "imatinib",
    "dasatinib",
    "erlotinib",
    "gefitinib",
    "sildenafil",
    "rosiglitazone",
    "oseltamivir-carboxylate",
    "ibuprofen",
    "indomethacin",
    "biotin",
    "methotrexate",
    # COD small-molecule crystal structures
    "tolbutamide",
    "aspirin",
    "fenofibrate",
    "ketoprofen",
    "ibuprofen-crystal",
]
COLORS = {"etkdg": "#1c1c1c", "reset": "#b9770e", "gauss-iid": "#7f8c8d", "gauss-corr": "#2471a3", "perlin": "#c0392b"}
CPK = {"C": "#404040", "O": "#e03b2c", "N": "#3050f8", "S": "#d4a017", "F": "#2ecc71", "Cl": "#27ae60", "P": "#ff8000"}

plt.rcParams.update(
    {
        "font.size": 8,
        "axes.titlesize": 8.5,
        "axes.labelsize": 8,
        "legend.fontsize": 7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.bbox": "tight",
        "savefig.dpi": 600,  # raster parts (imshow) of the vector PDFs: ACS asks >= 300 dpi
        "pdf.fonttype": 42,
    }
)


_PANEL = re.compile(r"^\(([a-z])\)\s*(.*)$", re.S)


def _bold_panel_labels(fig) -> None:
    """Scientific Reports style: '(a) title' -> bold 'a' followed by the title."""
    for ax in fig.axes:
        m = _PANEL.match(ax.get_title())
        if m:
            ax.title.set_text(rf"$\mathbf{{{m.group(1)}}}$  {m.group(2)}")  # keeps size and pad


def _save(fig, name: str) -> None:
    _bold_panel_labels(fig)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png", dpi=200)
    plt.close(fig)
    print("wrote", name)


def fig_rmsd(rep: dict) -> None:
    rows, methods = rep["rows"], rep["config"]["methods"]
    mols = list(rep["structures"])
    fig = plt.figure(figsize=(7.2, 3.0))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.2, 0.9, 0.95, 1.5], wspace=0.45)  # gs[2]: room for labels
    ax = fig.add_subplot(gs[0])
    data = [[r["min_rmsd"] for r in rows if r["method"] == m] for m in methods]
    bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False)
    for patch, m in zip(bp["boxes"], methods, strict=True):
        patch.set_facecolor(COLORS[m])
        patch.set_alpha(0.35)
    rng = np.random.default_rng(0)
    for i, (d, m) in enumerate(zip(data, methods, strict=True)):
        ax.scatter(i + 1 + rng.uniform(-0.18, 0.18, len(d)), d, s=3, color=COLORS[m], alpha=0.5, lw=0)
    ax.axhline(1.0, color="0.5", lw=0.6, ls="--")
    ax.set_xticks(range(1, len(methods) + 1), methods, rotation=40, ha="right", rotation_mode="anchor")
    ax.set_ylabel("min RMSD to experiment (Å)")
    ax.set_title("(a) closest conformer found")
    ax = fig.add_subplot(gs[1])
    x = np.arange(len(methods))
    s1 = [rep["summary"][m]["success_1"] for m in methods]
    s05 = [rep["summary"][m]["success_05"] for m in methods]
    ax.bar(x - 0.2, s1, 0.4, color=[COLORS[m] for m in methods], label="≤ 1.0 Å")
    ax.bar(x + 0.2, s05, 0.4, color=[COLORS[m] for m in methods], alpha=0.45, label="≤ 0.5 Å")
    ax.set_xticks(x, methods, rotation=40, ha="right", rotation_mode="anchor")
    ax.set_ylim(0, 1.25)
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_ylabel("fraction of runs")
    ax.set_title("(b) success rate")
    ax.legend(frameon=False, loc="upper center", ncol=2, fontsize=6, handlelength=1.2, columnspacing=0.8)
    ax = fig.add_subplot(gs[3])
    mat = np.array([[rep["summary"][m]["per_molecule"][n]["median_min_rmsd"] for m in methods] for n in mols])
    im = ax.imshow(mat, cmap="viridis_r", vmin=0, vmax=max(2.0, float(np.nanmax(mat))), aspect="auto")
    ax.set_xticks(range(len(methods)), methods, rotation=40, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(mols)), [f"{n} ({rep['structures'][n]['source']})" for n in mols], fontsize=6.5)
    for i in range(len(mols)):
        for j in range(len(methods)):
            ax.text(
                j, i, f"{mat[i, j]:.1f}", ha="center", va="center", fontsize=5.5, color="w" if mat[i, j] > 1.2 else "k"
            )
    ax.set_title("(c) median min RMSD per molecule (Å)")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    _save(fig, "experimental_rmsd")


def _draw(
    ax,
    mol: Chem.Mol,
    xyz: np.ndarray,
    frame: np.ndarray,
    centre: np.ndarray,
    color: str | None,
    lw: float,
    alpha: float,
    z: float,
) -> None:
    p = (xyz - centre) @ frame.T
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if color is None:
            mid = (p[i] + p[j]) / 2
            for a, q in ((i, p[i]), (j, p[j])):
                ax.plot(
                    [q[0], mid[0]],
                    [q[1], mid[1]],
                    color=CPK.get(mol.GetAtomWithIdx(a).GetSymbol(), "#888"),
                    lw=lw,
                    alpha=alpha,
                    solid_capstyle="round",
                    zorder=z,
                )
        else:
            ax.plot(
                [p[i, 0], p[j, 0]],
                [p[i, 1], p[j, 1]],
                color=color,
                lw=lw,
                alpha=alpha,
                solid_capstyle="round",
                zorder=z,
            )


def fig_overlays(rep: dict, overlays: dict, names: list[str]) -> None:
    cols = 3
    rows = int(np.ceil(len(names) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(7.2, 2.35 * rows))
    for ax, name in zip(axes.ravel(), names, strict=False):
        info = rep["structures"][name]
        mol = Chem.MolFromSmiles(info["smiles"])
        exp = overlays[name]["exp"]
        centre = exp.mean(0)
        _, _, vt = np.linalg.svd(exp - centre, full_matrices=False)
        _draw(ax, mol, exp, vt, centre, "#9a9a9a", 5.0, 0.55, 1)
        pred, r = overlays[name]["perlin"]
        _draw(ax, mol, pred, vt, centre, None, 1.6, 1.0, 3)
        rg = overlays[name]["gauss-corr"][1]
        re = overlays[name]["etkdg"][1]
        q = f"{info['resolution_A']:.2f} Å" if info.get("resolution_A") else f"R={info.get('r_factor', 0):.3f}"
        ax.set_title(
            f"{name} — {info['source']} {info['id']} ({q})\nRMSD perlin {r:.2f} · gauss-corr {rg:.2f}"
            f" · etkdg {re:.2f} Å",
            fontsize=6.8,
        )
        ax.set_aspect("equal")
        ax.axis("off")
    for ax in axes.ravel()[len(names) :]:
        ax.axis("off")
    fig.text(
        0.5,
        0.0,
        "grey: experimental structure; coloured: closest Perlin-search minimum after superposition",
        ha="center",
        fontsize=7,
    )
    fig.tight_layout()
    _save(fig, "experimental_overlays")


def fig_survival(rep: dict) -> None:
    fig, ax = plt.subplots(figsize=(3.6, 2.5))
    for m in rep["config"]["methods"]:
        s = np.array(rep["summary"][m]["survival_1"])
        ax.step(np.arange(len(s)), 1 - s, where="post", color=COLORS[m], lw=1.6 if m == "perlin" else 1.1, label=m)
    ax.set_xlabel("MMFF relaxations")
    ax.set_ylabel("runs with a ≤ 1 Å conformer")
    ax.set_ylim(0, 1)
    ax.set_title("time to reproduce the experimental conformation")
    ax.legend(frameon=False, fontsize=6.5, loc="lower right")
    _save(fig, "experimental_survival")


def fig_main(rep: dict, overlays: dict, names: list[str]) -> None:
    """Scientific Reports Fig. 4: predicted-vs-experimental overlays (a) and the quantitative comparison (b-d)."""
    methods = rep["config"]["methods"]
    fig = plt.figure(figsize=(7.2, 4.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.9, 1], hspace=0.3, wspace=0.42)
    for k, name in enumerate(names[:3]):
        ax = fig.add_subplot(gs[0, k])
        info = rep["structures"][name]
        mol = Chem.MolFromSmiles(info["smiles"])
        exp = overlays[name]["exp"]
        centre = exp.mean(0)
        _, _, vt = np.linalg.svd(exp - centre, full_matrices=False)
        _draw(ax, mol, exp, vt, centre, "#9a9a9a", 5.0, 0.55, 1)
        pred, r = overlays[name]["perlin"]
        _draw(ax, mol, pred, vt, centre, None, 1.5, 1.0, 3)
        ax.set_title(
            f"{name}, {info['source']} {info['id']}\nPerlin {r:.2f} Å, ETKDG {overlays[name]['etkdg'][1]:.2f} Å",
            fontsize=6.5,
        )
        ax.set_aspect("equal")
        ax.axis("off")
        if k == 0:
            ax.text(-0.08, 1.12, "a", transform=ax.transAxes, fontweight="bold", fontsize=10, va="bottom")
    rows = rep["rows"]
    ax = fig.add_subplot(gs[1, 0])
    data = [[r["min_rmsd"] for r in rows if r["method"] == m] for m in methods]
    bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False, whis=(0, 100))
    for patch, m in zip(bp["boxes"], methods, strict=True):
        patch.set_facecolor(COLORS[m])
        patch.set_alpha(0.35)
    ax.axhline(1.0, color="0.5", lw=0.6, ls="--")
    ax.set_xticks(range(1, len(methods) + 1), methods, rotation=40, ha="right", rotation_mode="anchor")
    ax.set_ylabel("min RMSD to experiment (Å)")
    ax.set_title("(b) closest conformer per run")
    ax = fig.add_subplot(gs[1, 1])
    x = np.arange(len(methods))
    ax.bar(
        x - 0.2,
        [rep["summary"][m]["success_1"] for m in methods],
        0.4,
        color=[COLORS[m] for m in methods],
        label="≤ 1.0 Å",
    )
    ax.bar(
        x + 0.2,
        [rep["summary"][m]["success_05"] for m in methods],
        0.4,
        color=[COLORS[m] for m in methods],
        alpha=0.45,
        label="≤ 0.5 Å",
    )
    ax.set_xticks(x, methods, rotation=40, ha="right", rotation_mode="anchor")
    ax.set_ylim(0, 1.25)
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_ylabel("fraction of runs")
    ax.set_title("(c) success rate")
    ax.legend(frameon=False, loc="upper center", ncol=2, fontsize=6, handlelength=1.2, columnspacing=0.8)
    ax = fig.add_subplot(gs[1, 2])
    for m in methods:
        sv = np.array(rep["summary"][m]["survival_1"])
        ax.step(np.arange(len(sv)), 1 - sv, where="post", color=COLORS[m], lw=1.6 if m == "perlin" else 1.0, label=m)
    ax.set_xlabel("MMFF94 relaxations")
    ax.set_ylabel("runs with a ≤ 1 Å conformer")
    ax.set_ylim(0, 1)
    ax.set_title("(d) time to reproduce")
    ax.legend(frameon=False, fontsize=5.5, loc="lower right")
    _save(fig, "experimental_main")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--molecules", default=",".join(DEFAULT))
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--steps", type=int, default=300)
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--out", default=str(ROOT / "results" / "experimental"))
    ap.add_argument("--figures-only", action="store_true")
    a = ap.parse_args()
    out = Path(a.out)
    manifest = json.loads((DATA / "manifest.json").read_text())
    names = [n for n in a.molecules.split(",") if n in manifest]

    if not a.figures_only:
        rep = evaluate({n: manifest[n] for n in names}, DATA, METHODS, a.seeds, a.steps, a.jobs)
        overlays = {}
        for n in names:
            exp = load_structure(DATA / manifest[n]["file"], manifest[n]["smiles"])
            perms = SymmetryRMSD(Chem.MolFromSmiles(manifest[n]["smiles"])).perms
            overlays[n] = {"exp": exp}
            for m in METHODS:
                rr = [r for r in rep["rows"] if r["molecule"] == n and r["method"] == m]
                best = min(rr, key=lambda r: r["min_rmsd"])
                overlays[n][m] = best_alignment(exp, best["best_match_xyz"], perms)
        save(rep, out)
        with open(out / "overlays.pkl", "wb") as fh:
            pickle.dump(overlays, fh)
    rep = json.loads((out / "experimental.json").read_text())
    with open(out / "overlays.pkl", "rb") as fh:
        overlays = pickle.load(fh)
    fig_rmsd(rep)
    pick = [
        n
        for n in ("imatinib", "oseltamivir-carboxylate", "methotrexate", "sildenafil", "fenofibrate", "tolbutamide")
        if n in overlays
    ]
    fig_overlays(rep, overlays, pick[:6])
    fig_survival(rep)
    fig_main(rep, overlays, [n for n in ("imatinib", "methotrexate", "sildenafil") if n in overlays])
    print((out / "summary.md").read_text())


if __name__ == "__main__":
    main()
