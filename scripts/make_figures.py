"""Generate every figure of the study from the saved results in results/ (written to figures/).

Usage: uv run --extra plot python scripts/make_figures.py
Figures whose inputs are missing are skipped with a message.
"""

from __future__ import annotations

import json
import pickle
import re
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from perlin_mole.benchmark import MOLECULES
from perlin_mole.chem import TorsionSystem
from perlin_mole.landscape import gradient_correlation
from perlin_mole.noise import FractalNoise
from perlin_mole.proposals import CorrelatedGaussianProposal, PerlinFieldProposal

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

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
        "figure.dpi": 150,
        "savefig.dpi": 600,  # raster parts (imshow) of the vector PDFs: ACS asks >= 300 dpi
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
    }
)
C_PERLIN, C_GCORR, C_IID, C_RESET, C_ETKDG = "#c0392b", "#2471a3", "#7f8c8d", "#b9770e", "#1c1c1c"
PALETTE_L = {2.0: "#1b9e77", 4.0: "#d95f02", 8.0: "#7570b3"}
NICE = {
    "ala1": "Ac-Ala-NHMe",
    "decane": "n-decane",
    "tetraglyme": "tetraglyme",
    "ala3": "Ac-(Ala)$_3$-NHMe",
    "imatinib": "imatinib",
}


_PANEL = re.compile(r"^\(([a-z])\)\s*(.*)$", re.S)


def _bold_panel_labels(fig) -> None:
    """Scientific Reports style: '(a) title' -> bold 'a' followed by the title."""
    for ax in fig.axes:
        m = _PANEL.match(ax.get_title())
        if m:
            ax.title.set_text(rf"$\mathbf{{{m.group(1)}}}$  {m.group(2)}")  # keeps size and pad


def save(fig, name: str) -> None:
    _bold_panel_labels(fig)
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png", dpi=200)
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- molecules (2-D)
def fig_molecules() -> None:
    mols, legends, hl_atoms, hl_bonds = [], [], [], []
    for name, smi in MOLECULES.items():
        s = TorsionSystem.from_smiles(smi)
        m = Chem.MolFromSmiles(smi)
        bonds = [m.GetBondBetweenAtoms(j, k).GetIdx() for _, j, k, _ in s.torsions]
        atoms = sorted({a for _, j, k, _ in s.torsions for a in (j, k)})
        mols.append(m)
        legends.append(f"{NICE[name].replace('$_3$', '3')}  ({s.n_dof} rotatable torsions)")
        hl_atoms.append(atoms)
        hl_bonds.append(bonds)
    d = rdMolDraw2D.MolDraw2DSVG(4 * 420, 330, 420, 330)
    opts = d.drawOptions()
    opts.legendFontSize = 22
    opts.highlightBondWidthMultiplier = 14
    opts.setHighlightColour((0.98, 0.72, 0.62, 1.0))
    d.DrawMolecules(mols, legends=legends, highlightAtoms=hl_atoms, highlightBonds=hl_bonds)
    d.FinishDrawing()
    svg = FIG / "molecules.svg"
    svg.write_text(d.GetDrawingText())
    for fmt in ("pdf", "png"):
        subprocess.run(["rsvg-convert", "-f", fmt, "-o", str(FIG / f"molecules.{fmt}"), str(svg)], check=True)
    print("wrote molecules")


# ---------------------------------------------------------------- Perlin field + kicks
def fig_field() -> None:
    s = TorsionSystem.from_smiles(MOLECULES["tetraglyme"])
    coords = s.dof_coordinates(3)
    L = 4.0
    rng = np.random.default_rng(7)
    fig = plt.figure(figsize=(7.0, 4.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1, 1], hspace=0.6, wspace=0.55)

    ax = fig.add_subplot(gs[:, 0])
    field = FractalNoise(2, octaves=2, period=64, seed=0)
    g = np.linspace(0, 8, 400)
    X, Y = np.meshgrid(g, g)
    Z = field(np.stack([X, Y], -1))
    ax.imshow(Z, origin="lower", extent=(0, 8, 0, 8), cmap="RdBu_r", vmin=-0.6, vmax=0.6)
    theta = 0.6
    pts = coords[:, 0] / L
    px, py = (
        1.2 + pts * np.cos(theta) - pts.min() * np.cos(theta),
        2.0 + pts * np.sin(theta) - pts.min() * np.sin(theta),
    )
    ax.plot(px, py, "-", color="k", lw=1)
    ax.scatter(
        px, py, s=22, c=field(np.stack([px, py], -1)), cmap="RdBu_r", vmin=-0.6, vmax=0.6, edgecolor="k", zorder=3
    )
    for i, off in ((0, (-14, -2)), (len(px) - 1, (4, -8))):
        ax.annotate(rf"$j={i + 1}$", (px[i], py[i]), textcoords="offset points", xytext=off, fontsize=7)
    ax.set_title("(a) fBm field $N$ read at\ntorsion sites $R\\mathbf{x}_j/L+\\mathbf{o}$")
    ax.set_xticks([])
    ax.set_yticks([])

    perlin = PerlinFieldProposal(coords, sigma=np.deg2rad(60), corr_length=L, seed=0)
    gauss = CorrelatedGaussianProposal.matched_to(perlin)
    idx = np.arange(1, s.n_dof + 1)
    for col, (prop, title, c) in enumerate(
        [(perlin, "(b) Perlin kick", C_PERLIN), (gauss, "(c) matched Gaussian kick", C_GCORR)]
    ):
        ax = fig.add_subplot(gs[0, 1 + col])
        for k in range(3):
            ax.plot(idx, np.rad2deg(prop.sample(rng)), "-o", ms=2.5, lw=1, color=c, alpha=1 - 0.3 * k)
        ax.axhline(0, color="0.6", lw=0.6)
        ax.set_ylim(-190, 190)
        ax.set_xlabel("torsion index $j$ (along chain)")
        ax.set_ylabel(r"$\Delta\phi_j$ (deg)")
        ax.set_title(title)

    ax = fig.add_subplot(gs[1, 1])
    dist = np.abs(coords[:, None, 0] - coords[None, :, 0])
    for Lc in (2.0, 4.0, 8.0):
        C = PerlinFieldProposal(coords, sigma=1.0, corr_length=Lc, seed=0).correlation()
        ds = np.arange(0, 8)
        ax.plot(ds, [C[np.isclose(dist, d)].mean() for d in ds], "-o", ms=3, color=PALETTE_L[Lc], label=f"L = {Lc:g}")
    ax.plot(ds, [0.0] * len(ds), color="0.7", lw=0.6)
    ax.set_xlabel("topological distance (bonds)")
    ax.set_ylabel(r"corr$(\Delta\phi_i,\Delta\phi_j)$")
    ax.set_title("(d) induced correlation kernel")
    ax.legend(frameon=False)

    ax = fig.add_subplot(gs[1, 2])
    cal = perlin._calibration.ravel()
    bins = np.linspace(-4, 4, 61)
    ax.hist(cal, bins=bins, density=True, color=C_PERLIN, alpha=0.55, label="Perlin")
    z = np.linspace(-4, 4, 200)
    ax.plot(z, np.exp(-z * z / 2) / np.sqrt(2 * np.pi), color=C_GCORR, lw=1.4, label="Gaussian")
    ax.set_xlabel(r"$\Delta\phi_j/\sigma$")
    ax.set_ylabel("density")
    ax.set_title("(e) marginal distribution")
    ax.set_ylim(0, 0.62)
    ax.legend(frameon=False, loc="upper center", fontsize=6)
    save(fig, "field")


# ---------------------------------------------------------------- formulation A
def fig_landscape() -> None:
    path = RES / "landscape.txt"
    if not path.exists():
        print("skip landscape")
        return
    rows = [ln.split() for ln in path.read_text().splitlines()[2:] if ln.strip()]
    n = np.array([int(r[0]) for r in rows])
    us = np.array([float(r[1]) for r in rows])
    corr = np.array([float(r[2]) if r[2] != "nan" else np.nan for r in rows])
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.3))
    ax = axes[0]
    ax.semilogy(n, us, "o-", color=C_PERLIN, ms=3, label="measured")
    ax.semilogy(n, us[5] * 2.0 ** (n - n[5]) * (n / n[5]), "--", color="0.5", lw=1, label=r"$\propto n\,2^n$")
    ax.axhline(14e3, color=C_GCORR, lw=0.8, ls=":")
    ax.text(1.2, 16e3, "one MMFF relaxation (decane)", fontsize=6.5, color=C_GCORR)
    ax.set_xlabel("number of torsions $n$")
    ax.set_ylabel(r"$\mu$s per point (value + gradient)")
    ax.set_title("(a) cost of $N(\\mathbf{q})$, $\\mathbf{q}\\in T^n$")
    ax.legend(frameon=False, loc="lower right")
    ax = axes[1]
    ax.plot(n, corr, "o-", color=C_PERLIN, ms=3, label=r"max$_{i\neq j}|$corr$(\partial_iN,\partial_jN)|$")
    ax.fill_between(n, 0, 3 / np.sqrt(4000), color="0.85", label=r"$3/\sqrt{N_s}$ sampling band")
    ax.set_ylim(0, 0.3)
    ax.set_xlabel("number of torsions $n$")
    ax.set_title("(b) gradient components uncorrelated")
    ax.legend(frameon=False, fontsize=6.5)
    ax = axes[2]
    C = gradient_correlation(6, n_points=20000)
    im = ax.imshow(C, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_title(r"(c) corr$(\partial_iN,\partial_jN)$, $n=6$")
    ax.set_xticks(range(6), range(1, 7))
    ax.set_yticks(range(6), range(1, 7))
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    save(fig, "landscape")


# ---------------------------------------------------------------- coverage curves (sigma=60)
def _report(name: str) -> dict | None:
    p = RES / name / "report.json"
    return json.loads(p.read_text()) if p.exists() else None


def fig_coverage() -> None:
    rep = _report("sigma60") or _report("main")
    if rep is None:
        print("skip coverage")
        return
    cps = rep["config"]["checkpoints"]
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.35), sharey=True)
    for ax, (mol, r) in zip(axes, rep["molecules"].items(), strict=True):
        for label, curve in r["coverage_curves"].items():
            if label in ("etkdg", "reset", "gauss-iid"):
                c = {"etkdg": C_ETKDG, "reset": C_RESET, "gauss-iid": C_IID}[label]
                ax.plot(cps, curve, color=c, lw=1.1, ls=":" if label == "gauss-iid" else "-", label=label)
            else:
                L = float(label.split("L=")[1].rstrip(")"))
                ax.plot(
                    cps,
                    curve,
                    color=PALETTE_L[L],
                    lw=1.5,
                    ls="-" if label.startswith("perlin") else "--",
                    label=label,
                )
        ax.set_title(f"{NICE[mol]} ({r['n_dof']} tors.)")
        ax.set_xlabel("MMFF relaxations")
        ax.set_ylim(0, 1)
    axes[0].set_ylabel("low-energy coverage")
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.17))
    fig.tight_layout()
    save(fig, "coverage")


# ---------------------------------------------------------------- sigma sensitivity (common reference)
def fig_sigma() -> None:
    p = RES / "discoveries" / "discoveries.json"
    if not p.exists():
        print("skip sigma")
        return
    d = json.loads(p.read_text())
    sigmas = [30, 60, 120]
    series = [
        ("gauss-iid", C_IID, ":", "gauss-iid"),
        ("gauss-corr(L=4)", PALETTE_L[4.0], "--", "gauss-corr L=4"),
        ("perlin(L=4)", PALETTE_L[4.0], "-", "Perlin L=4"),
        ("gauss-corr(L=8)", PALETTE_L[8.0], "--", "gauss-corr L=8"),
        ("perlin(L=8)", PALETTE_L[8.0], "-", "Perlin L=8"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.3))
    for ax, (mol, r) in zip(axes, d["molecules"].items(), strict=True):
        for key, c, ls, lab in series:
            pts = [
                (sg, *r["coverage_common"][f"{key} @{sg}"]) for sg in sigmas if f"{key} @{sg}" in r["coverage_common"]
            ]
            if not pts:
                continue
            sg, m, sd = map(np.array, zip(*pts, strict=True))
            ax.errorbar(sg, m, yerr=sd / np.sqrt(10), color=c, ls=ls, marker="o", ms=3, lw=1.3, capsize=2, label=lab)
        ax.set_xscale("log")
        ax.set_xticks(sigmas, [f"{s}°" for s in sigmas])
        ax.minorticks_off()
        ax.set_xlabel(r"kick size $\sigma$")
        ax.set_title(NICE[mol])
        ax.set_ylim(bottom=0)
    axes[0].set_ylabel("coverage (σ-pooled ref.)")
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.17))
    fig.tight_layout()
    save(fig, "sigma")


# ---------------------------------------------------------------- alanine dipeptide ground truth
def fig_alanine() -> None:
    p = RES / "validation" / "alanine_dipeptide.json"
    if not p.exists():
        print("skip alanine")
        return
    d = json.loads(p.read_text())
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.7), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    ph = np.array([m["phi"] for m in d["minima"]])
    ps = np.array([m["psi"] for m in d["minima"]])
    de = np.array([m["dE"] for m in d["minima"]])
    sc = ax.scatter(ph, ps, c=de, cmap="viridis", s=40, edgecolor="k", zorder=3)
    seen = set()
    for m in d["minima"]:
        if m["basin"] not in seen and m["basin"] != "other":
            seen.add(m["basin"])
            lab = {"alphaL": r"$\alpha_L$", "alpha'": r"$\alpha'$", "beta2": r"$\beta_2$", "alphaR": r"$\alpha_R$"}
            ax.annotate(
                lab.get(m["basin"], m["basin"]), (m["phi"], m["psi"]), xytext=(5, 4), textcoords="offset points"
            )
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks(range(-180, 181, 90))
    ax.set_yticks(range(-180, 181, 90))
    ax.axhline(0, color="0.85", lw=0.6)
    ax.axvline(0, color="0.85", lw=0.6)
    ax.set_xlabel(r"$\phi$ (deg)")
    ax.set_ylabel(r"$\psi$ (deg)")
    ax.set_title(f"(a) all {len(de)} MMFF94 minima ({d['n_starts']}-point grid)")
    fig.colorbar(sc, ax=ax, label=r"$\Delta E$ (kcal/mol)", fraction=0.046)
    ax = axes[1]
    names = ["etkdg", "reset", "gauss-iid", "gauss-corr", "perlin"]
    cols = [C_ETKDG, C_RESET, C_IID, C_GCORR, C_PERLIN]
    ax.bar(range(5), [d["methods"][n]["recall_mean"] for n in names], color=cols)
    ax.scatter(
        range(5), [d["methods"][n]["recall_min"] for n in names], marker="_", s=120, color="k", label="worst seed"
    )
    ax.set_xticks(range(5), names, rotation=30)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("recall of minima (≤ 5 kcal/mol)")
    ax.set_title("(b) recovery after 100 relaxations")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    save(fig, "alanine")


# ---------------------------------------------------------------- energy distributions (ETKDG vs BH)
def fig_energy_hist() -> None:
    p = RES / "sigma60" / "raw.pkl"
    if not p.exists():
        print("skip energy_hist")
        return
    runs = pickle.load(open(p, "rb"))
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.1))
    for ax, mol in zip(axes, MOLECULES, strict=True):
        es = {}
        for lab in ("etkdg", "gauss-iid", "perlin(L=4)"):
            es[lab] = np.concatenate(
                [[m[0] for m in r["minima"]] for r in runs if r["task"]["molecule"] == mol and r["label"] == lab]
            )
        e0 = min(v.min() for v in es.values())
        bins = np.linspace(0, 25, 51)
        for lab, c in (("etkdg", C_ETKDG), ("gauss-iid", C_IID), ("perlin(L=4)", C_PERLIN)):
            ax.hist(es[lab] - e0, bins=bins, histtype="step", density=True, color=c, lw=1.2, label=lab)
        ax.set_title(NICE[mol])
        ax.set_xlabel(r"$E-E^*$ (kcal/mol)")
    axes[0].set_ylabel("density of unique minima")
    axes[-1].legend(frameon=False)
    fig.tight_layout()
    save(fig, "energy_hist")


# ---------------------------------------------------------------- 3-D renders
CPK = {"C": "#404040", "O": "#e03b2c", "N": "#3050f8", "H": "#e8e8e8"}
RAD = {"C": 0.36, "O": 0.34, "N": 0.34, "H": 0.2}


def _hbonds(mol: Chem.Mol, xyz: np.ndarray) -> list[tuple[int, int]]:
    out = []
    for h in mol.GetAtoms():
        if h.GetSymbol() != "H":
            continue
        d = h.GetNeighbors()[0]
        if d.GetSymbol() not in ("N", "O"):
            continue
        for a in mol.GetAtoms():
            if a.GetSymbol() != "O" or a.GetIdx() == d.GetIdx():
                continue
            r = np.linalg.norm(xyz[h.GetIdx()] - xyz[a.GetIdx()])
            v1, v2 = xyz[d.GetIdx()] - xyz[h.GetIdx()], xyz[a.GetIdx()] - xyz[h.GetIdx()]
            ang = np.degrees(np.arccos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))
            if r < 2.5 and ang > 120:
                out.append((h.GetIdx(), a.GetIdx()))
    return out


def render(ax, mol: Chem.Mol, title: str = "", hide_h: bool = False) -> None:
    xyz = mol.GetConformer().GetPositions()
    xyz = xyz - xyz.mean(0)
    _, _, vt = np.linalg.svd(xyz, full_matrices=False)
    xyz = xyz @ vt.T  # principal axes: x longest, z = depth
    z = xyz[:, 2]
    zn = (z - z.min()) / (np.ptp(z) + 1e-9)
    keep = [a.GetIdx() for a in mol.GetAtoms() if not (hide_h and a.GetSymbol() == "H")]
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i not in keep or j not in keep:
            continue
        mid = (xyz[i] + xyz[j]) / 2
        for a, p in ((i, xyz[i]), (j, xyz[j])):
            sym = mol.GetAtomWithIdx(a).GetSymbol()
            ax.plot(
                [p[0], mid[0]],
                [p[1], mid[1]],
                color=CPK.get(sym, "#a0a0a0"),
                lw=1.2 + 1.6 * zn[a],
                solid_capstyle="round",
                zorder=2 + 10 * zn[a],
            )
    for h, a in _hbonds(mol, xyz):
        ax.plot([xyz[h, 0], xyz[a, 0]], [xyz[h, 1], xyz[a, 1]], ls=(0, (1.5, 1.5)), color="#2e86c1", lw=1.0, zorder=30)
    for a in sorted(keep, key=lambda k: z[k]):
        sym = mol.GetAtomWithIdx(a).GetSymbol()
        ax.add_patch(
            Circle(
                xyz[a, :2],
                RAD.get(sym, 0.3) * (0.75 + 0.45 * zn[a]),
                facecolor=CPK.get(sym, "#a0a0a0"),
                edgecolor="k",
                lw=0.4,
                zorder=3 + 10 * zn[a],
            )
        )
    pad = 0.8
    ax.set_xlim(xyz[:, 0].min() - pad, xyz[:, 0].max() + pad)
    ax.set_ylim(xyz[:, 1].min() - pad, xyz[:, 1].max() + pad)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=7.5)


def _sdf(mol: str) -> list[Chem.Mol]:
    p = RES / "discoveries" / "sdf" / f"{mol}_low_energy.sdf"
    return [m for m in Chem.SDMolSupplier(str(p), removeHs=False)] if p.exists() else []


def fig_global_minima() -> None:
    mols = {m: _sdf(m) for m in MOLECULES}
    if not all(mols.values()):
        print("skip global_minima")
        return
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.2), gridspec_kw={"width_ratios": [1.6, 1.2, 1, 1.5]})
    for ax, (name, ms) in zip(axes, mols.items(), strict=True):
        m = ms[0]
        grid = m.GetProp("in_grid_reference")
        render(ax, m, f"{NICE[name]}\n$E^*$ = {float(m.GetProp('E_MMFF94_kcal')):.2f} kcal/mol (grid: {grid})")
    fig.tight_layout()
    save(fig, "global_minima")


def fig_gallery(mol: str, n: int = 8) -> None:
    ms = _sdf(mol)
    if not ms:
        print("skip gallery", mol)
        return
    # skip mirror-image duplicates (identical energy) to show structural variety
    picked, last = [], None
    for m in ms:
        e = float(m.GetProp("E_MMFF94_kcal"))
        if last is not None and abs(e - last) < 1e-3:
            continue
        picked.append(m)
        last = e
        if len(picked) == n:
            break
    cols = 4
    rows = int(np.ceil(len(picked) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(7.2, 1.9 * rows))
    for ax, m in zip(axes.ravel(), picked, strict=False):
        tag = {"yes": "", "no": " ★ beyond grid", "n/a": ""}[m.GetProp("in_grid_reference")]
        render(
            ax,
            m,
            f"ΔE = {float(m.GetProp('dE_kcal')):.2f}{tag}\nfound in {m.GetProp('n_runs_found')} runs",
            hide_h=mol != "ala3",
        )
    for ax in axes.ravel()[len(picked) :]:
        ax.axis("off")
    fig.tight_layout()
    save(fig, f"gallery_{mol}")


# ---------------------------------------------------------------- ablation forest + time-to-known survival
def _ablation(path: Path | None = None) -> dict | None:
    p = path or (RES / "ablation" / "ablation.json")
    return json.loads(p.read_text()) if p.exists() else None


def fig_ablation_forest(path: Path | None = None) -> None:
    d = _ablation(path)
    if d is None:
        print("skip ablation_forest")
        return
    tests = d["tests"]
    qs = list(dict.fromkeys((t["question"], t["variant"], t["reference"]) for t in tests))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 0.32 * len(qs) + 1.0), sharey=True)
    for ax, metric, xl in (
        (axes[0], "recall_final", r"$\Delta$ final recall of known minima"),
        (axes[1], "t_global", r"$\Delta$ steps to known global min"),
    ):
        for i, key in enumerate(qs):
            t = next(t for t in tests if (t["question"], t["variant"], t["reference"]) == key and t["metric"] == metric)
            good = (t["mean_diff"] > 0) if metric == "recall_final" else (t["mean_diff"] < 0)
            sig = t["p_holm"] < 0.05
            c = (C_PERLIN if "perlin" in key[1] else C_GCORR) if sig else "0.55"
            ax.errorbar(
                t["mean_diff"],
                -i,
                xerr=[[t["mean_diff"] - t["ci95"][0]], [t["ci95"][1] - t["mean_diff"]]],
                fmt="o" if sig else "o",
                mfc=c if sig else "white",
                color=c,
                ms=4,
                capsize=2,
                lw=1,
            )
            if sig:
                ax.annotate(
                    "*" if good else "†", (t["ci95"][1], -i), xytext=(3, -3), textcoords="offset points", fontsize=8
                )
        ax.axvline(0, color="k", lw=0.6)
        ax.set_xlabel(xl)
        ax.grid(axis="x", alpha=0.3)
    axes[0].set_yticks([-i for i in range(len(qs))], [f"{q}\n({a} − {b})" for q, a, b in qs], fontsize=6)
    fig.suptitle(
        "Ablation: paired differences over molecule × seed (95% bootstrap CI; filled = Holm p < 0.05)", fontsize=8
    )
    fig.tight_layout()
    save(fig, "ablation_forest")


def fig_known_survival(path: Path | None = None) -> None:
    d = _ablation(path)
    if d is None:
        print("skip known_survival")
        return
    show = [
        ("reset", C_RESET, "-"),
        ("gauss-iid", C_IID, ":"),
        ("gauss-corr(L=4)", C_GCORR, "--"),
        ("perlin(L=4)", C_PERLIN, "-"),
        ("gauss-ar1(L=4)", "#5b2c6f", "--"),
        ("perlin-drift(L=4)", "#e67e22", "-"),
    ]
    mols = list(d["summary"])
    b = d["budget"]
    fig, axes = plt.subplots(1, len(mols), figsize=(7.2, 2.2), sharey=True)
    for ax, mol in zip(np.atleast_1d(axes), mols, strict=True):
        for lab, c, ls in show:
            if lab not in d["summary"][mol]:
                continue
            surv = d["summary"][mol][lab]["survival"]
            xs = np.linspace(0, b, len(surv))
            ax.step(xs, surv, where="post", color=c, ls=ls, lw=1.3, label=lab)
        ax.set_title(f"{NICE.get(mol, mol)} ({d['targets'][mol]['n_known_low']} known)", fontsize=7.5)
        ax.set_xlabel("relaxations")
        ax.set_ylim(0, 1.02)
    np.atleast_1d(axes)[0].set_ylabel("P(known global min\nnot yet found)")
    h, lab = np.atleast_1d(axes)[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=6, frameon=False, bbox_to_anchor=(0.5, -0.12), fontsize=6.5)
    fig.tight_layout()
    save(fig, "known_survival")


# ---------------------------------------------------------------- isomers / enantiomers / constitutional
def fig_isomers() -> None:
    from perlin_mole.conformational import compare

    iso = RES / "isomers"
    if not (iso / "isomers.json").exists():
        print("skip isomers")
        return
    d = json.loads((iso / "isomers.json").read_text())
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw={"width_ratios": [1.2, 1, 1.3]})
    ax = axes[0]
    rows = compare(iso / "grids")
    lab = [f"{r['molecule'].replace('1,2-', '').replace('ethane', '')}\n{r['quantity']}" for r in rows]
    y = np.arange(len(rows))
    ax.errorbar(
        [r["experiment"] for r in rows],
        y,
        xerr=[r["err"] for r in rows],
        fmt="s",
        color="k",
        ms=4,
        capsize=2,
        label="experiment (gas)",
    )
    ax.scatter([r["mmff94"] for r in rows], y, color=C_PERLIN, marker="o", s=22, zorder=3, label="MMFF94 (grid)")
    ax.axvline(0, color="0.7", lw=0.6)
    ax.set_yticks(y, lab, fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel(r"$\Delta E$ (kcal/mol)")
    ax.set_title("(a) conformational isomers")
    ax.legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
    ax = axes[1]
    pairs = d.get("enantiomers", {}).get("enantiomer_pairs", {})
    names = list(pairs)
    x = np.arange(len(names))
    ax.bar(x - 0.2, [pairs[n]["mirror_match_a_in_b"] for n in names], 0.4, color=C_GCORR, label="search: R in S")
    ax.bar(x + 0.2, [pairs[n]["grid_mirror_match"] for n in names], 0.4, color="0.6", label="grid")
    ax.set_xticks(x, [n.split(" / ")[0].rsplit("-", 1)[0] for n in names], rotation=30, fontsize=6)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("mirror-matched fraction")
    ax.set_title("(b) enantiomer consistency")
    ax.legend(frameon=True, framealpha=0.95, edgecolor="none", fontsize=6, loc="lower left")
    ax = axes[2]
    cs = d.get("constitutional", {}).get("ablation_summary", {})
    mols = list(cs)
    meths = [("reset", C_RESET), ("gauss-iid", C_IID), ("gauss-corr(L=4)", C_GCORR), ("perlin(L=4)", C_PERLIN)]
    x = np.arange(len(mols))
    for i, (m, c) in enumerate(meths):
        ax.bar(
            x + (i - 1.5) * 0.2,
            [cs[mol][m]["recall_final"] if m in cs[mol] else 0 for mol in mols],
            0.2,
            color=c,
            label=m,
        )
    short = {
        "n-decane": "n-C10",
        "2-methylnonane": "2-Me",
        "2,7-dimethyloctane": "2,7-Me2",
        "4-propylheptane": "4-Pr",
        "3,3-diethylhexane": "3,3-Et2",
    }
    ax.set_xticks(x, [short.get(m, m) for m in mols], fontsize=6.5)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("recall of known minima")
    ax.set_title(r"(c) C$_{10}$H$_{22}$ isomers")
    ax.legend(frameon=True, framealpha=0.95, edgecolor="none", fontsize=5.5, ncol=2, loc="lower left")
    fig.tight_layout()
    save(fig, "isomers")


# ---------------------------------------------------------------- regime hunt
def fig_hunt() -> None:
    reps = {
        sg: json.loads((RES / "hunt" / f"sigma{sg}" / "report.json").read_text())
        for sg in (60, 120)
        if (RES / "hunt" / f"sigma{sg}" / "report.json").exists()
    }
    if not reps:
        print("skip hunt")
        return
    order = [
        ("reset", C_RESET),
        ("gauss-iid", C_IID),
        ("gauss-corr(L=8)", C_GCORR),
        ("perlin(L=8)", C_PERLIN),
        ("gauss-anneal(L=8)", "#1b9e77"),
        ("perlin-anneal(L=8)", "#66c2a5"),
        ("gauss-ar1(L=8)", "#5b2c6f"),
        ("perlin-drift(L=8)", "#e67e22"),
    ]
    mols = list(next(iter(reps.values()))["molecules"])
    fig, axes = plt.subplots(len(reps), len(mols), figsize=(7.2, 2.1 * len(reps)), squeeze=False)
    for r_i, (sg, rep) in enumerate(reps.items()):
        for c_i, mol in enumerate(mols):
            ax = axes[r_i][c_i]
            sm = rep["molecules"][mol]["summary"]
            vals = [(lab, c, sm[lab]["coverage_mean"], sm[lab]["coverage_std"]) for lab, c in order if lab in sm]
            ax.bar(
                range(len(vals)),
                [v[2] for v in vals],
                yerr=[v[3] / np.sqrt(10) for v in vals],
                color=[v[1] for v in vals],
                capsize=1.5,
            )
            ax.set_xticks(range(len(vals)), [v[0].replace("(L=8)", "") for v in vals], rotation=60, fontsize=5)
            ax.set_title(f"{mol} ({rep['molecules'][mol]['n_dof']} tors.), $\\sigma$={sg}°", fontsize=7)
            if c_i == 0:
                ax.set_ylabel("coverage")
    fig.tight_layout()
    save(fig, "hunt")


# ---------------------------------------------------------------- 2-D bias (formulation A1)
def fig_bias2d() -> None:
    p = RES / "bias2d" / "bias2d.json"
    if not p.exists():
        print("skip bias2d")
        return
    from perlin_mole.bias2d import PerlinBias, Surface

    d = json.loads(p.read_text())
    grid = np.load(RES / "bias2d" / "surface.npy")
    surf = Surface(grid, 360.0 / grid.shape[0])
    a = np.linspace(-180, 180, 181)
    P, Q = np.meshgrid(a, a, indexing="ij")
    E = surf(P, Q)
    B = PerlinBias(4.0, cells=4, seed=0)(P, Q)
    fig, axes = plt.subplots(1, 4, figsize=(7.4, 2.0), gridspec_kw={"width_ratios": [1, 1, 1.1, 1.1]})
    for ax, Z, t in ((axes[0], E, "(a) MMFF94 $E(\\phi,\\psi)$"), (axes[1], E + B, "(b) $E+\\lambda N$, $\\lambda=4$")):
        im = ax.contourf(P, Q, np.minimum(Z, 15), levels=16, cmap="viridis")
        ax.set_xlabel(r"$\phi$")
        ax.set_ylabel(r"$\psi$")
        ax.set_title(t, fontsize=7)
        ax.set_xticks([-180, 0, 180])
        ax.set_yticks([-180, 0, 180])
    fig.colorbar(im, ax=axes[1], fraction=0.05, label="kcal/mol")
    s = d["summary"]
    lams = [2.0, 4.0, 6.0]
    for ax, key, yl in (
        (axes[2], "transitions_per_1e5", "transitions / $10^5$ steps"),
        (axes[3], "dF_rmse", r"$\Delta F$ RMSE (kcal/mol)"),
    ):
        for kind, c, ls in (
            ("perlin", C_PERLIN, "-"),
            ("grf", C_GCORR, "--"),
            ("perlin-drift", "#e67e22", "-"),
            ("grf-drift", "#5b2c6f", "--"),
        ):
            ys = [s.get(f"{kind}@{lam:g}", {}).get(key) for lam in lams]
            ax.plot(lams, [np.nan if v is None else v for v in ys], ls=ls, marker="o", ms=3, color=c, label=kind)
        for ref, c in (("none@0", "k"), ("metad@0", "#1b9e77")):
            v = s.get(ref, {}).get(key)
            if v is not None:
                ax.axhline(v, color=c, lw=0.8, ls=":", label=ref.split("@")[0])
        ax.set_xlabel(r"bias strength $\lambda$ (kcal/mol)")
        ax.set_ylabel(yl)
    axes[2].set_title("(c) barrier crossings", fontsize=7)
    axes[3].set_title(f"(d) free-energy error (exact {d['exact_dF']:.2f})", fontsize=7)
    axes[3].legend(frameon=False, fontsize=5.5)
    fig.tight_layout()
    save(fig, "bias2d")


# ---------------------------------------------------------------- Cartesian kicks
def fig_cartesian() -> None:
    p = RES / "cartesian" / "raw.pkl"
    a = RES / "cartesian" / "ablation.json"
    if not (p.exists() and a.exists()):
        print("skip cartesian")
        return
    runs = pickle.load(open(p, "rb"))
    abl = json.loads(a.read_text())
    meths = ["cart-iid", "cart-perlin@2", "cart-gauss-field@2", "cart-perlin@4", "cart-gauss-field@4"]
    cols = [C_IID, C_PERLIN, C_GCORR, "#e6a39b", "#8fb4d6"]
    mols = list(abl["summary"])
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.3))
    x = np.arange(len(mols))
    for i, (m, c) in enumerate(zip(meths, cols, strict=True)):
        integ = [
            1
            - np.mean([r["broken"] for r in runs if r["task"]["molecule"] == mol and r["label"] == m])
            / runs[0]["task"]["n_steps"]
            for mol in mols
        ]
        rec = [abl["summary"][mol].get(m, {}).get("recall_final", 0) for mol in mols]
        axes[0].bar(x + (i - 2) * 0.16, integ, 0.16, color=c, label=m.replace("cart-", ""))
        axes[1].bar(x + (i - 2) * 0.16, rec, 0.16, color=c)
    for ax, yl, t in (
        (axes[0], "chemically intact relaxations", "(a) integrity"),
        (axes[1], "recall of known minima", "(b) recall"),
    ):
        ax.set_xticks(x, [NICE.get(m, m) for m in mols], fontsize=6)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel(yl)
        ax.set_title(t)
    axes[0].legend(frameon=False, fontsize=6, ncol=2, loc="lower left")
    fig.tight_layout()
    save(fig, "cartesian")


# ---------------------------------------------------------------- DFT validation (GPU4PySCF)
def fig_dft() -> None:
    p = RES / "dft" / "dft.json"
    if not p.exists():
        print("skip dft")
        return
    d = json.loads(p.read_text())
    p2 = RES / "dft" / "dft_d2.json"
    if p2.exists():
        d["d2"] = {**d.get("d2", {}), **json.loads(p2.read_text()).get("d2", {})}
    mlj = (
        json.loads((RES / "ml_judge" / "ml_judge.json").read_text())
        if (RES / "ml_judge" / "ml_judge.json").exists()
        else {}
    )
    ani = {(r["molecule"], r["quantity"]): r["ani2x"] for r in mlj.get("conformational", [])}
    mols = list(d.get("d2", {}))
    fig, axes = plt.subplots(
        1, 1 + len(mols), figsize=(7.4, 2.5), gridspec_kw={"width_ratios": [1.6] + [1] * len(mols)}
    )
    ax = np.atleast_1d(axes)[0]
    rows = d.get("d1", [])
    y = np.arange(len(rows))
    ax.errorbar(
        [r["experiment"] for r in rows],
        y,
        xerr=[r["err"] for r in rows],
        fmt="s",
        color="k",
        ms=4,
        capsize=2,
        label="experiment",
    )
    ax.scatter([r["mmff94"] for r in rows], y, color=C_PERLIN, s=16, zorder=3, label="MMFF94")
    ya = [i for i, r in enumerate(rows) if (r["molecule"], r["quantity"]) in ani]
    ax.scatter(
        [ani[(rows[i]["molecule"], rows[i]["quantity"])] for i in ya],
        ya,
        color=C_GCORR,
        marker="D",
        s=14,
        zorder=3,
        label="ANI-2x",
    )
    ax.scatter([r["dft"] for r in rows], y, color="#1b9e77", marker="^", s=22, zorder=4, label="DFT")
    ax.set_yticks(y, [f"{r['molecule'].replace('1,2-', '')}\n{r['quantity']}" for r in rows], fontsize=5.5)
    ax.invert_yaxis()
    ax.axvline(0, color="0.85", lw=0.6)
    ax.set_xlabel(r"$\Delta E$ (kcal/mol)", fontsize=7)
    ax.set_title("(a) predicted vs real", fontsize=7.5)
    ax.legend(frameon=False, fontsize=5.5, loc="lower right")
    for ax, mol in zip(np.atleast_1d(axes)[1:], mols, strict=True):
        v = d["d2"][mol]
        x = np.array(v["dE_dft"])
        ax.scatter(x, v["dE_mmff"], s=6, color=C_PERLIN, alpha=0.7, label="MMFF94")
        if v["dE_ani"] and not np.isnan(v["dE_ani"][0]):
            ax.scatter(x, v["dE_ani"], s=6, color=C_GCORR, alpha=0.7, marker="D", label="ANI-2x")
        lim = max(x.max(), max(v["dE_mmff"])) * 1.05
        ax.plot([0, lim], [0, lim], color="0.6", lw=0.6)
        rho_a = v.get("spearman_dft_ani")
        ax.set_title(
            f"{NICE.get(mol, mol)}\n"
            + rf"$\rho_M$={v['spearman_dft_mmff']:.2f}"
            + (rf", $\rho_A$={rho_a:.2f}" if rho_a is not None else ""),
            fontsize=6.5,
        )
        ax.set_xlabel(r"DFT $\Delta E$", fontsize=6.5)
        ax.tick_params(labelsize=6)
    np.atleast_1d(axes)[1].set_ylabel(r"judge $\Delta E$", fontsize=6.5)
    np.atleast_1d(axes)[1].legend(frameon=False, fontsize=5.5, loc="upper left")
    fig.tight_layout()
    save(fig, "dft")


# ---------------------------------------------------------------- LJ38: exploratory -> pre-registered tests
def _wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def fig_lj38() -> None:
    ex = RES / "clusters" / "clusters.json"
    cf, rp = RES / "clusters_confirm" / "result.json", RES / "clusters_replication" / "result.json"
    if not (ex.exists() and cf.exists() and rp.exists()):
        print("skip lj38")
        return
    m = json.loads(ex.read_text())["clusters"]["38"]["methods"]
    stages = [
        ("exploratory\n(seeds 0-29)", {k: (round(v["success"] * v["n_runs"]), v["n_runs"]) for k, v in m.items()})
    ]
    for lab, f in (("confirmation$^\\dagger$\n(seeds 30-129)", cf), ("replication$^\\dagger$\n(seeds 130-429)", rp)):
        r = json.loads(f.read_text())
        stages.append((lab, {k: (v, r["n"]) for k, v in r["successes"].items()}))
    arms = [
        ("perlin@2", "Perlin field", C_PERLIN, "o"),
        ("gauss-field@2", "matched Gaussian field", C_GCORR, "s"),
        ("gauss-iid", "i.i.d. Gaussian", "0.45", "v"),
        ("iid-perlinmarg", "i.i.d., Perlin marginal", "0.2", "^"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.6), gridspec_kw={"width_ratios": [1.6, 1]})
    ax = axes[0]
    for j, (key, lab, col, mk) in enumerate(arms):
        xs, ys, lo, hi = [], [], [], []
        for i, (_, d) in enumerate(stages):
            if key in d:
                k, n = d[key]
                a, b = _wilson(k, n)
                xs.append(i + (j - 1.5) * 0.13)
                ys.append(k / n)
                lo.append(k / n - a)
                hi.append(b - k / n)
        ax.errorbar(xs, ys, yerr=[lo, hi], fmt=mk, color=col, ms=4, capsize=2, lw=0.9, label=lab)
    ax.set_xticks(range(len(stages)), [s[0] for s in stages], fontsize=6.5)
    ax.set_ylabel("runs reaching the LJ$_{38}$\nglobal minimum", fontsize=7)
    ax.set_ylim(0, 0.85)
    ax.legend(frameon=False, fontsize=6, ncol=2, loc="upper right")
    ax.set_title("a", loc="left", fontweight="bold", fontsize=9)
    ax = axes[1]
    ns, ds, cis = [], [], []
    for _, d in stages:
        (k1, n), (k2, _) = d["perlin@2"], d["gauss-field@2"]
        p1, p2 = k1 / n, k2 / n
        se = np.sqrt(p1 * (1 - p1) / n + p2 * (1 - p2) / n)
        ns.append(n)
        ds.append(p1 - p2)
        cis.append(1.96 * se)
    ax.errorbar(ns, ds, yerr=cis, fmt="o-", color=C_PERLIN, ms=4, capsize=2, lw=0.9)
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xscale("log")
    ax.set_xticks(ns, [str(n) for n in ns], fontsize=6.5)
    ax.minorticks_off()
    ax.set_xlabel("seeds per arm", fontsize=7)
    ax.set_ylabel("Perlin $-$ matched field\n(success rate)", fontsize=7)
    ax.set_title("b", loc="left", fontweight="bold", fontsize=9)
    for a in axes:
        a.tick_params(labelsize=6.5)
    fig.tight_layout()
    save(fig, "lj38")


def fig_dft_main() -> None:
    """Scientific Reports Fig. 7: (a) predicted vs experimental conformational energies for three judges;
    (b) def2-TZVP energy of each judge's global minimum relative to the def2-TZVP global minimum."""
    p, p3 = RES / "dft" / "dft.json", RES / "dft" / "dft_d3.json"
    if not (p.exists() and p3.exists()):
        print("skip dft_main")
        return
    rows = json.loads(p.read_text()).get("d1", [])
    d3 = json.loads(p3.read_text())["d3"]
    mlj = (
        json.loads((RES / "ml_judge" / "ml_judge.json").read_text())
        if (RES / "ml_judge" / "ml_judge.json").exists()
        else {}
    )
    ani = {(r["molecule"], r["quantity"]): r["ani2x"] for r in mlj.get("conformational", [])}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    y = np.arange(len(rows))
    ax.errorbar(
        [r["experiment"] for r in rows],
        y,
        xerr=[r["err"] for r in rows],
        fmt="s",
        color="k",
        ms=4,
        capsize=2,
        label="experiment",
    )
    ax.scatter([r["mmff94"] for r in rows], y, color=C_PERLIN, s=16, zorder=3, label="MMFF94")
    ya = [i for i, r in enumerate(rows) if (r["molecule"], r["quantity"]) in ani]
    ax.scatter(
        [ani[(rows[i]["molecule"], rows[i]["quantity"])] for i in ya],
        ya,
        color=C_GCORR,
        marker="D",
        s=14,
        zorder=3,
        label="ANI-2x",
    )
    ax.scatter([r["dft"] for r in rows], y, color="#1b9e77", marker="^", s=22, zorder=4, label="DFT")
    ax.set_yticks(y, [f"{r['molecule']}: {r['quantity']}" for r in rows], fontsize=6.5)
    ax.invert_yaxis()
    ax.axvline(0, color="0.85", lw=0.6)
    ax.set_xlabel(r"conformational $\Delta E$ (kcal/mol)")
    ax.set_title("(a) predicted vs experiment")
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    ax = axes[1]
    mols = [m for m in ("decane", "tetraglyme", "ala3", "imatinib") if m in d3]
    judges = [
        ("mmff_gm_dE_tzvp", "MMFF94", C_PERLIN),
        ("ani_gm_dE_tzvp", "ANI-2x", C_GCORR),
        ("svp_gm_dE_tzvp", "DFT/def2-SVP", "#1b9e77"),
    ]
    w = 0.26
    for j, (key, lab, col) in enumerate(judges):
        vals = [d3[m][key] if d3[m][key] is not None else np.nan for m in mols]
        xs = np.arange(len(mols)) + (j - 1) * w
        ax.bar(xs, vals, w, color=col, label=lab)
        for x, v in zip(xs, vals, strict=True):
            if v == 0:
                ax.text(x, 0.02, "✓", ha="center", va="bottom", fontsize=7, color=col)
    short = {"decane": "n-decane", "tetraglyme": "tetraglyme", "ala3": "Ac-(Ala)$_3$-\nNHMe", "imatinib": "imatinib"}
    ax.set_xticks(range(len(mols)), [short.get(m, m) for m in mols], fontsize=6.5)
    ax.set_ylabel("judge's global minimum:\nkcal/mol above DFT/def2-TZVP minimum", fontsize=7)
    ax.set_title("(b) who finds the DFT global minimum?")
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    fig.tight_layout()
    save(fig, "dft_main")


def fig_concept() -> None:
    """Concept figure, panels a-b: a Perlin field and the kicks it induces on a molecule.

    a: formulation A on two torsions -- a one-octave Perlin field over the torus with its gradient (drift)
       arrows, and the joint distribution of the two gradient components (Proposition 1).
    b: formulation B on tetraglyme -- one kick drawn from Perlin (L = 4 bonds), its matched Gaussian and
       independent kicks, each rotatable bond coloured by its torsion change.
    """
    from matplotlib.collections import LineCollection
    from matplotlib.colors import TwoSlopeNorm
    from rdkit.Chem import AllChem as AC

    from perlin_mole.chem import TorsionSystem
    from perlin_mole.noise.perlin import FractalNoise
    from perlin_mole.proposals import build_proposal

    fig = plt.figure(figsize=(7.2, 3.15))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.75], wspace=0.08)

    # ---- a: field over two torsion angles + drift arrows
    gsa = gs[0].subgridspec(1, 1)
    ax = fig.add_subplot(gsa[0])
    period = 4
    f = FractalNoise(2, octaves=1, period=period, seed=3)
    g = np.linspace(0, 360, 181)
    P, Q = np.meshgrid(g, g, indexing="xy")
    pts = np.stack([P, Q], -1) * period / 360.0
    v = f(pts)
    ax.imshow(
        v,
        origin="lower",
        extent=(0, 360, 0, 360),
        cmap="RdBu_r",
        vmin=-np.abs(v).max(),
        vmax=np.abs(v).max(),
        interpolation="bilinear",
    )
    c = np.linspace(15, 345, 12)
    CP, CQ = np.meshgrid(c, c, indexing="xy")
    _, grad = f(np.stack([CP, CQ], -1) * period / 360.0, grad=True)
    ax.quiver(CP, CQ, -grad[..., 0], -grad[..., 1], color="k", scale=28, width=0.006, headwidth=4)
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 360)
    ax.set_xticks([0, 180, 360])
    ax.set_yticks([0, 180, 360])
    ax.set_xlabel(r"torsion $\phi_1$ (deg)")
    ax.set_ylabel(r"torsion $\phi_2$ (deg)")
    ax.set_title("(a) formulation A: $-\\nabla N$ over torsion space", fontsize=7.5)
    # Proposition 1 is a statement about the ensemble of fields: pool gradients over 2000 independent fields
    rng = np.random.default_rng(0)
    gr = np.concatenate(
        [
            FractalNoise(2, octaves=1, period=period, seed=100 + k)(rng.uniform(0, period, (20, 2)), grad=True)[1]
            for k in range(2000)
        ]
    )
    r = float(np.corrcoef(gr[:, 0], gr[:, 1])[0, 1])
    ins = ax.inset_axes([0.69, 0.04, 0.28, 0.28])
    ins.scatter(gr[::8, 0], gr[::8, 1], s=0.3, color="0.15", alpha=0.35, lw=0)  # r uses all points
    ins.set_xticks([])
    ins.set_yticks([])
    ins.set_facecolor("white")
    ins.text(
        0.04,
        0.96,
        f"$\\partial_1N$ vs $\\partial_2N$\n$r={r:+.2f}$",
        transform=ins.transAxes,
        fontsize=5,
        va="top",
        ha="left",
        bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.8),
    )

    # ---- b: tetraglyme coloured by three kicks
    smi = "COCCOCCOCCOCCOC"
    sysm = TorsionSystem.from_smiles(smi)
    heavy = Chem.Mol(sysm.heavy_mol())
    AC.Compute2DCoords(heavy)
    xy = heavy.GetConformer().GetPositions()[:, :2]
    xy = xy - xy.mean(0)
    if np.ptp(xy[:, 1]) > np.ptp(xy[:, 0]):
        xy = xy[:, ::-1]
    coords = sysm.dof_coordinates()
    draws = {}
    for meth in ("perlin", "gauss-corr", "gauss-iid"):
        prop = build_proposal(meth, coords, 60.0, 4.0, seed=1)
        draws[meth] = prop.sample(np.random.default_rng(11))
    central = {frozenset((j, k)): t for t, (_, j, k, _) in enumerate(sysm.torsions)}
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-150, vmax=150)
    cmap = plt.get_cmap("RdBu_r")
    gsb = gs[1].subgridspec(3, 1, hspace=0.05)
    labels = {
        "perlin": "Perlin kick ($L=4$ bonds)",
        "gauss-corr": "matched Gaussian kick",
        "gauss-iid": "independent (i.i.d.) kick",
    }
    for row, meth in enumerate(("perlin", "gauss-corr", "gauss-iid")):
        ax = fig.add_subplot(gsb[row])
        thin, thick, cols = [], [], []
        for b in heavy.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            seg = [xy[i], xy[j]]
            t = central.get(frozenset((i, j)))
            if t is None:
                thin.append(seg)
            else:
                thick.append(seg)
                cols.append(cmap(norm(draws[meth][t])))
        ax.add_collection(LineCollection(thin, colors="0.55", linewidths=1.0, zorder=1))
        ax.add_collection(LineCollection(thick, colors=cols, linewidths=4.2, capstyle="round", zorder=2))
        for a in heavy.GetAtoms():
            if a.GetSymbol() == "O":
                ax.text(
                    *xy[a.GetIdx()],
                    "O",
                    ha="center",
                    va="center",
                    fontsize=6.5,
                    color="#c0392b",
                    fontweight="bold",
                    zorder=3,
                    bbox=dict(boxstyle="circle,pad=0.08", fc="white", ec="none"),
                )
        ax.set_aspect("equal")
        ax.set_xlim(xy[:, 0].min() - 0.8, xy[:, 0].max() + 0.8)
        ax.set_ylim(xy[:, 1].min() - 0.9, xy[:, 1].max() + 0.9)
        ax.axis("off")
        ax.text(0.0, 1.0, labels[meth], transform=ax.transAxes, fontsize=6.5, va="bottom", ha="left")
        if row == 0:
            ax.set_title("(b) formulation B on tetraglyme: one kick, coloured per torsion", fontsize=7.5, pad=20)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cax = fig.add_axes([0.93, 0.2, 0.012, 0.6])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label(r"$\Delta\phi_j$ (deg)", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    save(fig, "fig1_panels")


def fig_toc() -> None:
    """JCIM table-of-contents graphic: 3.25 x 1.75 in, sans-serif text >= 6 pt (ACS TOC guidelines)."""
    from matplotlib.collections import LineCollection
    from matplotlib.colors import TwoSlopeNorm
    from rdkit.Chem import AllChem as AC

    from perlin_mole.chem import TorsionSystem
    from perlin_mole.noise.perlin import FractalNoise
    from perlin_mole.proposals import build_proposal

    with plt.rc_context(
        {
            "savefig.bbox": "standard",
            "font.family": "sans-serif",
            "font.sans-serif": ["Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
        }
    ):
        fig = plt.figure(figsize=(3.25, 1.75))
        # left: a Perlin field patch
        ax_f = fig.add_axes([0.02, 0.30, 0.22, 0.42])
        f = FractalNoise(2, octaves=2, period=4, seed=3)
        g = np.linspace(0, 4, 120)
        P, Q = np.meshgrid(g, g)
        ax_f.imshow(f(np.stack([P, Q], -1)), cmap="RdBu_r", origin="lower")
        ax_f.set_xticks([])
        ax_f.set_yticks([])
        ax_f.set_title("Perlin field", fontsize=7, pad=2)
        for sp in ax_f.spines.values():
            sp.set_visible(True)
            sp.set_linewidth(0.6)
        # molecules: Perlin kick (top) and matched Gaussian kick (bottom) on tetraglyme
        s = TorsionSystem.from_smiles("COCCOCCOCCOCCOC")
        heavy = Chem.Mol(s.heavy_mol())
        AC.Compute2DCoords(heavy)
        xy = heavy.GetConformer().GetPositions()[:, :2]
        xy = xy - xy.mean(0)
        if np.ptp(xy[:, 1]) > np.ptp(xy[:, 0]):
            xy = xy[:, ::-1]
        coords = s.dof_coordinates()
        central = {frozenset((j, k)): t for t, (_, j, k, _) in enumerate(s.torsions)}
        norm = TwoSlopeNorm(vcenter=0.0, vmin=-150, vmax=150)
        cmap = plt.get_cmap("RdBu_r")
        rows = (("perlin", "Perlin kick", 0.56), ("gauss-corr", "matched Gaussian kick", 0.08))
        for meth, lab, y0 in rows:
            d = build_proposal(meth, coords, 60.0, 4.0, seed=1).sample(np.random.default_rng(11))
            ax = fig.add_axes([0.33, y0, 0.66, 0.36])
            thin, thick, cols = [], [], []
            for b in heavy.GetBonds():
                i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                t = central.get(frozenset((i, j)))
                (thin if t is None else thick).append([xy[i], xy[j]])
                if t is not None:
                    cols.append(cmap(norm(d[t])))
            ax.add_collection(LineCollection(thin, colors="0.5", linewidths=0.8))
            ax.add_collection(LineCollection(thick, colors=cols, linewidths=2.6, capstyle="round"))
            for a in heavy.GetAtoms():
                if a.GetSymbol() == "O":
                    ax.scatter(*xy[a.GetIdx()], s=9, color="#c0392b", zorder=3, edgecolors="w", linewidths=0.3)
            ax.set_xlim(xy[:, 0].min() - 0.5, xy[:, 0].max() + 0.5)
            ax.set_ylim(xy[:, 1].min() - 0.6, xy[:, 1].max() + 1.6)
            ax.set_aspect("equal")
            ax.axis("off")
            ax.text(0.5, 1.0, lab, transform=ax.transAxes, ha="center", va="top", fontsize=7)
        fig.text(0.255, 0.51, "→", fontsize=12, ha="center", va="center")
        fig.text(0.66, 0.50, "≈", fontsize=14, ha="center", va="center", fontweight="bold")
        fig.text(
            0.13, 0.12, "same correlation,\nsame search result", fontsize=6, ha="center", va="center", color="0.25"
        )
        fig.savefig(FIG / "toc_graphic.pdf", bbox_inches=None)
        fig.savefig(FIG / "toc_graphic.png", dpi=300, bbox_inches=None)
        fig.savefig(FIG / "toc_graphic.tif", dpi=300, bbox_inches=None, pil_kwargs={"compression": "tiff_lzw"})
        plt.close(fig)
    # flat RGB on white for the journal (matplotlib writes RGBA; ACS graphics should not carry transparency)
    from PIL import Image

    with Image.open(FIG / "toc_graphic.tif") as im:
        flat = Image.new("RGB", im.size, "white")
        flat.paste(im, mask=im.convert("RGBA").getchannel("A"))
    flat.save(FIG / "toc_graphic.tif", dpi=(300, 300), compression="tiff_lzw")
    print("wrote toc_graphic")


if __name__ == "__main__":
    which = sys.argv[1:]
    jobs = {
        "molecules": fig_molecules,
        "field": fig_field,
        "landscape": fig_landscape,
        "coverage": fig_coverage,
        "sigma": fig_sigma,
        "alanine": fig_alanine,
        "energy_hist": fig_energy_hist,
        "global_minima": fig_global_minima,
        "gallery": lambda: [fig_gallery(m) for m in ("ala3", "tetraglyme", "imatinib")],
        "ablation_forest": fig_ablation_forest,
        "known_survival": fig_known_survival,
        "isomers": fig_isomers,
        "hunt": fig_hunt,
        "bias2d": fig_bias2d,
        "cartesian": fig_cartesian,
        "dft": fig_dft,
        "lj38": fig_lj38,
        "concept": fig_concept,
        "toc": fig_toc,
        "dft_main": fig_dft_main,
    }
    for k, f in jobs.items():
        if not which or k in which:
            f()
