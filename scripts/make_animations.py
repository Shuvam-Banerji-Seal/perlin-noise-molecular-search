"""Animated GIFs for the README (docs/animations/). Every frame comes from the package's own code and data.

1. pes_bias.gif           the alanine-dipeptide (phi, psi) MMFF94 surface under a drifting Perlin bias (left) and under
                          well-tempered metadynamics (right), with Metropolis walkers (formulation A1)
2. perlin_to_molecule.gif a Perlin field read along tetraglyme's torsions becomes a coherent kick; MMFF94 relaxes it
                          and Metropolis accepts or rejects (formulation B, basin hopping)
3. lj38.gif               basin hopping with a Perlin displacement field on LJ38 reaching the known global minimum
4. molecule_breaking.gif  Cartesian Perlin kicks tear bonds and can invert stereocentres of Ac-(Ala)3-NHMe after
                          MMFF94 relaxation; the same field read on the bond graph (torsion kicks) cannot

Usage: uv run --extra plot python scripts/make_animations.py [pes] [molecule] [lj38] [breaking]
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "animations"
OUT.mkdir(parents=True, exist_ok=True)
C_PERLIN, C_GCORR = "#c0392b", "#2471a3"
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.spines.top": False, "axes.spines.right": False})


def _save(anim: FuncAnimation, name: str, fps: int) -> None:
    anim.save(OUT / name, writer=PillowWriter(fps=fps), dpi=90)
    plt.close(anim._fig)
    print(f"wrote {name} ({(OUT / name).stat().st_size / 1e6:.1f} MB)", flush=True)


# ---------------------------------------------------------------- 1. PES + Perlin bias vs metadynamics
def pes(frames: int = 90, steps_per_frame: int = 400, walkers: int = 24) -> None:
    from perlin_mole.bias2d import KT_300, PerlinBias, Surface

    surf = Surface(np.load(ROOT / "results" / "bias2d" / "surface.npy"), 10.0)
    g = np.linspace(-180, 180, 121)
    P, Q = np.meshgrid(g, g, indexing="xy")
    E = surf(P, Q)
    bias = PerlinBias(4.0, cells=4, drift=1.0, seed=3)
    kT, step = KT_300, 10.0
    rng = np.random.default_rng(0)
    # metadynamics on a periodic 2-degree grid: well-tempered hills (height 0.2, width 15 deg, bias factor 8)
    mg = np.linspace(-180, 180, 181)[:-1]
    MP, MQ = np.meshgrid(mg, mg, indexing="xy")
    V = np.zeros_like(MP)
    h0, w, gamma = 0.2, 15.0, 8.0

    def vgrid(p, q):
        i = np.round((p + 180) / 2).astype(int) % 180
        j = np.round((q + 180) / 2).astype(int) % 180
        return V[j, i]

    walk = {k: np.stack([rng.normal(-80, 5, walkers), rng.normal(75, 5, walkers)], 1) for k in ("perlin", "metad")}
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.6))
    state = {"t": 0.0, "cross": {"perlin": 0, "metad": 0}, "side": {k: walk[k][:, 0] > 0 for k in walk}}

    def energy(kind, x, t):
        if kind == "perlin":
            return surf(x[:, 0], x[:, 1]) + bias(x[:, 0], x[:, 1], t)
        return surf(x[:, 0], x[:, 1]) + vgrid(x[:, 0], x[:, 1])

    def mc(kind, t):
        nonlocal V
        x = walk[kind]
        e = energy(kind, x, t)
        for s in range(steps_per_frame):
            xt = (x + rng.normal(0, step, x.shape) + 180) % 360 - 180
            et = energy(kind, xt, t)
            acc = rng.random(len(x)) < np.exp(-(et - e) / kT)
            x[acc], e[acc] = xt[acc], et[acc]
            if kind == "metad" and s % 200 == 0:  # hill pace as in the study
                for p, q in x:
                    dp = (MP - p + 180) % 360 - 180
                    dq = (MQ - q + 180) % 360 - 180
                    h = h0 * np.exp(-vgrid(np.array([p]), np.array([q]))[0] / (kT * (gamma - 1)))
                    V += h * np.exp(-(dp**2 + dq**2) / (2 * w * w))
                e = energy(kind, x, t)
        core_b = (x[:, 0] > 40) & (x[:, 0] < 110)
        core_a = (x[:, 0] > -170) & (x[:, 0] < -40)
        state["cross"][kind] += int(np.sum(core_b & ~state["side"][kind]) + np.sum(core_a & state["side"][kind]))
        state["side"][kind] = np.where(core_b, True, np.where(core_a, False, state["side"][kind]))

    def draw(_):
        t = state["t"]
        mc("perlin", t)
        mc("metad", t)
        for ax, kind in zip(axes, ("perlin", "metad"), strict=True):
            ax.clear()
            Z = E + (bias(P, Q, t) if kind == "perlin" else vgrid(P, Q))
            ax.contourf(P, Q, np.clip(Z - Z.min(), 0, 16), levels=24, cmap="viridis")
            ax.scatter(walk[kind][:, 0], walk[kind][:, 1], s=9, c="white", edgecolors="k", linewidths=0.4)
            ax.set_xlim(-180, 180)
            ax.set_ylim(-180, 180)
            ax.set_xticks([-180, 0, 180])
            ax.set_yticks([-180, 0, 180])
            ax.set_xlabel(r"$\phi$ (deg)")
            ax.set_ylabel(r"$\psi$ (deg)")
            title = (
                r"$E+\lambda N(\phi,\psi,t)$: drifting Perlin bias"
                if kind == "perlin"
                else r"$E+V(t)$: well-tempered metadynamics (shared bias)"
            )
            ax.set_title(f"{title}\nbarrier crossings so far: {state['cross'][kind]}", fontsize=8)
        fig.suptitle(
            "Formulation A1 on alanine dipeptide: a random Perlin bias wanders; metadynamics fills the basins",
            fontsize=8.5,
        )
        fig.tight_layout()
        state["t"] += 0.02

    _save(FuncAnimation(fig, draw, frames=frames), "pes_bias.gif", fps=10)


# ---------------------------------------------------------------- 2. Perlin field -> kick -> molecule
def _kabsch(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Rotate/translate a onto b (both (n, 3)); returns the moved a."""
    ac, bc = a - a.mean(0), b - b.mean(0)
    u, _, vt = np.linalg.svd(ac.T @ bc)
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1, 1, d]) @ vt
    return ac @ r + b.mean(0)


def molecule(frames: int = 70, seed: int = 4) -> None:
    from rdkit import Chem

    from perlin_mole.chem import TorsionSystem
    from perlin_mole.proposals import build_proposal

    smi = "COCCOCCOCCOCCOC"
    s = TorsionSystem.from_smiles(smi)
    prop = build_proposal("perlin", s.dof_coordinates(), 60.0, 4.0, seed=seed)
    rng = np.random.default_rng(seed)
    xyz, e, _ = s.minimize(s.initial_xyz)
    heavy = s.heavy_idx
    hm = Chem.Mol(s.heavy_mol())
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in hm.GetBonds()]
    elem = [a.GetSymbol() for a in hm.GetAtoms()]
    central = {frozenset((j, k)): t for t, (_, j, k, _) in enumerate(s.torsions)}
    hist, best = [e], e
    ref = xyz[heavy]
    _, _, pcs = np.linalg.svd(ref - ref.mean(0))  # view frame: principal axes of the start conformer
    frames_data = []
    for _ in range(frames):
        d = prop.sample(rng)
        trial = s.apply_delta(xyz, np.deg2rad(d))
        xt, et, _ = s.minimize(trial)
        accept = et <= e or rng.random() < np.exp(-(et - e) / 1.0)
        frames_data.append((d.copy(), _kabsch(trial[heavy], ref), _kabsch(xt[heavy], ref), et, accept))
        if accept:
            xyz, e = xt, et
        best = min(best, et)
        hist.append(best)
    e0 = min(hist)

    fig = plt.figure(figsize=(7.6, 4.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 2.8], width_ratios=[2.4, 1], hspace=0.35, wspace=0.22)
    ax_k, ax_m, ax_e = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[:, 1])
    cmap = plt.get_cmap("RdBu_r")
    norm = plt.Normalize(-150, 150)

    def draw(i):
        d, xk, xr, et, acc = frames_data[i]
        ax_k.clear()
        ax_k.bar(np.arange(1, len(d) + 1), d, color=cmap(norm(d)))
        ax_k.axhline(0, color="0.6", lw=0.6)
        ax_k.set_ylim(-170, 170)
        ax_k.set_xlabel("torsion along the chain")
        ax_k.set_ylabel(r"$\Delta\phi_j$ (deg)")
        ax_k.set_title("Perlin field read at the torsion sites = one coherent kick", fontsize=8)
        ax_m.clear()
        ang = np.deg2rad(3 * i)  # tumble about the long (first principal) axis
        rot = np.array([[1, 0, 0], [0, np.cos(ang), -np.sin(ang)], [0, np.sin(ang), np.cos(ang)]])
        for xx, alpha, lw in ((xk, 0.3, 1.5), (xr, 1.0, 3.2)):
            p = ((xx - xx.mean(0)) @ pcs.T) @ rot.T
            for a, b in bonds:
                t = central.get(frozenset((a, b)))
                col = cmap(norm(d[t])) if (t is not None and alpha == 1.0) else "0.4"
                ax_m.plot(p[[a, b], 0], p[[a, b], 1], color=col, lw=lw, alpha=alpha, solid_capstyle="round")
            if alpha == 1.0:
                for k, el in enumerate(elem):
                    ax_m.scatter(
                        p[k, 0],
                        p[k, 1],
                        s=46 if el == "O" else 16,
                        zorder=3,
                        edgecolors="w",
                        linewidths=0.6,
                        color="#c0392b" if el == "O" else "0.25",
                    )
        ax_m.set_aspect("equal")
        ax_m.axis("off")
        verdict = "accepted" if acc else "rejected"
        ax_m.set_title(f"kick (faint) -> MMFF94 relaxation (bold): E = {et - e0:+.2f} kcal/mol, {verdict}", fontsize=8)
        ax_e.clear()
        ax_e.plot(np.arange(i + 2), np.array(hist[: i + 2]) - e0, color=C_PERLIN)
        ax_e.set_xlim(0, frames)
        ax_e.set_ylim(-0.3, max(hist) - e0 + 0.5)
        ax_e.set_xlabel("basin-hopping step")
        ax_e.set_ylabel("lowest energy found (kcal/mol)")
        ax_e.set_title("tetraglyme, MMFF94", fontsize=8)

    _save(FuncAnimation(fig, draw, frames=len(frames_data)), "perlin_to_molecule.gif", fps=4)


# ---------------------------------------------------------------- 3. LJ38 reaching its known global minimum
def lj38(seed: int = 5, budget: int = 260, stride: int = 3) -> None:
    from perlin_mole.clusters import KNOWN_GM, containment_radius, lj_energy, make_proposal, random_start, relax

    n, method, s, kT = 38, "perlin@2", 0.3, 0.8
    x = random_start(n, np.random.default_rng(10_000 + seed))
    rng = np.random.default_rng([seed, sum(map(ord, method))])
    prop = make_proposal(method, s, seed)
    e_gm, rc = KNOWN_GM[n]["energy"], containment_radius(n)
    e = lj_energy(x.ravel())
    best, best_x, trace, snaps, hit = e, x.copy(), [e], [], None
    for step in range(1, budget + 1):
        xt, et = relax(x + prop.sample(x, rng))
        xt = xt - xt.mean(axis=0)
        if np.max(np.linalg.norm(xt, axis=1)) <= rc:
            if et < best:
                best, best_x = et, xt.copy()
            if hit is None and et < e_gm + 1e-4:
                hit = step
            if et <= e or rng.random() < np.exp(-(et - e) / kT):
                x, e = xt, et
        trace.append(best)
        if step % stride == 0 or step == hit:
            snaps.append((step, best_x.copy(), best))
    print(f"LJ38 seed {seed}: global minimum at step {hit}", flush=True)

    fig = plt.figure(figsize=(7.2, 3.6))
    ax3 = fig.add_subplot(1, 2, 1, projection="3d")
    ax_e = fig.add_subplot(1, 2, 2)

    def coord(xx):
        dmat = np.linalg.norm(xx[:, None] - xx[None], axis=-1)
        return ((dmat > 0) & (dmat < 1.35)).sum(1)

    def draw(i):
        step, xx, eb = snaps[i]
        ax3.clear()
        c = coord(xx)
        ax3.scatter(*xx.T, c=c, cmap="plasma", vmin=4, vmax=12, s=420, edgecolors="k", linewidths=0.4, depthshade=True)
        lim = np.abs(xx).max() * 0.85
        ax3.set_xlim(-lim, lim)
        ax3.set_ylim(-lim, lim)
        ax3.set_zlim(-lim, lim)
        ax3.view_init(elev=20, azim=2 * i)
        ax3.set_axis_off()
        lab = "known global minimum" if (hit and step >= hit) else "lowest minimum so far"
        ax3.set_title(f"LJ$_{{38}}$ step {step}: {lab}", fontsize=8)
        ax_e.clear()
        ax_e.plot(np.arange(step + 1), trace[: step + 1], color=C_PERLIN, label="Perlin field ($L=2\\sigma$)")
        ax_e.axhline(e_gm, color="k", ls="--", lw=0.8, label=f"known global minimum ({e_gm:.3f})")
        ax_e.set_xlim(0, budget)
        ax_e.set_ylim(e_gm - 0.5, trace[min(10, len(trace) - 1)] + 1)
        ax_e.set_xlabel("relaxations")
        ax_e.set_ylabel(r"lowest energy ($\epsilon$)")
        ax_e.legend(frameon=False, fontsize=7, loc="upper right")
        ax_e.set_title("atom colour: coordination number", fontsize=7)
        fig.subplots_adjust(left=0.02, right=0.98, wspace=0.25, bottom=0.15)

    _save(FuncAnimation(fig, draw, frames=len(snaps)), "lj38.gif", fps=8)


# ---------------------------------------------------------------- 4. Cartesian kicks tear molecules; torsion kicks cannot
def _measured_break_rate(mol: str, method: str) -> float | None:
    """Fraction of relaxations that changed connectivity or stereochemistry in the Cartesian suite (raw.pkl)."""
    import pickle

    p = ROOT / "results" / "cartesian" / "raw.pkl"
    if not p.exists():
        return None
    runs = [r for r in pickle.load(open(p, "rb")) if r["task"]["molecule"] == mol and r["task"]["method"] == method]
    return float(np.mean([r["broken"] / r["task"]["n_steps"] for r in runs])) if runs else None


def breaking(seed: int = 4, kicks: int = 6, ramp: int = 7, hold: int = 3, relax: int = 6) -> None:
    """Ac-(Ala)3-NHMe: independent kicks from the same MMFF94 minimum, with the study's settings.

    Left: the literal idea, a Perlin field displacing every atom (L = 2 A, 0.5 A per coordinate). Bonds stretch
    and tear; after MMFF94 relaxation the molecule is checked for the same connectivity and CIP labels.
    Right: formulation B, a Perlin field read on the bond graph and applied as torsion changes
    (sigma = 60 deg, L = 4): bond lengths and stereocentres are untouched by construction.
    The seed is chosen so that the clip contains a chemistry change; the measured rate is printed on the frame.
    """
    from perlin_mole.benchmark import MOLECULES
    from perlin_mole.cartesian import _heavy_graph, intact
    from perlin_mole.chem import TorsionSystem
    from perlin_mole.clusters import make_proposal
    from perlin_mole.isomers import cip_labels
    from perlin_mole.proposals import build_proposal

    smi = MOLECULES["ala3"]
    s = TorsionSystem.from_smiles(smi, seed=0)
    x0, _, _ = s.minimize(s.initial_xyz)
    ref_graph, ref_cip = _heavy_graph(s.mol), cip_labels(smi)
    heavy = list(s.heavy_idx)
    stereo = {heavy[i]: lab for i, lab in ref_cip.items()}  # all-atom index -> CIP label of the input
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in s.mol.GetBonds()]
    r0 = np.array([np.linalg.norm(x0[a] - x0[b]) for a, b in bonds])
    elem = [a.GetSymbol() for a in s.mol.GetAtoms()]
    colour = {"C": "0.25", "N": "#2471a3", "O": "#c0392b", "H": "0.75"}
    _, _, pcs = np.linalg.svd(x0[heavy] - x0[heavy].mean(0))

    rng = np.random.default_rng(seed)
    cart = make_proposal("perlin@2", 0.5, seed)
    tors = build_proposal("perlin", s.dof_coordinates(), 60.0, 4.0, seed=seed)
    data = []
    for _ in range(kicks):
        d_c = cart.sample(x0, rng)
        d_t = np.deg2rad(tors.sample(rng))
        rel_c = s.minimize(x0 + d_c)[0]
        rel_t = s.minimize(s.apply_delta(x0, d_t))[0]
        verdict = []
        for rel in (rel_c, rel_t):
            if not intact(s, rel, ref_graph, {}):
                verdict.append(("connectivity changed", None))
            else:
                cip = cip_labels(smi, rel[heavy])
                verdict.append(("intact" if cip == ref_cip else "stereocentre inverted", cip))
        data.append((d_c, d_t, rel_c, rel_t, verdict))
    rate = _measured_break_rate("ala3", "cart-perlin@2")

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 4.0))
    per_kick = ramp + hold + relax

    def project(xx):
        return (xx - xx[heavy].mean(0)) @ pcs.T

    def draw_mol(ax, xx, faint=None, labels=None):
        p = project(xx)
        if faint is not None:
            q = project(faint)
            for a, b in bonds:
                ax.plot(q[[a, b], 0], q[[a, b], 1], color="0.8", lw=1.0, zorder=1)
        stretch = np.array([np.linalg.norm(xx[a] - xx[b]) for a, b in bonds]) / r0 - 1
        torn = 0
        for (a, b), st in zip(bonds, stretch, strict=True):
            if st > 0.4:
                torn += 1
                ax.plot(p[[a, b], 0], p[[a, b], 1], color="#c0392b", lw=1.0, ls=":", zorder=2)
            else:
                col = "#e67e22" if st > 0.15 else "0.35"
                ax.plot(p[[a, b], 0], p[[a, b], 1], color=col, lw=2.2, solid_capstyle="round", zorder=2)
        for k, el in enumerate(elem):
            ax.scatter(p[k, 0], p[k, 1], s=8 if el == "H" else 30, color=colour.get(el, "0.25"), zorder=3,
                       edgecolors="w", linewidths=0.4)
        if labels is not None:
            for k, lab in labels.items():
                ok = lab == stereo[k]
                ax.text(p[k, 0] + 0.25, p[k, 1] + 0.25, lab, fontsize=8, fontweight="bold", zorder=4,
                        color="#1e8449" if ok else "#c0392b")
        ax.set_xlim(-7.5, 7.5)
        ax.set_ylim(-4.6, 4.6)
        ax.set_aspect("equal")
        ax.axis("off")
        return torn

    def draw(f):
        kick, ph = divmod(f, per_kick)
        d_c, d_t, rel_c, rel_t, verdict = data[kick]
        frac = min(1.0, (ph + 1) / ramp)
        start_labels = dict(stereo)
        for side, ax in zip(("cart", "tors"), axes, strict=True):
            ax.clear()
            kicked = x0 + frac * d_c if side == "cart" else s.apply_delta(x0, frac * d_t)
            if ph < ramp + hold:
                torn = draw_mol(ax, kicked, labels=start_labels if ph == 0 else None)
                status = f"kick {kick + 1}: {100 * frac:.0f}% applied, bonds torn (>40% stretch): {torn}"
            else:
                rel = rel_c if side == "cart" else rel_t
                text, cip = verdict[0 if side == "cart" else 1]
                labels = {heavy[i]: lab for i, lab in (cip or {}).items()}
                draw_mol(ax, rel, faint=kicked, labels=labels)
                status = f"kick {kick + 1} after MMFF94 relaxation: {text}"
            head = (
                "Perlin field on Cartesian atom positions\n(the literal idea: L = 2 Å, 0.5 Å per coordinate)"
                if side == "cart"
                else "Perlin field on the bond graph (formulation B)\napplied as torsion changes (σ = 60°, L = 4)"
            )
            j = 0 if side == "cart" else 1
            done = kick + 1 if ph >= ramp + hold else kick  # kicks whose relaxation has been shown
            changed = sum(dd[4][j][0] != "intact" for dd in data[:done])
            ax.set_title(f"{head}\n{status}\nchemistry changed so far: {changed} of {done} kicks", fontsize=7.5)
        note = f"; in the study {100 * rate:.1f}% of Cartesian Perlin relaxations of this molecule" if rate else ""
        fig.suptitle(
            "Ac-(Ala)$_3$-NHMe: Cartesian kicks tear bonds and can invert stereocentres, torsion kicks cannot"
            f"\nCIP labels: green = as in the input, red = inverted{note}"
            + (" changed its chemistry" if rate else ""),
            fontsize=7.5,
        )
        fig.subplots_adjust(left=0.01, right=0.99, top=0.76, bottom=0.01, wspace=0.04)

    _save(FuncAnimation(fig, draw, frames=kicks * per_kick), "molecule_breaking.gif", fps=6)


if __name__ == "__main__":
    which = sys.argv[1:] or ["pes", "molecule", "lj38", "breaking"]
    if "pes" in which:
        pes()
    if "molecule" in which:
        molecule()
    if "lj38" in which:
        lj38()
    if "breaking" in which:
        breaking()
