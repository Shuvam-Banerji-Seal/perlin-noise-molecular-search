"""LJ-cluster benchmark: time to reach the known global minimum, Perlin vs matched Gaussian field vs i.i.d.

Usage: uv run --extra plot python scripts/run_clusters.py [--jobs N] [--seeds 30] [--quick]
Writes results/clusters/{clusters.json, summary.md} and figures/clusters_{survival,structures}.{pdf,png}.
"""

import os

# One BLAS thread per process: with multi-threaded BLAS on a busy node, an
# LJ38 relaxation took 7.8 s instead of 67 ms (thread oversubscription).
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402

from perlin_mole.clusters import KNOWN_GM, basin_hopping, load_known, write_known_minima  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "clusters"
FIG = ROOT / "figures"
C_PERLIN, C_FIELD, C_IID = "#c0392b", "#2471a3", "#7f8c8d"

# budget (relaxations) and methods per cluster; L in sigma units. At L = 1 the Perlin
# correlation at the nearest-neighbour distance (~1.12 sigma) is ~ -0.07, i.e. no coherence,
# so L in {2, 4} is the informative range.
PLAN = {
    13: (100, ["gauss-iid", "perlin@2", "gauss-field@2", "perlin@4", "gauss-field@4"]),
    19: (400, ["gauss-iid", "perlin@2", "gauss-field@2", "perlin@4", "gauss-field@4"]),
    26: (600, ["gauss-iid", "perlin@2", "gauss-field@2", "perlin@4", "gauss-field@4"]),
    38: (3000, ["gauss-iid", "perlin@2", "gauss-field@2", "perlin@4", "gauss-field@4"]),
}
STEP = 0.3  # per-coordinate displacement std (sigma); gives 30-45% acceptance with i.i.d. moves
KT = 0.8


def _run(args):
    n, method, seed, budget = args
    return basin_hopping(n, method, seed, budget, s=STEP, kT=KT)


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return 0.0, 1.0
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, c - h), min(1.0, c + h)


def holm(ps):
    order = np.argsort(ps)
    m, adj, run = len(ps), [1.0] * len(ps), 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (m - rank) * ps[i]))
        adj[i] = run
    return adj


def analyse(runs: list[dict]) -> dict:
    res: dict = {"step_sigma": STEP, "kT": KT, "clusters": {}, "tests": []}
    for n in sorted({r["n"] for r in runs}):
        rr = [r for r in runs if r["n"] == n]
        budget = rr[0]["budget"]
        cl = {"budget": budget, "e_gm": KNOWN_GM[n]["energy"], "methods": {}}
        for m in dict.fromkeys(r["method"] for r in rr):
            mr = sorted([r for r in rr if r["method"] == m], key=lambda r: r["seed"])
            hits = [r["t_hit"] for r in mr]
            k = sum(h is not None for h in hits)
            t = np.array([budget + 1 if h is None else h for h in hits])
            cl["methods"][m] = {
                "n_runs": len(mr),
                "success": k / len(mr),
                "success_ci95": wilson(k, len(mr)),
                "rmst": float(np.minimum(t, budget).mean()),
                "median_t": float(np.median(t)) if k / len(mr) >= 0.5 else None,
                "mean_gap": float(np.mean([r["gap"] for r in mr])),
                "acceptance": float(np.mean([r["acceptance"] for r in mr])),
                "survival_t": list(range(0, budget + 1, max(1, budget // 100))),
                "survival": [float((t > s).mean()) for s in range(0, budget + 1, max(1, budget // 100))],
                "by_seed_t": {r["seed"]: (budget + 1 if r["t_hit"] is None else r["t_hit"]) for r in mr},
                "by_seed_gap": {r["seed"]: r["gap"] for r in mr},
                "wall_s_mean": float(np.mean([r["wall_s"] for r in mr])),
            }
        res["clusters"][n] = cl
        for m in cl["methods"]:
            if not m.startswith("perlin"):
                continue
            L = m.split("@")[1]
            for ref in (f"gauss-field@{L}", "gauss-iid"):
                if ref not in cl["methods"]:
                    continue
                for metric in ("t", "gap"):
                    a, b = cl["methods"][m][f"by_seed_{metric}"], cl["methods"][ref][f"by_seed_{metric}"]
                    d = np.array([a[s] - b[s] for s in sorted(set(a) & set(b))], float)
                    nz = d[np.abs(d) > 1e-9]
                    p = float(wilcoxon(nz).pvalue) if nz.size >= 5 else 1.0
                    res["tests"].append(
                        {
                            "cluster": n,
                            "variant": m,
                            "reference": ref,
                            "metric": metric,
                            "n": int(d.size),
                            "mean_diff": float(d.mean()),
                            "wins": int((d < -1e-9).sum()),
                            "losses": int((d > 1e-9).sum()),
                            "p": p,
                        }
                    )
    for metric in ("t", "gap"):
        fam = [t for t in res["tests"] if t["metric"] == metric]
        for t, ph in zip(fam, holm([t["p"] for t in fam]), strict=True):
            t["p_holm"] = ph
    return res


def markdown(res: dict) -> str:
    lines = [
        "# LJ clusters: time to the known global minimum",
        "",
        f"Basin hopping, kT = {res['kT']}, displacement std = {res['step_sigma']} sigma per coordinate "
        "for every method. Known GM energies: Cambridge Cluster Database. Hit = |E - E_GM| < 1e-4. "
        "Times right-censored at the budget.",
        "",
    ]
    for n, cl in res["clusters"].items():
        lines += [
            f"## LJ{n} (E_GM = {cl['e_gm']}, budget {cl['budget']} relaxations)",
            "",
            "| method | success [95% Wilson CI] | restricted mean steps | median steps | mean final gap | acceptance |",
            "|---|---|---|---|---|---|",
        ]
        for m, v in cl["methods"].items():
            lo, hi = v["success_ci95"]
            med = "-" if v["median_t"] is None else f"{v['median_t']:.0f}"
            lines.append(
                f"| {m} | {v['success']:.2f} [{lo:.2f}, {hi:.2f}] | {v['rmst']:.0f} | {med} | {v['mean_gap']:.3f} | "
                f"{v['acceptance']:.2f} |"
            )
        lines.append("")
    lines += [
        "## Paired tests (by seed; same start per seed)",
        "",
        "metric t = steps to GM (censored = budget+1); gap = best E - E_GM. Wins = variant better. Holm within metric.",
        "",
        "| cluster | variant − reference | metric | mean Δ | wins/losses | p | p_Holm |",
        "|---|---|---|---|---|---|---|",
    ]
    for t in res["tests"]:
        lines.append(
            f"| LJ{t['cluster']} | {t['variant']} − {t['reference']} | {t['metric']} | {t['mean_diff']:+.3f} | "
            f"{t['wins']}/{t['losses']} | {t['p']:.3g} | {t['p_holm']:.3g} |"
        )
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------------------- figures
def _style(m: str) -> dict:
    if m == "gauss-iid":
        return {"color": C_IID, "ls": ":", "lw": 1.4}
    L = m.split("@")[1]
    return {"color": C_PERLIN if m.startswith("perlin") else C_FIELD, "ls": "-" if L == "4" else "--", "lw": 1.6}


def fig_survival(res: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cls = list(res["clusters"].items())
    fig, axes = plt.subplots(1, len(cls), figsize=(7.2, 2.4), sharey=True)
    for ax, (n, cl) in zip(axes, cls, strict=True):
        for m, v in cl["methods"].items():
            lab = m.replace("@", " L=") + "σ" if "@" in m else m
            ax.step(v["survival_t"], v["survival"], where="post", label=lab, **_style(m))
        ax.set_title(f"LJ{n} ({KNOWN_GM[n]['point_group']} GM)", fontsize=8.5)
        ax.set_xlabel("relaxations", fontsize=8)
        ax.set_ylim(-0.02, 1.02)
        ax.tick_params(labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("P(known GM not yet found)", fontsize=8)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=5, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.12))
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"clusters_survival.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def _render(ax, x: np.ndarray, title: str) -> None:
    from matplotlib.patches import Circle

    x = x - x.mean(0)
    _, _, vt = np.linalg.svd(x, full_matrices=False)
    x = x @ vt.T
    z = x[:, 2]
    zn = (z - z.min()) / (np.ptp(z) + 1e-9)
    d = np.linalg.norm(x[:, None] - x[None], axis=-1)
    coord = ((d < 1.35) & (d > 0)).sum(1)
    cmap = {12: "#b03a2e", 9: "#d4ac0d", 8: "#28b463", 7: "#2e86c1", 6: "#7d3c98"}
    for i in np.argsort(z):
        c = cmap.get(int(coord[i]), "#5d6d7e")
        ax.add_patch(
            Circle(
                x[i, :2], 0.36 * (0.75 + 0.4 * zn[i]), facecolor=c, edgecolor="k", lw=0.4, alpha=0.95, zorder=2 + zn[i]
            )
        )
    r = np.abs(x[:, :2]).max() + 0.7
    ax.set_xlim(-r, r)
    ax.set_ylim(-r, r)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=7.5)


def fig_structures(runs: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = [n for n in (13, 38) if any(r["n"] == n for r in runs)]
    fig, axes = plt.subplots(len(rows), 3, figsize=(6.4, 2.3 * len(rows)), squeeze=False)
    for ax_row, n in zip(axes, rows, strict=True):
        known = load_known(n)
        e_gm = KNOWN_GM[n]["energy"]
        _render(ax_row[0], known, f"LJ{n} known GM ({KNOWN_GM[n]['point_group']})\nE = {e_gm:.4f}")
        for ax, fam in zip(ax_row[1:], ("perlin", "gauss-field"), strict=True):
            cand = [r for r in runs if r["n"] == n and r["method"].startswith(fam)]
            best = min(cand, key=lambda r: r["best_e"])
            ok = "= GM" if abs(best["best_e"] - e_gm) < 1e-4 else f"gap {best['best_e'] - e_gm:+.3f}"
            nh = sum(r["t_hit"] is not None for r in cand)
            _render(
                ax,
                np.array(best["best_x"]),
                f"best {fam} (L={best['method'].split('@')[1]}σ)\nE = {best['best_e']:.4f} ({ok}); "
                f"hits {nh}/{len(cand)} (L=2 and 4 pooled)",
            )
    fig.text(
        0.5,
        -0.01,
        "colour = coordination number (red 12, yellow 9, green 8, blue 7, purple 6)",
        ha="center",
        fontsize=7,
    )
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"clusters_structures.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--seeds", type=int, default=30)
    ap.add_argument("--quick", action="store_true", help="tiny budgets for a smoke test")
    ap.add_argument("--analyse-only", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    raw_path = OUT / "runs.json"
    if not a.analyse_only:
        write_known_minima()
        tasks = [
            (n, m, sd, max(10, b // 20) if a.quick else b)
            for n, (b, methods) in PLAN.items()
            for m in methods
            for sd in range(a.seeds)
        ]
        tasks.sort(key=lambda t: -t[0] * t[3])  # biggest first
        t0 = time.perf_counter()
        with ProcessPoolExecutor(max_workers=a.jobs) as pool:
            runs = list(pool.map(_run, tasks, chunksize=1))
        print(
            f"{len(runs)} runs in {time.perf_counter() - t0:.0f} s wall, "
            f"{sum(r['wall_s'] for r in runs) / 3600:.2f} CPU-h",
            flush=True,
        )
        raw_path.write_text(json.dumps(runs))
    runs = json.loads(raw_path.read_text())
    res = analyse(runs)
    (OUT / "clusters.json").write_text(json.dumps(res, indent=1))
    (OUT / "summary.md").write_text(markdown(res))
    fig_survival(res)
    fig_structures(runs)
    print(markdown(res))


if __name__ == "__main__":
    main()
