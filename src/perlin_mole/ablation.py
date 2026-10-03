"""Ablation and time-to-known-structure analysis against exhaustive grid references.

Question 1 (ablation): which ingredient of the Perlin proposal, if any, changes
search performance? Every Perlin variant is compared with (a) the full Perlin
proposal and (b) a Gaussian matched to that variant's own correlation matrix.

Question 2 (speed on known structures): how many MMFF relaxations does a method
need to find the *independently known* global minimum (from the exhaustive
grid), and how fast does it recover the known low-energy set? Runs that never
hit the target are right-censored at the budget, so the statistics are
survival-type (Kaplan-Meier, restricted mean time) rather than medians over
successful runs only.
"""

from __future__ import annotations

import json
import os
import pickle
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

from perlin_mole.benchmark import KNOWN, MOLECULES, _system
from perlin_mole.search import SymmetryRMSD
from perlin_mole.validate import reference_targets

E_TOL, RMSD_TOL = 0.1, 0.25

# Pairs that answer one ablation question each: (variant, reference, question).
COMPARISONS = [
    ("perlin(L=4)", "gauss-corr(L=4)", "Perlin vs matched Gaussian (baseline)"),
    ("gauss-corr(L=4)", "gauss-iid", "correlation vs none"),
    ("perlin(L=4)", "gauss-iid", "Perlin vs i.i.d."),
    ("perlin(L=4)", "reset", "Perlin vs MCMM reset"),
    ("perlin-gaussmarg(L=4)", "perlin(L=4)", "marginals: Gaussianize Perlin"),
    ("gauss-perlinmarg(L=4)", "gauss-corr(L=4)", "marginals: Perlin-shaped Gaussian"),
    ("perlin-oct1(L=4)", "perlin(L=4)", "octaves: 1 vs 2"),
    ("perlin-oct4(L=4)", "perlin(L=4)", "octaves: 4 vs 2"),
    ("perlin-oct1(L=4)", "gauss-corr-oct1(L=4)", "1 octave vs its matched Gaussian"),
    ("perlin-oct4(L=4)", "gauss-corr-oct4(L=4)", "4 octaves vs its matched Gaussian"),
    ("perlin-norot(L=4)", "perlin(L=4)", "random rotation off"),
    ("perlin-norot(L=4)", "gauss-corr-norot(L=4)", "no rotation vs its matched Gaussian"),
    ("perlin-shuffled(L=4)", "perlin(L=4)", "topology: shuffled sites"),
    ("gauss-corr-shuffled(L=4)", "gauss-corr(L=4)", "topology: shuffled (Gaussian)"),
    ("perlin-drift(L=4)", "perlin(L=4)", "temporal coherence on"),
    ("perlin-drift(L=4)", "gauss-ar1(L=4)", "temporal coherence vs matched AR(1) Gaussian"),
    ("gauss-ar1(L=4)", "gauss-corr(L=4)", "temporal coherence (Gaussian)"),
]

_CTX: dict = {}


def _init(targets: dict) -> None:
    _CTX.update(targets)


def _cov_mat(tx: np.ndarray, gx: np.ndarray, rmsd: SymmetryRMSD) -> dict:
    """COV/MAT ensemble metrics (as in the ML conformer-generation literature) between a reference set
    ``tx`` and a generated set ``gx``. R = recall direction (per reference), P = precision (per generated)."""
    if len(gx) == 0:
        nan = float("nan")
        return {"cov_r_05": 0.0, "cov_r_10": 0.0, "mat_r": nan, "cov_p_05": nan, "cov_p_10": nan, "mat_p": nan}
    d = np.stack([rmsd.many(x, gx) for x in tx])  # (n_ref, n_gen)
    r, p = d.min(axis=1), d.min(axis=0)
    return {
        "cov_r_05": float((r < 0.5).mean()),
        "cov_r_10": float((r < 1.0).mean()),
        "mat_r": float(r.mean()),
        "cov_p_05": float((p < 0.5).mean()),
        "cov_p_10": float((p < 1.0).mean()),
        "mat_p": float(p.mean()),
    }


def _score_run(r: dict) -> dict:
    """First step at which a run found each target minimum (None if never), plus COV/MAT."""
    mol = r["task"]["molecule"]
    t = _CTX[mol]
    te, tx, rmsd = t["e"], t["x"], t["rmsd"]
    first = [None] * len(te)
    mins = sorted(r["minima"], key=lambda m: m[2])  # by step
    me = np.array([m[0] for m in mins])
    mx = np.stack([m[1] for m in mins])
    ms = np.array([m[2] for m in mins])
    for k in range(len(te)):
        cand = np.flatnonzero(np.abs(me - te[k]) < E_TOL)
        if cand.size:
            hit = cand[rmsd.many(tx[k], mx[cand]) < RMSD_TOL]
            if hit.size:
                first[k] = int(ms[hit].min())
    gen = mx[me <= te[0] + t["window"]]  # generated ensemble: run minima inside the known energy window
    return {"molecule": mol, "label": r["label"], "seed": r["task"]["seed"], "first": first, **_cov_mat(tx, gen, rmsd)}


def kaplan_meier(times: list[int | None], horizon: int) -> np.ndarray:
    """P(not yet hit) at steps 0..horizon; None = censored at horizon."""
    n = len(times)
    t = np.array([horizon + 1 if x is None else x for x in times])
    steps = np.arange(horizon + 1)
    return np.array([(t > s).sum() / n for s in steps])  # no intermediate censoring: KM = empirical survival


def _nanmean(v: list[float]) -> float:
    v = [x for x in v if x == x]
    return float(np.mean(v)) if v else float("nan")


def _paired(a: dict, b: dict) -> tuple[np.ndarray, list]:
    keys = sorted(set(a) & set(b))
    return np.array([a[k] - b[k] for k in keys]), keys


def _boot_ci(d: np.ndarray, n: int = 4000, seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    means = d[rng.integers(0, len(d), size=(n, len(d)))].mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def _holm(ps: list[float]) -> list[float]:
    order = np.argsort(ps)
    m, adj, run = len(ps), [1.0] * len(ps), 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (m - rank) * ps[i]))
        adj[i] = run
    return adj


def analyze(
    raw_dir: str | Path = "results/ablation",
    grid_dir: str | Path = "results/validation",
    out_dir: str | Path | None = None,
    window: float = 3.0,
    jobs: int | None = None,
) -> dict:
    raw_dir, grid_dir = Path(raw_dir), Path(grid_dir)
    out = Path(out_dir) if out_dir else raw_dir
    with open(raw_dir / "raw.pkl", "rb") as fh:
        runs = pickle.load(fh)
    budget = runs[0]["task"]["n_steps"]
    smiles = {**MOLECULES, **KNOWN}

    targets, meta = {}, {}
    for mol in dict.fromkeys(r["task"]["molecule"] for r in runs):
        if not (grid_dir / f"{mol}_systematic.npz").exists():
            print(f"no grid reference for {mol}; skipped")
            continue
        g = np.load(grid_dir / f"{mol}_systematic.npz")
        low, e0, n_excl = reference_targets(g["energies"], g["heavy"], smiles[mol], window)
        targets[mol] = {
            "e": g["energies"][low],
            "x": g["heavy"][low],
            "rmsd": SymmetryRMSD(_system(smiles[mol]).heavy_mol()),
            "window": window,
        }
        meta[mol] = {"n_known_low": int(len(low)), "e_known_min": e0, "n_excluded_fixed_features": n_excl}

    runs = [r for r in runs if r["task"]["molecule"] in targets]
    jobs = jobs or os.cpu_count() or 1
    with ProcessPoolExecutor(max_workers=jobs, initializer=_init, initargs=(targets,)) as pool:
        scored = list(pool.map(_score_run, runs, chunksize=4))

    # per-run metrics
    rows = []
    for s in scored:
        f = s["first"]
        n = len(f)
        te = targets[s["molecule"]]["e"]
        # the global minimum can be degenerate (e.g. a mirror-image pair of an achiral molecule):
        # reaching any member counts
        glob = [f[k] for k in range(n) if te[k] - te[0] < 1e-3 and f[k] is not None]
        found = sorted(x for x in f if x is not None)
        t50 = found[int(np.ceil(0.5 * n)) - 1] if len(found) >= np.ceil(0.5 * n) else None
        rows.append(
            {
                **{k: s[k] for k in ("molecule", "label", "seed")},
                "t_global": min(glob) if glob else None,
                "n_degenerate_global": int(sum(te - te[0] < 1e-3)),
                "t50": t50,
                "recall": {c: sum(x <= c for x in found) / n for c in (50, 100, 250, budget)},
                "recall_final": len(found) / n,
                **{k: s[k] for k in ("cov_r_05", "cov_r_10", "mat_r", "cov_p_05", "cov_p_10", "mat_p")},
            }
        )

    labels = list(dict.fromkeys(r["label"] for r in rows))
    mols = list(targets)
    summary: dict = {}
    for mol in mols:
        summary[mol] = {}
        for lab in labels:
            rr = [r for r in rows if r["molecule"] == mol and r["label"] == lab]
            tg = [r["t_global"] for r in rr]
            cens = np.array([budget + 1 if t is None else t for t in tg])
            summary[mol][lab] = {
                "hit_rate": float(np.mean([t is not None for t in tg])),
                "rmst_global": float(np.minimum(cens, budget).mean()),  # restricted mean time-to-known-global-min
                "median_t_global": None if np.mean([t is not None for t in tg]) < 0.5 else float(np.median(cens)),
                "recall_final": float(np.mean([r["recall_final"] for r in rr])),
                "recall_100": float(np.mean([r["recall"][100] for r in rr])),
                **{k: _nanmean([r[k] for r in rr]) for k in ("cov_r_05", "cov_r_10", "mat_r", "cov_p_05", "mat_p")},
                "survival": kaplan_meier(tg, budget)[:: max(1, budget // 50)].tolist(),
            }

    # paired comparisons, pooled over (molecule, seed)
    def per(lab: str, key: str) -> dict:
        out = {}
        for r in rows:
            if r["label"] != lab:
                continue
            if key == "t_global":
                v = budget + 1 if r["t_global"] is None else r["t_global"]
            else:
                v = r["recall_final"]
            out[(r["molecule"], r["seed"])] = v
        return out

    tests = []
    for a, b, q in COMPARISONS:
        if a not in labels or b not in labels:
            continue
        for key in ("recall_final", "t_global"):
            d, keys = _paired(per(a, key), per(b, key))
            nz = d[np.abs(d) > 1e-12]
            p = float(wilcoxon(d).pvalue) if nz.size >= 5 else 1.0
            lo, hi = _boot_ci(d)
            tests.append(
                {
                    "variant": a,
                    "reference": b,
                    "question": q,
                    "metric": key,
                    "n": len(d),
                    "mean_diff": float(d.mean()),
                    "ci95": [lo, hi],
                    "wins": int((d > 1e-12).sum() if key == "recall_final" else (d < -1e-12).sum()),
                    "losses": int((d < -1e-12).sum() if key == "recall_final" else (d > 1e-12).sum()),
                    "p": p,
                }
            )
    for key in ("recall_final", "t_global"):
        fam = [t for t in tests if t["metric"] == key]
        for t, ph in zip(fam, _holm([t["p"] for t in fam]), strict=True):
            t["p_holm"] = ph

    report = {"budget": budget, "window": window, "targets": meta, "summary": summary, "tests": tests, "rows": rows}
    (out / "ablation.json").write_text(json.dumps(report, indent=1))
    (out / "ablation.md").write_text(format_markdown(report))
    return report


def format_markdown(rep: dict) -> str:
    b = rep["budget"]
    lines = [
        "# Ablation and time-to-known-structure analysis",
        "",
        f"Targets are the exhaustive-grid minima within {rep['window']} kcal/mol of the grid's global minimum. "
        f"Budget = {b} relaxations. Times are right-censored at the budget (counted as {b + 1}).",
        "",
        "## Paired comparisons (pooled over molecule × seed)",
        "",
        "Δrecall > 0 and Δt < 0 favour the *variant*. 95% CI: paired bootstrap. p: Wilcoxon; Holm within metric.",
        "",
        "| question | variant − reference | Δ final recall [95% CI] | p (Holm) "
        "| Δ steps to known global min [95% CI] | p (Holm) |",
        "|---|---|---|---|---|---|",
    ]
    by_q: dict = {}
    for t in rep["tests"]:
        by_q.setdefault((t["question"], t["variant"], t["reference"]), {})[t["metric"]] = t
    for (q, a, bref), m in by_q.items():
        r, tg = m["recall_final"], m["t_global"]
        lines.append(
            f"| {q} | {a} − {bref} | {r['mean_diff']:+.3f} [{r['ci95'][0]:+.3f}, {r['ci95'][1]:+.3f}] | "
            f"{r['p']:.3g} ({r['p_holm']:.2g}) | {tg['mean_diff']:+.1f} [{tg['ci95'][0]:+.1f}, {tg['ci95'][1]:+.1f}] | "
            f"{tg['p']:.3g} ({tg['p_holm']:.2g}) |"
        )
    lines.append("")
    for mol, s in rep["summary"].items():
        t = rep["targets"][mol]
        lines += [
            f"## {mol}: {t['n_known_low']} known low-energy minima",
            "",
            "| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 "
            "| COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
        for lab, v in sorted(s.items(), key=lambda kv: kv[1]["rmst_global"]):
            lines.append(
                f"| {lab} | {v['hit_rate']:.2f} | {v['rmst_global']:.0f} | "
                f"{v['recall_final']:.3f} | {v['recall_100']:.3f} | {v['cov_r_05']:.3f} | {v['mat_r']:.3f} | "
                f"{v['cov_p_05']:.3f} | {v['mat_p']:.3f} |"
            )
        lines.append("")
    return "\n".join(lines)


if __name__ == "__main__":  # pragma: no cover
    import sys

    print(format_markdown(analyze(*(sys.argv[1:2] or ["results/ablation"]))))
