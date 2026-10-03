"""Equal-budget comparison of conformer-search proposal schemes.

Budget is counted in MMFF local minimisations (the expensive step that a
QM/ML judge would replace). Every method starts each seed from the same
relaxed ETKDG geometry.

Metrics per run, relative to the union of all minima found by all runs:

* ``gap``       best energy found minus the global best (kcal/mol)
* ``coverage``  fraction of the reference low-energy set (within ``window``
                of the global best) that this run found
* ``n_low``     unique minima found within the window
* ``hit_step``  first step reaching the global best (within 0.1 kcal/mol)
"""

from __future__ import annotations

import json
import os
import pickle
import zlib
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

from perlin_mole.chem import TorsionSystem
from perlin_mole.proposals import build_proposal
from perlin_mole.search import ConformerArchive, SymmetryRMSD, basin_hopping, deduplicate, random_embedding

MOLECULES = {
    "decane": "CCCCCCCCCC",
    "tetraglyme": "COCCOCCOCCOCCOC",
    "ala3": "CC(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)NC",
    "imatinib": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
}
# Small molecules with complete, independently known minima (used for time-to-known-structure tests).
KNOWN = {
    "ala1": "CC(=O)N[C@@H](C)C(=O)NC",  # alanine dipeptide, Ac-Ala-NHMe
}

METHODS = ("etkdg", "reset", "gauss-iid", "gauss-corr", "perlin")
UNCORRELATED = ("etkdg", "reset", "gauss-iid")


@dataclass(frozen=True)
class Task:
    molecule: str
    smiles: str
    method: str
    seed: int
    n_steps: int
    sigma_deg: float = 60.0
    corr_length: float = 4.0
    kT: float = 1.0
    stereo_guard: bool = False  # reject relaxed minima whose CIP labels differ from the input

    @property
    def label(self) -> str:
        if self.method in UNCORRELATED:
            return self.method
        return f"{self.method}(L={self.corr_length:g})"


@lru_cache(maxsize=8)
def _system(smiles: str) -> TorsionSystem:
    return TorsionSystem.from_smiles(smiles, seed=0)


def _start(system: TorsionSystem, seed: int) -> np.ndarray:
    return system.embed_random(10_000 + seed)


def run_task(task: Task) -> dict:
    system = _system(task.smiles)
    archive = ConformerArchive(SymmetryRMSD(system.heavy_mol()))
    start = _start(system, task.seed)
    if task.method == "etkdg":
        res = random_embedding(system, task.n_steps, task.seed, start, archive)
    else:
        proposal = build_proposal(
            task.method,
            system.dof_coordinates(3),
            np.deg2rad(task.sigma_deg),
            task.corr_length,
            seed=task.seed,
            n_steps=task.n_steps,
        )
        rng = np.random.default_rng([task.seed, zlib.crc32(task.method.encode())])
        guard = None
        if task.stereo_guard:
            from perlin_mole.isomers import cip_labels

            ref = cip_labels(task.smiles)
            if ref:
                guard = lambda xyz: cip_labels(task.smiles, xyz[system.heavy_idx]) == ref  # noqa: E731
        res = basin_hopping(system, proposal, task.n_steps, task.kT, rng, start, archive, guard=guard)
    return {
        "task": asdict(task),
        "label": task.label,
        "best_trace": res.best_trace.tolist(),
        "minima": [(m.energy, m.heavy_xyz, m.step) for m in res.minima],
        "acceptance": res.acceptance,
        "unconverged": res.unconverged,
        "wall_time": res.wall_time,
        "rejected_chem": res.rejected_chem,
    }


def _score_molecule(runs: list[dict], window: float, checkpoints: list[int], jobs: int = 1) -> tuple[list[dict], dict]:
    smiles = runs[0]["task"]["smiles"]
    rmsd = SymmetryRMSD(_system(smiles).heavy_mol())

    energies = np.array([e for r in runs for e, _, _ in r["minima"]])
    xyz = np.stack([x for r in runs for _, x, _ in r["minima"]])
    run_of = np.repeat(np.arange(len(runs)), [len(r["minima"]) for r in runs])
    step_of = np.array([s for r in runs for _, _, s in r["minima"]])
    order = np.argsort(energies, kind="stable")
    energies, xyz, run_of, step_of = energies[order], xyz[order], run_of[order], step_of[order]

    rep = deduplicate(energies, xyz, rmsd, jobs=jobs)
    unique = rep == np.arange(len(rep))
    e_star = float(energies[0])
    low = np.flatnonzero(unique & (energies <= e_star + window))
    low_set = set(low.tolist())

    # first step at which each run reached each low-energy reference minimum
    first: list[dict[int, int]] = [{} for _ in runs]
    for i in range(len(rep)):
        k = int(rep[i])
        if k in low_set:
            d = first[run_of[i]]
            d[k] = min(d.get(k, step_of[i]), step_of[i])

    scored = []
    for r, found_d in zip(runs, first, strict=True):
        found = list(found_d.values())
        best = np.asarray(r["best_trace"])
        hit = np.flatnonzero(best <= e_star + 0.1)
        scored.append(
            {
                "label": r["label"],
                "seed": r["task"]["seed"],
                "gap": max(float(best[-1] - e_star), 0.0) + 0.0,
                "coverage": len(found) / len(low),
                "coverage_curve": [sum(s <= c for s in found) / len(low) for c in checkpoints],
                "n_unique": len(r["minima"]),
                "n_low": int(sum(e <= e_star + window for e, _, _ in r["minima"])),
                "hit_step": int(hit[0]) if hit.size else None,
                "acceptance": r["acceptance"],
                "unconverged": r["unconverged"],
                "wall_time": r["wall_time"],
            }
        )
    meta = {"e_star": e_star, "n_reference_low": int(len(low)), "n_reference_total": int(unique.sum())}
    return scored, meta


def _summarise(scored: list[dict]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for label in dict.fromkeys(s["label"] for s in scored):
        rows = [s for s in scored if s["label"] == label]
        cov = np.array([s["coverage"] for s in rows])
        gap = np.array([s["gap"] for s in rows])
        hits = [s["hit_step"] for s in rows]
        out[label] = {
            "coverage_mean": float(cov.mean()),
            "coverage_std": float(cov.std(ddof=1)) if len(cov) > 1 else 0.0,
            "gap_mean": float(gap.mean()),
            "gap_median": float(np.median(gap)),
            "n_low_mean": float(np.mean([s["n_low"] for s in rows])),
            "hit_rate": float(np.mean([h is not None for h in hits])),
            "hit_step_median": float(np.median([h for h in hits if h is not None]))
            if any(h is not None for h in hits)
            else None,
            "acceptance": float(np.nanmean([s["acceptance"] for s in rows])) if label != "etkdg" else None,
            "wall_time_mean": float(np.mean([s["wall_time"] for s in rows])),
            "coverage_by_seed": {s["seed"]: s["coverage"] for s in rows},
        }
    return out


def _paired_test(summary: dict[str, dict], a: str, b: str) -> dict | None:
    """Wilcoxon signed-rank on per-seed coverage (same starting geometry per seed)."""
    if a not in summary or b not in summary:
        return None
    sa, sb = summary[a]["coverage_by_seed"], summary[b]["coverage_by_seed"]
    seeds = sorted(set(sa) & set(sb))
    diff = np.array([sa[s] - sb[s] for s in seeds])
    if len(seeds) < 5 or np.allclose(diff, 0):
        return {"a": a, "b": b, "mean_diff": float(diff.mean()) if diff.size else 0.0, "p": None, "n": len(seeds)}
    return {"a": a, "b": b, "mean_diff": float(diff.mean()), "p": float(wilcoxon(diff).pvalue), "n": len(seeds)}


def run_benchmark(
    molecules: list[str],
    methods: list[str],
    seeds: int,
    n_steps: int,
    sigma_deg: float = 60.0,
    corr_lengths: tuple[float, ...] = (4.0,),
    kT: float = 1.0,
    window: float = 3.0,
    jobs: int | None = None,
    out_dir: str | Path = "results/latest",
    stereo_guard: bool = False,
) -> dict:
    tasks = []
    for mol in molecules:
        smiles = {**MOLECULES, **KNOWN}.get(mol, mol)
        for method in methods:
            lengths = (corr_lengths[0],) if method in UNCORRELATED else corr_lengths
            for L in lengths:
                for seed in range(seeds):
                    tasks.append(Task(mol, smiles, method, seed, n_steps, sigma_deg, L, kT, stereo_guard))

    jobs = jobs or os.cpu_count() or 1
    # Longest molecules first keeps the pool busy until the end.
    order = sorted(range(len(tasks)), key=lambda i: -len(tasks[i].smiles))
    results: list[dict | None] = [None] * len(tasks)
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        for i, res in zip(order, pool.map(run_task, [tasks[i] for i in order], chunksize=1), strict=True):
            results[i] = res

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    # Persist raw runs first so scoring can be redone with rescore() without re-minimising.
    with open(out / "raw.pkl", "wb") as fh:
        pickle.dump(results, fh)
    return score_results(results, molecules, n_steps, seeds, sigma_deg, corr_lengths, kT, window, out, jobs)


def rescore(out_dir: str | Path, window: float = 3.0, jobs: int | None = None) -> dict:
    """Re-run scoring on a saved ``raw.pkl`` without repeating any minimisation."""
    out = Path(out_dir)
    with open(out / "raw.pkl", "rb") as fh:
        results = pickle.load(fh)
    t = results[0]["task"]
    molecules = list(dict.fromkeys(r["task"]["molecule"] for r in results))
    seeds = len({r["task"]["seed"] for r in results})
    lengths = tuple(sorted({r["task"]["corr_length"] for r in results if r["task"]["method"] == "perlin"})) or (
        t["corr_length"],
    )
    jobs = jobs or os.cpu_count() or 1
    return score_results(results, molecules, t["n_steps"], seeds, t["sigma_deg"], lengths, t["kT"], window, out, jobs)


def score_results(
    results: list[dict],
    molecules: list[str],
    n_steps: int,
    seeds: int,
    sigma_deg: float,
    corr_lengths: tuple[float, ...],
    kT: float,
    window: float,
    out: Path,
    jobs: int = 1,
) -> dict:
    checkpoints = sorted({max(1, n_steps * f // 10) for f in range(1, 11)})
    report: dict = {
        "config": {
            "seeds": seeds,
            "n_steps": n_steps,
            "sigma_deg": sigma_deg,
            "corr_lengths": list(corr_lengths),
            "kT": kT,
            "window": window,
            "checkpoints": checkpoints,
        },
        "molecules": {},
    }
    for mol in molecules:
        runs = [r for r in results if r["task"]["molecule"] == mol]
        system = _system(runs[0]["task"]["smiles"])
        scored, meta = _score_molecule(runs, window, checkpoints, jobs)
        summary = _summarise(scored)
        tests = [_paired_test(summary, f"perlin(L={L:g})", f"gauss-corr(L={L:g})") for L in corr_lengths] + [
            _paired_test(summary, f"perlin(L={corr_lengths[0]:g})", "gauss-iid")
        ]
        curves = {
            label: np.mean([s["coverage_curve"] for s in scored if s["label"] == label], axis=0).tolist()
            for label in summary
        }
        report["molecules"][mol] = {
            "smiles": system.smiles,
            "n_dof": system.n_dof,
            "n_atoms": system.n_atoms,
            **meta,
            "summary": summary,
            "tests": [t for t in tests if t],
            "coverage_curves": curves,
            "runs": scored,
        }

    report["pooled"] = pooled_tests(report)
    write_report(report, out)
    return report


def write_report(report: dict, out: Path) -> None:
    (out / "report.json").write_text(json.dumps(report, indent=1))
    (out / "summary.md").write_text(format_markdown(report))
    _try_plot(report, out / "coverage.png")


def pooled_tests(report: dict) -> list[dict]:
    """Wilcoxon over all (molecule, seed) pairs: is a method better *across* molecules?"""

    def by_seed(mol: dict, label: str) -> dict[str, float]:
        return {str(k): v for k, v in mol["summary"].get(label, {}).get("coverage_by_seed", {}).items()}

    pairs = []
    for L in report["config"]["corr_lengths"]:
        pairs += [
            (f"perlin(L={L:g})", f"gauss-corr(L={L:g})"),
            (f"gauss-corr(L={L:g})", "gauss-iid"),
            (f"perlin(L={L:g})", "gauss-iid"),
        ]
    pairs += [("gauss-iid", "reset"), ("etkdg", "gauss-iid")]
    out = []
    for a, b in pairs:
        diff = []
        for mol in report["molecules"].values():
            sa, sb = by_seed(mol, a), by_seed(mol, b)
            diff += [sa[k] - sb[k] for k in sorted(set(sa) & set(sb))]
        diff = np.array(diff)
        if diff.size < 5:
            continue
        p = None if np.allclose(diff, 0) else float(wilcoxon(diff).pvalue)
        out.append(
            {
                "a": a,
                "b": b,
                "n": int(diff.size),
                "mean_diff": float(diff.mean()),
                "wins": int((diff > 1e-12).sum()),
                "losses": int((diff < -1e-12).sum()),
                "p": p,
            }
        )
    return out


def format_markdown(report: dict) -> str:
    cfg = report["config"]
    lines = [
        "# Benchmark summary",
        "",
        f"seeds={cfg['seeds']}, steps={cfg['n_steps']} minimisations/run, sigma={cfg['sigma_deg']} deg, "
        f"kT={cfg['kT']} kcal/mol, window={cfg['window']} kcal/mol",
        "",
    ]
    for mol, r in report["molecules"].items():
        lines += [
            f"## {mol}  ({r['n_dof']} torsions, {r['n_atoms']} atoms)",
            "",
            f"Reference: {r['n_reference_total']} unique minima, {r['n_reference_low']} within window, "
            f"E* = {r['e_star']:.3f} kcal/mol",
            "",
            "| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |",
            "|---|---|---|---|---|---|---|",
        ]
        for label, s in sorted(r["summary"].items(), key=lambda kv: -kv[1]["coverage_mean"]):
            hs = "-" if s["hit_step_median"] is None else f"{s['hit_step_median']:.0f}"
            acc = "-" if s["acceptance"] is None else f"{s['acceptance']:.2f}"
            gap = max(s["gap_median"], 0.0) + 0.0
            lines.append(
                f"| {label} | {s['coverage_mean']:.3f} ± {s['coverage_std']:.3f} | {gap:.3f} | "
                f"{s['n_low_mean']:.1f} | {s['hit_rate']:.2f} | {hs} | {acc} |"
            )
        lines.append("")
        for t in r["tests"]:
            p = "n/a" if t["p"] is None else f"{t['p']:.3g}"
            lines.append(f"- {t['a']} − {t['b']}: Δcoverage = {t['mean_diff']:+.3f} (Wilcoxon p = {p}, n = {t['n']})")
        lines.append("")
    if report.get("pooled"):
        lines += [
            "## Pooled over molecules (paired by molecule and seed)",
            "",
            "| comparison | Δcoverage | wins / losses | Wilcoxon p | n |",
            "|---|---|---|---|---|",
        ]
        for t in report["pooled"]:
            p = "n/a" if t["p"] is None else f"{t['p']:.3g}"
            lines.append(
                f"| {t['a']} − {t['b']} | {t['mean_diff']:+.3f} | {t['wins']} / {t['losses']} | {p} | {t['n']} |"
            )
        lines.append("")
    return "\n".join(lines)


def _style(label: str, lengths: list[float]) -> dict:
    baseline = {"etkdg": "0.15", "reset": "0.45", "gauss-iid": "0.7"}
    if label in baseline:
        return {"color": baseline[label], "lw": 1.3, "ls": ":" if label == "gauss-iid" else "-"}
    L = float(label.split("L=")[1].rstrip(")"))
    palette = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
    color = palette[sorted(lengths).index(L) % len(palette)] if L in lengths else "k"
    return {"color": color, "lw": 2.0, "ls": "-" if label.startswith("perlin") else "--"}


def _try_plot(report: dict, path: Path) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    mols = list(report["molecules"])
    cps = report["config"]["checkpoints"]
    lengths = [float(x) for x in report["config"]["corr_lengths"]]
    fig, axes = plt.subplots(1, len(mols), figsize=(4.2 * len(mols), 3.9), squeeze=False)
    for ax, mol in zip(axes[0], mols, strict=True):
        for label, curve in report["molecules"][mol]["coverage_curves"].items():
            ax.plot(cps, curve, label=label, **_style(label, lengths))
        ax.set_title(f"{mol} ({report['molecules'][mol]['n_dof']} torsions)")
        ax.set_xlabel("MMFF minimisations")
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.3)
    axes[0][0].set_ylabel("low-energy coverage")
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=min(len(labels), 9), fontsize=8, frameon=False)
    fig.suptitle(
        f"solid = Perlin, dashed = Gaussian with the same correlation (colour = L); "
        f"sigma = {report['config']['sigma_deg']:g} deg",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    fig.savefig(path, dpi=140)
    plt.close(fig)
