"""Re-score every grid-referenced analysis with the corrected ground truth.

Grid relaxations of clashed starts occasionally invert a stereocentre or flip a secondary amide to cis. Torsion
moves cannot do either, so such minima are outside the searched space (an inverted one is a different
stereoisomer). ``validate.reference_targets`` now excludes them; this script re-runs the analyses that score
against grids from the saved runs (no new searches):

* results/ablation                         five-molecule ablation (Ac-(Ala)3-NHMe: one cis-amide target removed)
* results/isomers/enantiomers[-guarded]    propranolol grids: 27 (S) and 2 (R) stereo-inverted targets removed
* results/hunt/kT0.3, kT3                  temperature probe
* results/cartesian                        literal Cartesian kicks
* results/lowT_confirm/filtered            sensitivity analysis only: the pre-registered result.json is kept

Usage: uv run python scripts/rescore_ground_truth.py
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402

from perlin_mole import ablation, benchmark  # noqa: E402
from perlin_mole.isomers import ENANTIOMERS, cip_labels, mirror_smiles  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
VALID = RES / "validation"
DEFAULT = list(ablation.COMPARISONS)


def analyze(raw_dir: Path, comparisons: list | None, grid_dir: Path = VALID, out_dir: Path | None = None) -> dict:
    ablation.COMPARISONS[:] = comparisons or DEFAULT
    try:
        rep = ablation.analyze(raw_dir, grid_dir, out_dir or raw_dir)
    finally:
        ablation.COMPARISONS[:] = DEFAULT
    excl = {m: v.get("n_excluded_fixed_features", 0) for m, v in rep["targets"].items()}
    print(f"{raw_dir.relative_to(ROOT)}: excluded targets {excl}", flush=True)
    return rep


def main() -> None:
    analyze(RES / "ablation", None)

    # enantiomer suites: molecule names resolve through benchmark.KNOWN, as in scripts/run_isomers.py
    mols = {}
    for n, s in ENANTIOMERS.items():
        for smi in (s, mirror_smiles(s)):
            mols[f"{n}-{''.join(cip_labels(smi).values())}"] = smi
    benchmark.KNOWN.update(mols)
    iso_path = RES / "isomers" / "isomers.json"
    iso = json.loads(iso_path.read_text())
    for key, sub in (("enantiomers", "enantiomers"), ("enantiomers_guarded", "enantiomers-guarded")):
        rep = analyze(RES / "isomers" / sub, None, RES / "isomers" / "grids")
        iso[key]["ablation_summary"] = rep["summary"]
        iso[key]["ablation_tests"] = rep["tests"]
        iso[key]["ground_truth_excluded"] = {m: v["n_excluded_fixed_features"] for m, v in rep["targets"].items()}
    iso_path.write_text(json.dumps(iso, indent=1))

    for kt in ("0.3", "3"):
        analyze(
            RES / "hunt" / f"kT{kt}",
            [
                ("perlin(L=4)", "gauss-corr(L=4)", f"Perlin vs matched (L=4), kT={kt}"),
                ("perlin(L=8)", "gauss-corr(L=8)", f"Perlin vs matched (L=8), kT={kt}"),
                ("gauss-corr(L=8)", "gauss-iid", f"correlation vs none, kT={kt}"),
            ],
        )

    analyze(
        RES / "cartesian",
        [
            ("cart-perlin@2", "cart-gauss-field@2", "Cartesian Perlin vs matched field (L=2 A)"),
            ("cart-perlin@4", "cart-gauss-field@4", "Cartesian Perlin vs matched field (L=4 A)"),
            ("cart-perlin@2", "cart-iid", "Cartesian Perlin (L=2 A) vs i.i.d."),
            ("cart-gauss-field@2", "cart-iid", "Cartesian correlated field vs i.i.d."),
        ],
    )

    # pre-registered low-T confirmation: the registered analysis (result.json) stays as it was; this is a
    # sensitivity analysis with the corrected targets
    low = RES / "lowT_confirm"
    out = low / "filtered"
    out.mkdir(exist_ok=True)
    rep = analyze(low, [("perlin(L=4)", "gauss-corr(L=4)", "low-T confirmation")], out_dir=out)
    rec = {(r["molecule"], r["seed"], r["label"]): r["recall_final"] for r in rep["rows"]}
    keys = sorted({(m, s) for m, s, _ in rec})
    d = np.array([rec[(m, s, "perlin(L=4)")] - rec[(m, s, "gauss-corr(L=4)")] for m, s in keys])
    res = {
        "note": "sensitivity analysis, stereo/amide-consistent targets; pre-registered result: ../result.json",
        "n_pairs": len(d),
        "mean_diff": float(d.mean()),
        "wins": int((d > 0).sum()),
        "losses": int((d < 0).sum()),
        "wilcoxon_p_two_sided": float(wilcoxon(d).pvalue),
        "bootstrap_ci95": next(t["ci95"] for t in rep["tests"] if t["metric"] == "recall_final"),
    }
    (low / "result_filtered.json").write_text(json.dumps(res, indent=1))
    print("low-T sensitivity:", res, flush=True)


if __name__ == "__main__":
    main()
