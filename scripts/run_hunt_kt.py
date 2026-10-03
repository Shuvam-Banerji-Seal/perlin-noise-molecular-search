"""Temperature sweep: does proposal quality (Perlin vs matched Gaussian) matter more when basin hopping is greedy?

Usage: uv run python scripts/run_hunt_kt.py   (all cores; writes results/hunt/kT*/)
"""

from pathlib import Path

from perlin_mole import ablation
from perlin_mole.benchmark import run_benchmark

ROOT = Path(__file__).resolve().parents[1]
MOLS = ["decane", "tetraglyme", "ala3"]
for kT in (0.3, 3.0):
    out = ROOT / "results" / "hunt" / f"kT{kT:g}"
    run_benchmark(MOLS, ["gauss-iid", "gauss-corr", "perlin"], 10, 500, 60.0, (4.0, 8.0), kT, 3.0, None, out)
    ablation.COMPARISONS[:] = [
        ("perlin(L=4)", "gauss-corr(L=4)", f"Perlin vs matched (L=4), kT={kT:g}"),
        ("perlin(L=8)", "gauss-corr(L=8)", f"Perlin vs matched (L=8), kT={kT:g}"),
        ("gauss-corr(L=8)", "gauss-iid", f"correlation vs none, kT={kT:g}"),
    ]
    rep = ablation.analyze(out, ROOT / "results" / "validation", out)
    print(ablation.format_markdown(rep)[:2500], flush=True)
