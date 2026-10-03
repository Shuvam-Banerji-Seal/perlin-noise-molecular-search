"""Command-line entry point: ``perlin-mole {bench,rescore,report,landscape,molecules}``."""

from __future__ import annotations

import argparse

import numpy as np

from perlin_mole.benchmark import METHODS, MOLECULES, run_benchmark


def _bench(args: argparse.Namespace) -> None:
    report = run_benchmark(
        molecules=args.molecules.split(","),
        methods=args.methods.split(","),
        seeds=args.seeds,
        n_steps=args.steps,
        sigma_deg=args.sigma,
        corr_lengths=tuple(float(x) for x in args.corr_lengths.split(",")),
        kT=args.kT,
        window=args.window,
        jobs=args.jobs,
        out_dir=args.out,
    )
    from perlin_mole.benchmark import format_markdown

    print(format_markdown(report))
    print(f"wrote {args.out}/report.json, summary.md, coverage.png")


def _report(args: argparse.Namespace) -> None:
    import json
    from pathlib import Path

    from perlin_mole.benchmark import format_markdown, pooled_tests, write_report

    out = Path(args.out)
    report = json.loads((out / "report.json").read_text())
    report["pooled"] = pooled_tests(report)
    write_report(report, out)
    print(format_markdown(report))


def _rescore(args: argparse.Namespace) -> None:
    from perlin_mole.benchmark import format_markdown, rescore

    print(format_markdown(rescore(args.out, window=args.window, jobs=args.jobs)))


def _landscape(args: argparse.Namespace) -> None:
    from perlin_mole.landscape import evaluation_cost, gradient_correlation

    print("Formulation A: N(q) over all n torsions")
    print(f"{'n':>3} {'us/point':>10} {'max |corr(dN/dq_i, dN/dq_j)|':>32}")
    for n in range(1, args.max_dim + 1):
        cost = evaluation_cost(n, n_points=max(50, 4000 >> max(0, n - 6)))
        c = gradient_correlation(n, n_points=4000) if n > 1 else np.ones((1, 1))
        off = np.abs(c[~np.eye(n, dtype=bool)]).max() if n > 1 else float("nan")
        print(f"{n:>3} {cost * 1e6:>10.1f} {off:>32.3f}")


def _molecules(_: argparse.Namespace) -> None:
    from perlin_mole.chem import TorsionSystem

    for name, smi in MOLECULES.items():
        s = TorsionSystem.from_smiles(smi)
        print(f"{name:12s} {s.n_dof:3d} torsions {s.n_atoms:4d} atoms  {smi}")


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(prog="perlin-mole", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("bench", help="equal-budget conformer search benchmark")
    b.add_argument("--molecules", default="decane,tetraglyme,ala3", help=f"names from {list(MOLECULES)} or SMILES")
    b.add_argument("--methods", default=",".join(METHODS))
    b.add_argument("--seeds", type=int, default=8)
    b.add_argument("--steps", type=int, default=400, help="MMFF minimisations per run")
    b.add_argument("--sigma", type=float, default=60.0, help="per-torsion kick std (deg)")
    b.add_argument("--corr-lengths", default="4", help="Perlin correlation lengths in bonds, comma-separated")
    b.add_argument("--kT", type=float, default=1.0, help="basin-hopping temperature (kcal/mol)")
    b.add_argument("--window", type=float, default=3.0, help="low-energy window (kcal/mol)")
    b.add_argument("--jobs", type=int, default=None)
    b.add_argument("--out", default="results/latest")
    b.set_defaults(func=_bench)

    r = sub.add_parser("rescore", help="re-score a saved benchmark (raw.pkl) without re-running minimisations")
    r.add_argument("--out", default="results/latest")
    r.add_argument("--window", type=float, default=3.0)
    r.add_argument("--jobs", type=int, default=None)
    r.set_defaults(func=_rescore)

    rp = sub.add_parser("report", help="regenerate summary.md / coverage.png / pooled tests from report.json")
    rp.add_argument("--out", default="results/latest")
    rp.set_defaults(func=_report)

    s = sub.add_parser("landscape", help="cost / gradient-correlation check of a full-dimensional Perlin landscape")
    s.add_argument("--max-dim", type=int, default=12)
    s.set_defaults(func=_landscape)

    m = sub.add_parser("molecules", help="list built-in benchmark molecules")
    m.set_defaults(func=_molecules)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
