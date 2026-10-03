"""Pre-registered confirmatory LJ38 test (see results/clusters_confirm/PREREGISTRATION.md)."""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # one BLAS thread per worker process; set before numpy is imported

import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact, wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_clusters import KT, STEP  # noqa: E402

from perlin_mole.clusters import basin_hopping  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results" / os.environ.get("LJ38_OUT", "clusters_confirm")
_lo, _hi = (int(x) for x in os.environ.get("LJ38_SEEDS", "30-130").split("-"))
SEEDS = range(_lo, _hi)
METHODS = os.environ.get("LJ38_METHODS", "perlin@2,gauss-field@2").split(",")
BUDGET = 3000


def run(args):
    m, sd = args
    return basin_hopping(38, m, sd, BUDGET, s=STEP, kT=KT)


if __name__ == "__main__":
    jobs = int(sys.argv[1]) if len(sys.argv) > 1 else None
    # resumable: every finished run is appended to runs.jsonl; finished (method, seed) pairs are skipped
    log = OUT / "runs.jsonl"
    done = [json.loads(line) for line in log.read_text().splitlines()] if log.exists() else []
    todo = [(m, s) for m in METHODS for s in SEEDS if (m, s) not in {(r["method"], r["seed"]) for r in done}]
    with ProcessPoolExecutor(max_workers=jobs) as pool, open(log, "a") as fh:
        futs = [pool.submit(run, t) for t in todo]
        for f in as_completed(futs):
            r = f.result()
            fh.write(json.dumps({k: r[k] for k in ("method", "seed", "t_hit", "best_e")}) + "\n")
            fh.flush()
            done.append(r)
    runs = done
    (OUT / "runs.json").write_text(json.dumps([{k: r[k] for k in ("method", "seed", "t_hit", "best_e")} for r in runs]))
    hit = {m: {r["seed"]: r["t_hit"] is not None for r in runs if r["method"] == m} for m in METHODS}
    k = {m: sum(hit[m].values()) for m in METHODS}
    n = len(SEEDS)
    p = fisher_exact([[k[METHODS[0]], n - k[METHODS[0]]], [k[METHODS[1]], n - k[METHODS[1]]]])[1]
    gap = {m: {r["seed"]: r["best_e"] - (-173.928427) for r in runs if r["method"] == m} for m in METHODS}
    d = np.array([gap[METHODS[0]][s] - gap[METHODS[1]][s] for s in SEEDS])
    res = {
        "successes": k,
        "n": n,
        "fisher_p_two_sided": float(p),
        "mean_gap_diff": float(d.mean()),
        "gap_wilcoxon_p": float(wilcoxon(d).pvalue) if np.any(d) else None,
    }
    (OUT / "result.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
