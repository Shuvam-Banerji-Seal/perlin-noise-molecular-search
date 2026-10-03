"""Exhaustive torsion-grid reference minima for the benchmark molecules (ground truth for recall).

Usage: uv run python scripts/systematic_reference.py [molecule ...]
Writes results/validation/<molecule>_systematic.npz and systematic.json.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

from perlin_mole.benchmark import KNOWN, MOLECULES
from perlin_mole.validate import systematic_search

STAGGERED = [60.0, 180.0, 300.0]
GRIDS = {
    "decane": STAGGERED,  # 3^7  = 2,187 starts
    "imatinib": [0.0, 90.0, 180.0, 270.0],  # 4^7  = 16,384 (aryl torsions are 2-fold)
    "ala3": [-180.0, -120.0, -60.0, 0.0, 60.0, 120.0],  # 6^6  = 46,656 (phi/psi need finer grid)
    "ala1": [float(a) for a in range(-180, 180, 15)],  # 24^2 = 576 (alanine dipeptide)
    "tetraglyme": STAGGERED,  # 3^12 = 531,441
}

out = Path("results/validation")
out.mkdir(parents=True, exist_ok=True)
summary_path = out / "systematic.json"
summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
for mol in sys.argv[1:] or list(GRIDS):
    t0 = time.perf_counter()
    e, heavy, full, n, n_unconv = systematic_search({**MOLECULES, **KNOWN}[mol], GRIDS[mol], return_unconverged=True)
    np.savez_compressed(out / f"{mol}_systematic.npz", energies=e, heavy=heavy, full=full)
    summary[mol] = {
        "grid_deg": GRIDS[mol],
        "n_starts": n,
        "n_unconverged_discarded": n_unconv,
        "n_unique": int(len(e)),
        "e_min": float(e[0]),
        "n_within_3": int((e <= e[0] + 3.0).sum()),
        "wall_s": round(time.perf_counter() - t0, 1),
    }
    print(mol, summary[mol], flush=True)
    summary_path.write_text(json.dumps(summary, indent=1))
