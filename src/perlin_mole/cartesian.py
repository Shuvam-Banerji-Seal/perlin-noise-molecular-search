"""Perlin noise applied literally to Cartesian coordinates of a molecule.

The original idea, ``R -> R + eta * N(R)``, applied to atoms and then relaxed with
MMFF94. Reuses the 3-D displacement fields from :mod:`perlin_mole.clusters`
(Perlin; a Gaussian random field with Perlin's isotropic correlation k(r/L);
i.i.d.). Unlike torsion moves, Cartesian kicks can change chemistry after
relaxation: bonds can break or form, and stereocentres can invert. Every relaxed
structure is therefore checked for **chemical integrity** (same heavy-atom
connectivity and same CIP labels as the input). Only intact minima enter the
archive.
"""

from __future__ import annotations

import time

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

from perlin_mole.benchmark import Task, _start, _system
from perlin_mole.clusters import make_proposal
from perlin_mole.isomers import cip_labels
from perlin_mole.search import ConformerArchive, SymmetryRMSD


def _heavy_graph(mol: Chem.Mol) -> frozenset:
    return frozenset(
        tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
        for b in mol.GetBonds()
        if b.GetBeginAtom().GetAtomicNum() > 1 and b.GetEndAtom().GetAtomicNum() > 1
    )


def intact(system, xyz: np.ndarray, ref_graph: frozenset, ref_cip: dict) -> bool:
    """Same heavy-atom connectivity (from distances) and same stereo labels as the input molecule."""
    m = Chem.RWMol(system.mol)
    m.RemoveAllConformers()
    for b in list(m.GetBonds()):
        m.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
    conf = Chem.Conformer(m.GetNumAtoms())
    for i, p in enumerate(xyz):
        conf.SetAtomPosition(i, p.tolist())
    m.AddConformer(conf, assignId=True)
    try:
        rdDetermineBonds.DetermineConnectivity(m, useHueckel=False)
    except Exception:  # noqa: BLE001
        return False
    if _heavy_graph(m) != ref_graph:
        return False
    return not ref_cip or cip_labels(system.smiles, xyz[system.heavy_idx]) == ref_cip


def run_cartesian_task(task: Task) -> dict:
    """Same output format as benchmark.run_task. ``task.method`` is e.g. 'cart-perlin@2', 'cart-iid'.

    ``task.sigma_deg`` is reused as the per-coordinate displacement std in units of 0.01 Angstrom.
    """
    system = _system(task.smiles)
    ref_graph = _heavy_graph(system.mol)
    ref_cip = cip_labels(system.smiles)
    s = task.sigma_deg / 100.0
    method = task.method.removeprefix("cart-").replace("iid", "gauss-iid")
    prop = make_proposal(method, s, task.seed)
    archive = ConformerArchive(SymmetryRMSD(system.heavy_mol()))
    rng = np.random.default_rng([task.seed, 977])
    t0 = time.perf_counter()
    xyz, e, _ = system.minimize(_start(system, task.seed))
    archive.add(xyz[system.heavy_idx], e, 0)
    best, accepted, broken, unconverged = [e], 0, 0, 0
    for step in range(1, task.n_steps + 1):
        trial, e_t, conv = system.minimize(xyz + prop.sample(xyz, rng))
        unconverged += not conv
        if not intact(system, trial, ref_graph, ref_cip):
            broken += 1
            best.append(best[-1])
            continue
        archive.add(trial[system.heavy_idx], e_t, step)
        if e_t <= e or rng.random() < np.exp(-(e_t - e) / task.kT):
            xyz, e = trial, e_t
            accepted += 1
        best.append(min(best[-1], e_t))
    from dataclasses import asdict

    return {
        "task": asdict(task),
        "label": task.method,
        "best_trace": best,
        "minima": [(m.energy, m.heavy_xyz, m.step) for m in archive.minima],
        "acceptance": accepted / task.n_steps,
        "unconverged": unconverged,
        "broken": broken,
        "wall_time": time.perf_counter() - t0,
    }
