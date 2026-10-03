"""Unrelaxed (rigid-rotor) Metropolis MC in torsion space: where higher-order kick statistics could matter.

Basin hopping relaxes every proposal, which discards most of the kick's fine
structure. Here there is no relaxation. Each proposal sets the torsions and is
accepted on the *instantaneous* MMFF94 energy, so the chain samples
exp(-E/kT) restricted to the torsion manifold (bond lengths and angles fixed at
the start geometry). All proposals are symmetric and state-independent, so all
chains target the same distribution, which we check via the mean energy.
Efficiency metrics:

* acceptance rate;
* ESJD: expected squared (circular) jump distance in torsion space per step;
* rotamer transitions: changes of any torsion's {g+, t, g-} bin per 1000 steps.
"""

from __future__ import annotations

import time

import numpy as np
from rdkit.Chem import rdForceFieldHelpers as ffh
from rdkit.Chem import rdMolTransforms

from perlin_mole.benchmark import Task, _start, _system
from perlin_mole.proposals import build_proposal


def _bins(phi: np.ndarray) -> np.ndarray:
    d = np.degrees(phi) % 360.0
    return np.where(d < 120, 0, np.where(d < 240, 1, 2))


def run_rigid_task(task: Task) -> dict:
    s = _system(task.smiles)
    xyz, _, _ = s.minimize(_start(s, task.seed))  # relax once: the rigid bond/angle geometry
    s._load(xyz)
    conf = s.mol.GetConformer()
    ff = ffh.MMFFGetMoleculeForceField(s.mol, s.props)
    prop = build_proposal(
        task.method,
        s.dof_coordinates(3),
        np.deg2rad(task.sigma_deg),
        task.corr_length,
        seed=task.seed,
        n_steps=task.n_steps,
    )
    rng = np.random.default_rng([task.seed, 4242])
    kT = task.kT
    phi = np.array([rdMolTransforms.GetDihedralRad(conf, *t) for t in s.torsions])
    pos = conf.GetPositions()
    e = ff.CalcEnergy(pos.ravel().tolist())
    b = _bins(phi)
    acc = 0
    esjd = 0.0
    trans = 0
    e_sum = 0.0
    t0 = time.perf_counter()
    for _ in range(task.n_steps):
        d = prop.sample(rng)
        for a, p in enumerate(pos):
            conf.SetAtomPosition(a, p.tolist())
        for t, phi_new in zip(s.torsions, phi + d, strict=True):
            rdMolTransforms.SetDihedralRad(conf, *t, float(phi_new))
        new_pos = conf.GetPositions()
        e_new = ff.CalcEnergy(new_pos.ravel().tolist())
        if e_new <= e or rng.random() < np.exp(-(e_new - e) / kT):
            acc += 1
            esjd += float(np.sum(np.angle(np.exp(1j * d)) ** 2))
            phi = np.angle(np.exp(1j * (phi + d)))
            pos, e = new_pos, e_new
            nb = _bins(phi)
            trans += int(np.any(nb != b))
            b = nb
        e_sum += e
    n = task.n_steps
    return {
        "molecule": task.molecule,
        "method": task.method,
        "seed": task.seed,
        "sigma_deg": task.sigma_deg,
        "acceptance": acc / n,
        "esjd": esjd / n,
        "transitions_per_1000": 1000.0 * trans / n,
        "mean_energy": e_sum / n,
        "wall_time": time.perf_counter() - t0,
    }
