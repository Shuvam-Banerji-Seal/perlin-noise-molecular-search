"""Validation against independently known minima.

The benchmark's reference set is the union of what the compared methods
found, so it cannot reveal minima that *every* method missed. Here the
reference comes from an exhaustive systematic search instead: relax from every
point of a torsion grid. For small molecules this enumerates essentially all
minima, giving ground truth to measure recall against.
"""

from __future__ import annotations

import itertools
import os
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache

import numpy as np

from perlin_mole.chem import TorsionSystem
from perlin_mole.search import SymmetryRMSD, deduplicate

ALANINE_DIPEPTIDE = "CC(=O)N[C@@H](C)C(=O)NC"
# Backbone atoms in the SMILES above: phi = C(=O)-N-CA-C(=O), psi = N-CA-C(=O)-N
ALA_PHI, ALA_PSI = (1, 3, 4, 6), (3, 4, 6, 8)
# Canonical gas-phase basins of Ac-Ala-NHMe (approximate centres, degrees).
ALA_BASINS = {
    "C7eq": (-83.0, 73.0),
    "C5": (-155.0, 160.0),
    "C7ax": (74.0, -60.0),
    "beta2": (-140.0, 25.0),
    "alphaR": (-65.0, -40.0),
    "alphaL": (65.0, 35.0),
    "alpha'": (-165.0, -45.0),
}


# Features the torsion search does not change by design: stereocentres and secondary-amide cis/trans.
_AMIDE = "[OX1]=[CX3][NX3;H1][#6]"  # O=C-N-C: |dihedral| < 90 deg is the usual trans (Z) peptide-like amide


def _dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


@lru_cache(maxsize=64)
def _fixed_features(smiles: str) -> tuple[dict, tuple]:
    from rdkit import Chem

    from perlin_mole.isomers import cip_labels

    mol = Chem.MolFromSmiles(smiles)
    return cip_labels(smiles), tuple(mol.GetSubstructMatches(Chem.MolFromSmarts(_AMIDE)))


def keeps_fixed_features(smiles: str, heavy_xyz: np.ndarray) -> bool:
    """True if a minimum has the input's stereocentres (CIP labels) and only trans secondary amides.

    Torsion moves cannot invert a stereocentre or flip an amide, so minima that differ in these features are
    outside the searched space. Grid relaxations of clashed starts occasionally produce them; they are excluded
    from the ground truth (a stereo-inverted minimum belongs to a different stereoisomer).
    """
    from perlin_mole.isomers import cip_labels

    ref, amides = _fixed_features(smiles)
    if ref and cip_labels(smiles, heavy_xyz) != ref:
        return False
    return all(abs(_dihedral(*heavy_xyz[list(m)])) < 90.0 for m in amides)


def reference_targets(
    energies: np.ndarray, heavy: np.ndarray, smiles: str, window: float
) -> tuple[np.ndarray, float, int]:
    """Ground-truth targets of a grid: minima that keep the input's fixed features and lie within ``window`` of
    the lowest such minimum. Returns (indices, reference minimum energy, number of excluded minima in range)."""
    order = np.argsort(energies, kind="stable")
    keep, e0, excluded = [], None, 0
    for i in order:
        if e0 is not None and energies[i] > e0 + window:
            break
        if keeps_fixed_features(smiles, heavy[i]):
            e0 = energies[i] if e0 is None else e0
            keep.append(int(i))
        else:
            excluded += 1
    return np.asarray(keep, dtype=int), float(e0), excluded


@lru_cache(maxsize=4)
def _sys(smiles: str) -> TorsionSystem:
    return TorsionSystem.from_smiles(smiles, seed=0)


def _relax_from(args: tuple[str, np.ndarray]) -> tuple[float, np.ndarray, np.ndarray, bool]:
    smiles, angles = args
    s = _sys(smiles)
    start = s.initial_xyz  # fixed reference geometry: every grid point starts from the same embedding
    xyz = s.apply_delta(start, angles - s.torsion_values(start))
    xyz, e, conv = s.minimize(xyz)
    return e, xyz[s.heavy_idx], xyz, conv


def systematic_search(smiles: str, values_deg: list[float], jobs: int | None = None, return_unconverged: bool = False):
    """Relax from every grid point; return unique (energies, heavy_xyz, full_xyz, n_starts) sorted by energy.

    Relaxations that hit the iteration limit are discarded: an unconverged geometry is
    not a minimum and would otherwise show up as a spurious "unique" entry.
    """
    s = _sys(smiles)
    grid = [np.deg2rad(np.array(p)) for p in itertools.product(values_deg, repeat=s.n_dof)]
    jobs = jobs or os.cpu_count() or 1
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        out = list(pool.map(_relax_from, [(smiles, g) for g in grid], chunksize=max(1, len(grid) // (jobs * 8))))
    n_unconverged = sum(not o[3] for o in out)
    out = [o for o in out if o[3]]
    e = np.array([o[0] for o in out])
    order = np.argsort(e, kind="stable")
    e = e[order]
    heavy = np.stack([out[i][1] for i in order])
    full = np.stack([out[i][2] for i in order])
    rep = deduplicate(e, heavy, SymmetryRMSD(s.heavy_mol()), jobs=jobs)
    keep = rep == np.arange(len(rep))
    if return_unconverged:
        return e[keep], heavy[keep], full[keep], int(len(grid)), int(n_unconverged)
    return e[keep], heavy[keep], full[keep], int(len(grid))


def recall(
    ref_e: np.ndarray, ref_xyz: np.ndarray, found: list[tuple[float, np.ndarray]], rmsd: SymmetryRMSD, window: float
) -> tuple[float, list[int]]:
    """Fraction of reference minima within ``window`` of the reference best that appear in ``found``."""
    low = np.flatnonzero(ref_e <= ref_e.min() + window)
    fe = np.array([f[0] for f in found])
    fx = np.stack([f[1] for f in found])
    hit = []
    for i in low:
        cand = np.flatnonzero(np.abs(fe - ref_e[i]) < 0.1)
        if cand.size and (rmsd.many(ref_xyz[i], fx[cand]) < 0.25).any():
            hit.append(int(i))
    return len(hit) / len(low), hit


def classify_ala(phi_deg: float, psi_deg: float) -> str:
    def d(a: float, b: float) -> float:
        return abs((a - b + 180.0) % 360.0 - 180.0)

    name, (p0, s0) = min(ALA_BASINS.items(), key=lambda kv: np.hypot(d(phi_deg, kv[1][0]), d(psi_deg, kv[1][1])))
    return name if np.hypot(d(phi_deg, p0), d(psi_deg, s0)) < 50.0 else "other"
