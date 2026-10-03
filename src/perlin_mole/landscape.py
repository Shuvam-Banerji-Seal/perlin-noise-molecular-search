"""Formulation A: a Perlin landscape N(q) over the full torsion space.

This is the ``E_total(q) = E_chem(q) + lambda * N(q)`` / ``dq = -eta grad N(q)``
version from the original proposal. Two properties decide whether it is
practical, and both are measured here:

1. **Cost.** Lattice gradient noise visits 2^n hypercube corners per point,
   so a landscape over n torsions costs O(2^n n). See ``evaluation_cost``.
2. **What "correlated" means.** For a stationary, statistically isotropic
   field, the gradient components at a single point are uncorrelated:
   corr(dN/dq_i, dN/dq_j) ~ 0 for i != j. So ``-grad N`` does *not* make
   different torsions move together. It is correlated *in q* (nearby states
   get similar drifts), which is a smooth random bias, not a collective
   move. See ``gradient_correlation``. The collective-move behaviour the
   proposal wants comes from formulation B
   (``proposals.PerlinFieldProposal``): a low-dimensional field read at each
   DOF's position in the molecule.
"""

from __future__ import annotations

import time

import numpy as np

from perlin_mole.noise import FractalNoise, PerlinNoise


class TorsionPerlinBias:
    """lambda * N(phi * P / 2pi): a smooth random bias on the k-torus of torsions."""

    def __init__(self, n_dof: int, strength: float, cells_per_turn: int = 4, octaves: int = 1, seed: int = 0) -> None:
        self.strength = strength
        self.cells = cells_per_turn
        self.field = FractalNoise(n_dof, octaves=octaves, period=cells_per_turn, seed=seed)

    def __call__(self, phi: np.ndarray, grad: bool = False):
        x = np.mod(phi, 2 * np.pi) * (self.cells / (2 * np.pi))
        if not grad:
            return self.strength * self.field(x)
        v, g = self.field(x, grad=True)
        return self.strength * v, self.strength * g * (self.cells / (2 * np.pi))


def gradient_correlation(ndim: int, n_points: int = 20000, seed: int = 0) -> np.ndarray:
    """Correlation matrix of grad N components at random points."""
    noise = PerlinNoise(ndim, period=64, seed=seed)
    x = np.random.default_rng(seed).uniform(0, 64, size=(n_points, ndim))
    _, g = noise(x, grad=True)
    return np.corrcoef(g, rowvar=False)


def evaluation_cost(ndim: int, n_points: int = 2000, seed: int = 0) -> float:
    """Seconds per point for one value+gradient evaluation of n-D Perlin noise."""
    noise = PerlinNoise(ndim, period=64, seed=seed)
    x = np.random.default_rng(seed).uniform(0, 64, size=(n_points, ndim))
    t0 = time.perf_counter()
    noise(x, grad=True)
    return (time.perf_counter() - t0) / n_points
