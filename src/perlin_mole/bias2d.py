"""Formulation A1 tested where it is affordable: a Perlin *bias* on a 2-D torsional free-energy surface.

Alanine dipeptide's (phi, psi) surface is built by a relaxed MMFF94 scan (all
other coordinates relaxed under torsion constraints). Metropolis Monte Carlo
on the interpolated surface then compares, at 300 K:

* unbiased MC;
* static bias  ``lambda * N(phi, psi)``, with N either Perlin or a Gaussian
  random field whose power spectrum equals Perlin's (second-order matched);
* drifting bias, where N moves slowly in time (Perlin: a third, time-like
  noise coordinate; GRF: AR(1) Fourier coefficients with matched lag
  correlation);
* well-tempered metadynamics as the standard enhanced-sampling baseline.

Metrics: committed transitions between the phi<0 and phi>0 basins (sampling
speed), and the error of the reweighted free-energy difference against the
exact value from the surface itself (sampling correctness).
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache

import numpy as np
from rdkit.Chem import rdForceFieldHelpers as ffh
from scipy.interpolate import RegularGridInterpolator

from perlin_mole.chem import TorsionSystem
from perlin_mole.noise import FractalNoise
from perlin_mole.validate import ALA_PHI, ALA_PSI, ALANINE_DIPEPTIDE

KT_300 = 0.0019872 * 300.0  # kcal/mol


@lru_cache(maxsize=2)
def _sys(seed: int) -> TorsionSystem:
    return TorsionSystem.from_smiles(ALANINE_DIPEPTIDE, seed=seed)


def _constrained(args: tuple[float, float, int]) -> float:
    phi, psi, seed = args
    s = _sys(seed)
    x = s.initial_xyz.copy()
    s._load(x)
    from rdkit.Chem import rdMolTransforms

    conf = s.mol.GetConformer()
    rdMolTransforms.SetDihedralDeg(conf, *ALA_PHI, phi)
    rdMolTransforms.SetDihedralDeg(conf, *ALA_PSI, psi)
    ff = ffh.MMFFGetMoleculeForceField(s.mol, s.props)
    for t, v in ((ALA_PHI, phi), (ALA_PSI, psi)):
        ff.MMFFAddTorsionConstraint(*t, False, v, v, 1.0e4)
    ff.Minimize(maxIts=4000)
    return ffh.MMFFGetMoleculeForceField(s.mol, s.props).CalcEnergy()


def relaxed_scan(
    step_deg: float = 10.0, seeds: tuple[int, ...] = tuple(range(16)), jobs: int | None = None
) -> np.ndarray:
    """E(phi, psi) on a periodic grid; minimum over several start embeddings.

    Constrained relaxations depend strongly on the start geometry. With 3 start embeddings, 174 of the
    1296 grid points were left >1 kcal/mol above their converged value (and the exact dF was 1.06 instead
    of 2.40 kcal/mol). 8 and 16 embeddings agree on dF to 1e-3, so 16 is the default.
    """
    ang = np.arange(-180.0, 180.0, step_deg)
    tasks = [(p, q, s) for s in seeds for p in ang for q in ang]
    with ProcessPoolExecutor(max_workers=jobs or os.cpu_count()) as pool:
        e = np.array(list(pool.map(_constrained, tasks, chunksize=16)))
    e = e.reshape(len(seeds), len(ang), len(ang)).min(axis=0)
    return e - e.min()


class Surface:
    """Periodic cubic interpolation of a grid E[i, j] over (phi, psi) in degrees."""

    def __init__(self, grid: np.ndarray, step_deg: float) -> None:
        n = grid.shape[0]
        pad = 3
        g = np.pad(grid, pad, mode="wrap")
        ax = -180.0 + step_deg * (np.arange(n + 2 * pad) - pad)
        self._f = RegularGridInterpolator((ax, ax), g, method="cubic")

    def __call__(self, phi: np.ndarray, psi: np.ndarray) -> np.ndarray:
        w = lambda a: (a + 180.0) % 360.0 - 180.0  # noqa: E731
        return self._f(np.stack([w(phi), w(psi)], axis=-1))


def exact_dF(surface: Surface, kT: float = KT_300, n: int = 720) -> float:
    """F(phi>0) - F(phi<0) by direct Boltzmann integration of the surface."""
    a = np.linspace(-180, 180, n, endpoint=False) + 180.0 / n  # cell midpoints: symmetric, none on phi = 0
    P, Q = np.meshgrid(a, a, indexing="ij")
    w = np.exp(-surface(P, Q) / kT)
    return float(-kT * np.log(w[P > 0].sum() / w[P < 0].sum()))


# ------------------------------------------------------------------ bias fields on the torus
class PerlinBias:
    """lambda * N(phi, psi[, t]) with N fBm Perlin, periodic on the torus; optionally drifting in time."""

    def __init__(self, strength: float, cells: int = 4, drift: float = 0.0, seed: int = 0) -> None:
        self.strength, self.cells, self.drift = strength, cells, drift
        self.field = FractalNoise(3, octaves=2, period=cells, seed=seed)
        rng = np.random.default_rng(seed + 7)
        self.z0 = rng.uniform(0, cells)
        a = rng.uniform(0, 360, size=(40000, 2))
        self.scale = float(self.field(self._pts(a[:, 0], a[:, 1], 0)).std())

    def _pts(self, phi, psi, t):
        u = (np.asarray(phi) + 180.0) / 360.0 * self.cells
        v = (np.asarray(psi) + 180.0) / 360.0 * self.cells
        z = np.full_like(u, self.z0 + self.drift * t)
        return np.stack([u, v, z], axis=-1)

    def __call__(self, phi, psi, t: float = 0.0):
        return self.strength * self.field(self._pts(phi, psi, t)) / self.scale


class MatchedGRFBias:
    """Gaussian random field on the torus with Perlin's empirical power spectrum (unit variance), times lambda.

    ``rho`` > 0 makes the Fourier coefficients an AR(1) process in time (drifting bias).
    """

    def __init__(self, perlin: PerlinBias, n: int = 128, rho: float = 0.0, seed: int = 0, n_spec: int = 64) -> None:
        self.n, self.rho, self.strength = n, rho, perlin.strength
        a = np.linspace(-180, 180, n, endpoint=False)
        P, Q = np.meshgrid(a, a, indexing="ij")
        spec = np.zeros((n, n))
        rng = np.random.default_rng(seed + 11)
        for _ in range(n_spec):  # average |FFT|^2 over random time slices of the Perlin field
            z = rng.uniform(0, 1e3)
            f = perlin.field(
                np.stack([(P + 180) / 360 * perlin.cells, (Q + 180) / 360 * perlin.cells, np.full_like(P, z)], -1)
            )
            spec += np.abs(np.fft.fft2(f - f.mean())) ** 2
        self.amp = np.sqrt(spec / n_spec)
        self.rng = np.random.default_rng(seed)
        self._white = np.fft.fft2(self.rng.standard_normal((n, n)))
        self._make()

    def _make(self) -> None:
        g = np.real(np.fft.ifft2(self._white * self.amp))
        g /= g.std()
        self.grid = Surface(g, 360.0 / self.n)

    def step(self) -> None:
        if self.rho > 0:
            new = np.fft.fft2(self.rng.standard_normal((self.n, self.n)))
            self._white = self.rho * self._white + np.sqrt(1 - self.rho**2) * new
            self._make()

    def __call__(self, phi, psi, t: float = 0.0):
        return self.strength * self.grid(phi, psi)


# ------------------------------------------------------------------ Monte Carlo
def run_mc(
    surface: Surface,
    bias=None,
    n_steps: int = 200_000,
    n_walkers: int = 32,
    step_deg: float = 15.0,
    kT: float = KT_300,
    update_every: int = 200,
    metad: dict | None = None,
    seed: int = 0,
) -> dict:
    """Vectorised Metropolis walkers on E + bias. Returns transitions and reweighted dF."""
    rng = np.random.default_rng(seed)
    phi = np.full(n_walkers, -80.0) + rng.normal(0, 5, n_walkers)  # start in C7eq
    psi = np.full(n_walkers, 75.0) + rng.normal(0, 5, n_walkers)
    # metadynamics: hills are deposited onto a periodic per-walker grid (2-degree cells) and read back by
    # bilinear interpolation, so the cost per step is O(1) instead of O(number of hills)
    ng, cell = 180, 2.0
    mid = -180.0 + cell * (np.arange(ng) + 0.5)
    grid = np.zeros((n_walkers, ng, ng))
    walkers = np.arange(n_walkers)

    def grid_value(p, q):
        x = (np.asarray(p) + 180.0) / cell - 0.5
        y = (np.asarray(q) + 180.0) / cell - 0.5
        i0, j0 = np.floor(x).astype(int), np.floor(y).astype(int)
        fx, fy = x - i0, y - j0
        i0, j0 = i0 % ng, j0 % ng
        i1, j1 = (i0 + 1) % ng, (j0 + 1) % ng
        return (
            (1 - fx) * (1 - fy) * grid[walkers, i0, j0]
            + fx * (1 - fy) * grid[walkers, i1, j0]
            + (1 - fx) * fy * grid[walkers, i0, j1]
            + fx * fy * grid[walkers, i1, j1]
        )

    def V(p, q, t):
        v = np.zeros_like(p)
        if bias is not None:
            v = v + bias(p, q, t)
        if metad is not None:
            v = v + grid_value(p, q)
        return v

    state = np.where(phi < 0, 0, 1)
    core_a = lambda p: (p > -170) & (p < -40)  # noqa: E731
    core_b = lambda p: (p > 40) & (p < 110)  # noqa: E731
    transitions = np.zeros(n_walkers, int)
    w_pos = np.zeros(n_walkers)
    w_neg = np.zeros(n_walkers)
    e_cur = surface(phi, psi) + V(phi, psi, 0)
    for t in range(n_steps):
        tb = (t // update_every) * update_every  # bias is piecewise constant: consistent Metropolis ratios
        if t % update_every == 0 and t:
            if hasattr(bias, "step"):
                bias.step()
            if metad is not None:
                h = metad["height"] * np.exp(-grid_value(phi, psi) / (kT * (metad["gamma"] - 1)))
                dp = (mid[None, :] - phi[:, None] + 180) % 360 - 180
                dq = (mid[None, :] - psi[:, None] + 180) % 360 - 180
                w2 = 2 * metad["width"] ** 2
                grid += h[:, None, None] * np.exp(-(dp[:, :, None] ** 2) / w2) * np.exp(-(dq[:, None, :] ** 2) / w2)
            e_cur = surface(phi, psi) + V(phi, psi, tb)
        p2 = (phi + rng.normal(0, step_deg, n_walkers) + 180) % 360 - 180
        q2 = (psi + rng.normal(0, step_deg, n_walkers) + 180) % 360 - 180
        e2 = surface(p2, q2) + V(p2, q2, tb)
        acc = rng.random(n_walkers) < np.exp(-(e2 - e_cur) / kT)
        phi, psi, e_cur = np.where(acc, p2, phi), np.where(acc, q2, psi), np.where(acc, e2, e_cur)
        new_a, new_b = core_a(phi), core_b(phi)
        transitions += (state == 0) & new_b
        transitions += (state == 1) & new_a
        state = np.where(new_a, 0, np.where(new_b, 1, state))
        # reweighting with the instantaneous bias (exact for static bias; approximate for slow drift)
        w = np.exp(V(phi, psi, tb) / kT) if (bias is not None and metad is None) else np.ones(n_walkers)
        if metad is None:
            w_pos += w * (phi > 0)
            w_neg += w * (phi < 0)
    if metad is not None:
        # well-tempered estimator: F = -gamma/(gamma-1) V on the grid at the end of the run
        P = np.broadcast_to(mid[:, None], (ng, ng))
        fes = []
        for k in range(n_walkers):
            f = -metad["gamma"] / (metad["gamma"] - 1) * grid[k]
            wb = np.exp(-(f - f.min()) / kT)
            fes.append(-kT * np.log(wb[P > 0].sum() / wb[P < 0].sum()))
        df = np.array(fes)
    else:
        with np.errstate(divide="ignore"):
            df = -kT * np.log(w_pos / w_neg)
    return {"transitions": transitions, "dF": df}
