"""Torsion-space move generators for basin hopping.

The default proposals here are *state independent* and *symmetric*
(P(delta) == P(-delta)), so plain Metropolis acceptance is valid for them and the
comparison is apples to apples. The time-dependent ablation variants (drifting field,
its AR(1) Gaussian twin, coarse-to-fine schedules) are symmetric only marginally:
they are conditioned on the previous kick or on the step number.

The question the benchmark asks is whether a Perlin field helps, and there
are two separate effects to tell apart:

* correlation between DOFs (neighbouring torsions rotate together), and
* everything else about Perlin noise (bounded, non-Gaussian marginals,
  lattice structure).

``CorrelatedGaussianProposal.matched_to(perlin)`` has exactly the Perlin
field's correlation matrix and marginal scale. If Perlin cannot beat it,
whatever helps is the correlation, and a Gaussian with any reasonable kernel
would do the same job.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.special import ndtr, ndtri
from scipy.stats import special_ortho_group

from perlin_mole.noise import FractalNoise


@dataclass
class RandomResetProposal:
    """MCMM-style move: re-draw 1..max_changed random torsions uniformly on the circle."""

    n_dof: int
    max_changed: int = 3
    name: str = "reset"

    def sample(self, rng: np.random.Generator) -> np.ndarray:
        delta = np.zeros(self.n_dof)
        m = rng.integers(1, min(self.max_changed, self.n_dof) + 1)
        idx = rng.choice(self.n_dof, size=m, replace=False)
        delta[idx] = rng.uniform(-np.pi, np.pi, size=m)
        return delta


@dataclass
class GaussianProposal:
    """Independent N(0, sigma^2) kick on every torsion."""

    n_dof: int
    sigma: float
    name: str = "gauss-iid"

    def sample(self, rng: np.random.Generator) -> np.ndarray:
        return rng.normal(0.0, self.sigma, self.n_dof)


@dataclass
class CorrelatedGaussianProposal:
    """N(0, sigma^2 * C) with a given correlation matrix C.

    Ablation options (defaults reproduce the plain matched Gaussian exactly):

    * ``marginal_quantiles``: push each coordinate through Phi and then the
      given empirical quantile function (a Gaussian copula with Perlin
      marginals).
    * ``ar``: AR(1) coefficient across successive calls (a time-coherent
      Gaussian, the matched control for ``PerlinFieldProposal(drift>0)``).
    """

    corr: np.ndarray
    sigma: float
    name: str = "gauss-corr"
    marginal_quantiles: np.ndarray | None = None
    ar: float = 0.0
    schedule: list[np.ndarray] | None = None  # correlation matrices for successive stages of the run
    n_steps: int = 0
    _root: np.ndarray = field(init=False, repr=False)
    _prev: np.ndarray | None = field(default=None, init=False, repr=False)
    _t: int = field(default=0, init=False, repr=False)

    @staticmethod
    def _sqrt(c: np.ndarray) -> np.ndarray:
        w, v = np.linalg.eigh(c)
        return v * np.sqrt(np.clip(w, 0.0, None))

    def __post_init__(self) -> None:
        self._root = self._sqrt(self.corr)
        self._roots = [self._sqrt(c) for c in self.schedule] if self.schedule else None

    def sample(self, rng: np.random.Generator) -> np.ndarray:
        root = self._root
        if self._roots:
            root = self._roots[_stage(self._t, self.n_steps, len(self._roots))]
            self._t += 1
        z = root @ rng.standard_normal(root.shape[1])
        if self.ar > 0.0:
            z = z if self._prev is None else self.ar * self._prev + np.sqrt(1.0 - self.ar**2) * z
            self._prev = z
        if self.marginal_quantiles is not None:
            q = self.marginal_quantiles
            z = np.interp(ndtr(z) * (len(q) - 1), np.arange(len(q)), q)
        return self.sigma * z

    @classmethod
    def matched_to(cls, perlin: PerlinFieldProposal, **kw) -> CorrelatedGaussianProposal:
        return cls(perlin.correlation(), perlin.sigma, **kw)


@dataclass
class PerlinFieldProposal:
    """Spatially coherent kicks from a fractal Perlin field over DOF positions.

    Each torsion j has a fixed position ``coords[j]`` (topological embedding,
    in bonds). A proposal samples a random translate and rotation of the field
    plus a random sign, then reads it at every torsion:

        delta_j = sigma * s * N(R coords[j] / L + o) / std(N)

    Torsions closer than ~L bonds move together. Random (R, o, s) makes the
    draw independent of the current state and symmetric.

    Ablation options (defaults reproduce the main-benchmark proposal exactly):

    * ``rotate=False``: no random rotation (lattice-aligned field).
    * ``shuffle_sites=True``: randomly permute which torsion sits where, so
      correlation no longer follows the molecular topology.
    * ``marginal="gaussian"``: map each value through the field's empirical
      CDF and then Phi^-1. This keeps Perlin's dependence structure (copula)
      but gives Gaussian marginals.
    * ``drift>0``: keep (R, s) fixed and random-walk the offset by ``drift``
      field units per call, so successive kicks are correlated in time (the
      animation use of Perlin noise). ``lag1`` gives the resulting lag-1
      autocorrelation.
    """

    coords: np.ndarray
    sigma: float
    corr_length: float = 2.0
    octaves: int = 2
    persistence: float = 0.5
    period: int = 64
    seed: int = 0
    n_calibration: int = 8192
    name: str = "perlin"
    rotate: bool = True
    shuffle_sites: bool = False
    marginal: str = "perlin"
    drift: float = 0.0
    anneal_to: float = 0.0  # >0: shrink L geometrically from corr_length to anneal_to over n_steps
    n_steps: int = 0
    n_stages: int = 8
    _field: FractalNoise = field(init=False, repr=False)
    _scale: float = field(init=False, repr=False)
    _state: list | None = field(default=None, init=False, repr=False)
    _t: int = field(default=0, init=False, repr=False)

    def __post_init__(self) -> None:
        self.coords = np.atleast_2d(np.asarray(self.coords, dtype=np.float64))
        if self.shuffle_sites:
            self.coords = self.coords[np.random.default_rng(self.seed + 2).permutation(self.coords.shape[0])]
        self._field = FractalNoise(
            self.coords.shape[1],
            octaves=self.octaves,
            persistence=self.persistence,
            period=self.period,
            seed=self.seed,
        )
        cal = self._raw(np.random.default_rng(self.seed + 1), self.n_calibration)
        self._scale = float(cal.std())
        self._calibration = cal / self._scale
        self._sorted = np.sort(self._calibration.ravel())
        self.lag1 = self._lag1() if self.drift > 0 else 0.0

    @property
    def n_dof(self) -> int:
        return self.coords.shape[0]

    def _rot(self, rng: np.random.Generator, n: int) -> np.ndarray:
        d = self.coords.shape[1]
        if d > 1 and self.rotate:
            return special_ortho_group.rvs(d, size=n, random_state=rng).reshape(n, d, d)
        return np.broadcast_to(np.eye(d), (n, d, d))

    def stage_lengths(self) -> list[float]:
        """Correlation length used in each stage of an annealed run (geometric from corr_length to anneal_to)."""
        return list(np.geomspace(self.corr_length, self.anneal_to, self.n_stages))

    def _raw(self, rng: np.random.Generator, n: int, length: float | None = None) -> np.ndarray:
        d = self.coords.shape[1]
        rot = self._rot(rng, n)
        offset = rng.uniform(0.0, self.period, size=(n, 1, d))
        pts = np.einsum("nij,kj->nki", rot, self.coords / (length or self.corr_length)) + offset
        sign = rng.choice([-1.0, 1.0], size=(n, 1))
        return sign * self._field(pts)

    def _drifting(self, rng: np.random.Generator) -> np.ndarray:
        d = self.coords.shape[1]
        if self._state is None:
            self._state = [self._rot(rng, 1)[0], rng.uniform(0.0, self.period, size=d), rng.choice([-1.0, 1.0])]
        else:
            step = rng.standard_normal(d)
            self._state[1] = self._state[1] + self.drift * step / np.linalg.norm(step)
        rot, off, sign = self._state
        return sign * self._field((self.coords / self.corr_length) @ rot.T + off)

    def _lag1(self, n: int = 3000) -> float:
        saved, self._state = self._state, None
        rng = np.random.default_rng(self.seed + 3)
        chain = np.array([self._drifting(rng) for _ in range(n)])
        self._state = saved
        return float(np.mean([np.corrcoef(chain[:-1, j], chain[1:, j])[0, 1] for j in range(self.n_dof)]))

    def sample(self, rng: np.random.Generator) -> np.ndarray:
        if self.anneal_to > 0:
            length = self.stage_lengths()[_stage(self._t, self.n_steps, self.n_stages)]
            self._t += 1
            return self.sigma * self._raw(rng, 1, length)[0] / self._scale
        raw = self._drifting(rng) if self.drift > 0 else self._raw(rng, 1)[0]
        if self.marginal == "gaussian":
            u = (np.searchsorted(self._sorted, raw / self._scale) + 0.5) / (len(self._sorted) + 1)
            return self.sigma * ndtri(u)
        # keep this exact operation order: published runs are bit-reproducible with it
        return self.sigma * raw / self._scale

    def correlation(self) -> np.ndarray:
        """Empirical correlation matrix of the (unit-scaled) field across DOFs."""
        if self.n_dof == 1:
            return np.ones((1, 1))
        c = np.corrcoef(self._calibration, rowvar=False)
        return np.nan_to_num(c, nan=0.0) + 1e-9 * np.eye(self.n_dof)


def _stage(t: int, n_steps: int, n_stages: int) -> int:
    return min(n_stages - 1, t * n_stages // max(n_steps, 1))


# Ablation registry: name -> (Perlin field options, which proposal to return).
# "gauss-corr*" variants are Gaussians matched to the correlation of the Perlin field with the same options.
ABLATIONS: dict[str, tuple[dict, str]] = {
    "perlin": ({}, "perlin"),
    "gauss-corr": ({}, "gauss"),
    "perlin-gaussmarg": ({"marginal": "gaussian"}, "perlin"),
    "gauss-perlinmarg": ({}, "gauss-perlin-marginals"),
    "perlin-oct1": ({"octaves": 1}, "perlin"),
    "gauss-corr-oct1": ({"octaves": 1}, "gauss"),
    "perlin-oct4": ({"octaves": 4}, "perlin"),
    "gauss-corr-oct4": ({"octaves": 4}, "gauss"),
    "perlin-norot": ({"rotate": False}, "perlin"),
    "gauss-corr-norot": ({"rotate": False}, "gauss"),
    "perlin-shuffled": ({"shuffle_sites": True}, "perlin"),
    "gauss-corr-shuffled": ({"shuffle_sites": True}, "gauss"),
    "perlin-drift": ({"drift": 0.25}, "perlin"),
    "gauss-ar1": ({"drift": 0.25}, "gauss-ar1"),
    # coarse-to-fine: L shrinks from the task's corr_length to 1 bond over the run (8 stages)
    "perlin-anneal": ({"anneal_to": 1.0}, "perlin"),
    "gauss-anneal": ({"anneal_to": 1.0}, "gauss-anneal"),
}


def build_proposal(
    method: str,
    dof_coords: np.ndarray,
    sigma: float,
    corr_length: float,
    seed: int = 0,
    n_steps: int = 0,
):
    """Factory used by the benchmark. ``dof_coords`` is (n_dof, d). ``n_steps`` is needed by annealed variants."""
    n_dof = dof_coords.shape[0]
    if method == "reset":
        return RandomResetProposal(n_dof)
    if method == "gauss-iid":
        return GaussianProposal(n_dof, sigma)
    if method in ABLATIONS:
        opts, kind = ABLATIONS[method]
        perlin = PerlinFieldProposal(
            dof_coords, sigma, corr_length=corr_length, seed=seed, name=method, n_steps=n_steps, **opts
        )
        if kind == "perlin":
            return perlin
        if kind == "gauss":
            return CorrelatedGaussianProposal.matched_to(perlin, name=method)
        if kind == "gauss-perlin-marginals":
            return CorrelatedGaussianProposal.matched_to(perlin, name=method, marginal_quantiles=perlin._sorted)
        if kind == "gauss-ar1":
            return CorrelatedGaussianProposal.matched_to(perlin, name=method, ar=perlin.lag1)
        if kind == "gauss-anneal":
            stages = [
                PerlinFieldProposal(dof_coords, sigma, corr_length=L, seed=seed).correlation()
                for L in perlin.stage_lengths()
            ]
            return CorrelatedGaussianProposal.matched_to(perlin, name=method, schedule=stages, n_steps=n_steps)
    raise ValueError(f"unknown proposal method {method!r}")
