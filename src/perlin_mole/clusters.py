"""Lennard-Jones clusters: can Perlin-correlated moves find *known* global minima faster?

LJ_N clusters in reduced units (epsilon = sigma = 1) are the standard test bed
for structure prediction. Their global minima (GM) are tabulated in the
Cambridge Cluster Database; LJ38 is the fcc truncated octahedron, a small
crystal fragment. That gives a "predicted vs real" test with known answers.

Basin hopping (Wales & Doye 1997) moves all atoms, relaxes with L-BFGS, and
accepts on the relaxed energy. Three displacement proposals share the same
marginal std ``s`` per Cartesian component:

* ``gauss-iid``: independent N(0, s^2) for every coordinate.
* ``perlin``: three independent fBm Perlin fields (one per Cartesian
  component), each randomly rotated, offset and signed, read at the *current*
  atom positions scaled by 1/L. Reading the field at current positions makes
  the proposal state-dependent. That is fine here, because basin hopping is
  used as an optimiser, not a sampler.
* ``gauss-field``: a Gaussian random field whose covariance
  K_ij = s^2 k(|r_i - r_j| / L) uses the Perlin field's *empirical* isotropic
  correlation function k. It matches Perlin's second-order structure exactly
  at the current geometry, which isolates any Perlin-specific effect.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.stats import special_ortho_group

from perlin_mole.noise import FractalNoise

DATA = Path(__file__).resolve().parents[2] / "results" / "clusters"
CCD_URL = "https://www-wales.ch.cam.ac.uk/~jon/structures/LJ/"

# Global-minimum energies (reduced units) and point groups, Cambridge Cluster Database
# table of LJ_N global minima (tables.150.html), checked 2026-10-02.
KNOWN_GM = {
    13: {"energy": -44.326801, "point_group": "Ih"},
    19: {"energy": -72.659782, "point_group": "D5h"},
    26: {"energy": -108.315616, "point_group": "Td"},
    38: {"energy": -173.928427, "point_group": "Oh"},
}


# ----------------------------------------------------------------------------- potential
def lj_energy_grad(x: np.ndarray) -> tuple[float, np.ndarray]:
    """Total LJ energy and gradient for flat coordinates x (3N,)."""
    p = x.reshape(-1, 3)
    d = p[:, None, :] - p[None, :, :]
    r2 = np.einsum("ijk,ijk->ij", d, d)
    np.fill_diagonal(r2, np.inf)
    ir2 = 1.0 / r2
    ir6 = ir2 * ir2 * ir2
    ir12 = ir6 * ir6
    e = 2.0 * np.sum(ir12 - ir6)  # 4 * sum_{i<j} = 2 * sum_{i != j}
    coef = (-48.0 * ir12 + 24.0 * ir6) * ir2  # dE/dr * 1/r for each pair
    g = np.einsum("ij,ijk->ik", coef, d)
    return float(e), g.ravel()


def lj_energy(x: np.ndarray) -> float:
    return lj_energy_grad(x)[0]


def relax(x: np.ndarray, gtol: float = 1e-6, maxiter: int = 2000) -> tuple[np.ndarray, float]:
    res = minimize(
        lj_energy_grad,
        np.asarray(x, float).ravel(),
        jac=True,
        method="L-BFGS-B",
        options={"gtol": gtol, "ftol": 1e-14, "maxiter": maxiter},
    )
    return res.x.reshape(-1, 3), float(res.fun)


def load_known(n: int) -> np.ndarray | None:
    p = DATA / f"gm_LJ{n}.xyz.txt"
    if not p.exists():
        return None
    return np.loadtxt(p).reshape(-1, 3)


def random_start(n: int, rng: np.random.Generator) -> np.ndarray:
    """Uniform random points in a sphere at roughly liquid density, then relaxed."""
    radius = 1.1 * (3 * n / (4 * np.pi)) ** (1 / 3)
    pts = []
    while len(pts) < n:
        p = rng.uniform(-radius, radius, 3)
        if np.linalg.norm(p) <= radius and all(np.linalg.norm(p - q) > 0.7 for q in pts):
            pts.append(p)
    x, _ = relax(np.array(pts))
    return x - x.mean(axis=0)


def containment_radius(n: int) -> float:
    return (3 * n / (4 * np.pi)) ** (1 / 3) + 2.0


# ----------------------------------------------------------------------------- Perlin correlation function
@lru_cache(maxsize=8)
def perlin_correlation_function(
    octaves: int = 2, seed: int = 0, r_max: float = 6.0, dr: float = 0.05, n: int = 20000
) -> tuple[np.ndarray, np.ndarray]:
    """Isotropic correlation k(r) of the fBm field, averaged over random positions and directions."""
    f = FractalNoise(3, octaves=octaves, period=64, seed=seed)
    rng = np.random.default_rng(seed + 11)
    y = rng.uniform(0, 64, size=(n, 3))
    u = rng.standard_normal((n, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    r = np.arange(0.0, r_max + dr / 2, dr)
    f0 = f(y)
    f0c = f0 - f0.mean()
    k = np.empty_like(r)
    for i, ri in enumerate(r):
        fi = f(y + ri * u)
        fic = fi - fi.mean()
        k[i] = np.sum(f0c * fic) / np.sqrt(np.sum(f0c**2) * np.sum(fic**2))
    return r, k


# ----------------------------------------------------------------------------- proposals
@dataclass
class IIDDisplacement:
    s: float
    name: str = "gauss-iid"

    def sample(self, x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        return rng.normal(0.0, self.s, x.shape)


@dataclass
class PerlinDisplacement:
    """Displacement d_i^c = s * sign_c * F_c(R x_i / L + o_c) / std(F), c = x, y, z.

    State-dependent: the field is read at the current positions x_i.
    """

    s: float
    corr_length: float = 1.0
    octaves: int = 2
    seed: int = 0
    period: int = 64
    name: str = "perlin"
    _fields: list = field(init=False, repr=False)
    _scale: float = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self._fields = [
            FractalNoise(3, octaves=self.octaves, period=self.period, seed=self.seed * 3 + c) for c in range(3)
        ]
        rng = np.random.default_rng(self.seed + 7)
        pts = rng.uniform(0, self.period, size=(20000, 3))
        self._scale = float(np.std(np.concatenate([f(pts) for f in self._fields])))

    def sample(self, x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        rot = special_ortho_group.rvs(3, random_state=rng)
        y = (x / self.corr_length) @ rot.T
        out = np.empty_like(x)
        for c, f in enumerate(self._fields):
            o = rng.uniform(0, self.period, 3)
            out[:, c] = rng.choice([-1.0, 1.0]) * f(y + o)
        return self.s * out / self._scale


@dataclass
class GaussianFieldDisplacement:
    """Gaussian random field with the Perlin field's isotropic correlation k(r/L), per component."""

    s: float
    corr_length: float = 1.0
    octaves: int = 2
    seed: int = 0
    jitter: float = 1e-6
    name: str = "gauss-field"

    def kernel(self, x: np.ndarray) -> np.ndarray:
        r, k = perlin_correlation_function(self.octaves, 0)
        dist = np.linalg.norm(x[:, None] - x[None], axis=-1) / self.corr_length
        return np.interp(dist, r, k, right=0.0)

    def sample(self, x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        K = self.kernel(x)
        w, v = np.linalg.eigh(K)
        root = v * np.sqrt(np.clip(w, 0.0, None))  # robust square root (K is PSD up to MC noise)
        z = rng.standard_normal((x.shape[0], 3))
        return self.s * (root @ z)


@dataclass
class IIDPerlinMarginalDisplacement:
    """Independent per-coordinate kicks drawn from the Perlin field's own (unit-variance) marginal distribution.

    Separates "Perlin's light-tailed marginal" from "Perlin's spatial structure": no correlation at all.
    """

    s: float
    seed: int = 0
    name: str = "iid-perlinmarg"
    _values: np.ndarray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        f = FractalNoise(3, octaves=2, period=64, seed=self.seed)
        v = f(np.random.default_rng(self.seed + 5).uniform(0, 64, size=(200000, 3)))
        self._values = v / v.std()

    def sample(self, x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        return self.s * rng.choice(self._values, size=x.shape) * rng.choice([-1.0, 1.0], size=x.shape)


def make_proposal(method: str, s: float, seed: int):
    if method == "gauss-iid":
        return IIDDisplacement(s)
    if method == "iid-perlinmarg":
        return IIDPerlinMarginalDisplacement(s, seed=seed)
    kind, L = method.split("@")
    L = float(L)
    if kind == "perlin":
        return PerlinDisplacement(s, corr_length=L, seed=seed, name=method)
    if kind == "gauss-field":
        return GaussianFieldDisplacement(s, corr_length=L, seed=seed, name=method)
    raise ValueError(method)


# ----------------------------------------------------------------------------- basin hopping
def basin_hopping(
    n: int, method: str, seed: int, budget: int, s: float = 0.35, kT: float = 0.8, tol: float = 1e-4
) -> dict:
    """Run BH from the seed's common start; returns time-to-GM (None if not reached) and traces."""
    start = random_start(n, np.random.default_rng(10_000 + seed))  # same start for every method
    rng = np.random.default_rng([seed, sum(map(ord, method))])
    prop = make_proposal(method, s, seed)
    e_gm = KNOWN_GM[n]["energy"]
    rc = containment_radius(n)
    x = start
    e = lj_energy(x.ravel())
    best_e, best_x = e, x.copy()
    t_hit = 0 if e < e_gm + tol else None
    acc = rejected_evap = 0
    trace = [best_e]
    t0 = time.perf_counter()
    for step in range(1, budget + 1):
        xt, et = relax(x + prop.sample(x, rng))
        xt = xt - xt.mean(axis=0)
        if np.max(np.linalg.norm(xt, axis=1)) > rc:
            rejected_evap += 1
            trace.append(best_e)
            continue
        if et < best_e:
            best_e, best_x = et, xt.copy()
        if t_hit is None and et < e_gm + tol:
            t_hit = step
        if et <= e or rng.random() < np.exp(-(et - e) / kT):
            x, e = xt, et
            acc += 1
        trace.append(best_e)
    return {
        "n": n,
        "method": method,
        "seed": seed,
        "budget": budget,
        "t_hit": t_hit,
        "best_e": best_e,
        "gap": best_e - e_gm,
        "best_x": best_x.tolist(),
        "acceptance": acc / budget,
        "evaporation_rejects": rejected_evap,
        "trace": trace[:: max(1, budget // 100)],
        "wall_s": time.perf_counter() - t0,
    }


def write_known_minima() -> dict:
    info = {
        "source": CCD_URL + "tables.150.html",
        "coordinates": {str(n): CCD_URL + f"points/{n}" for n in KNOWN_GM},
        "units": "reduced LJ units (epsilon = sigma = 1)",
        "retrieved": "2026-10-02",
        "minima": {},
    }
    for n, v in KNOWN_GM.items():
        xk = load_known(n)
        rec = dict(v)
        if xk is not None:
            rec["energy_of_downloaded_coords"] = lj_energy(xk.ravel())
            rec["energy_after_relax"] = relax(xk)[1]
        info["minima"][str(n)] = rec
    (DATA / "known_minima.json").write_text(json.dumps(info, indent=1))
    return info
