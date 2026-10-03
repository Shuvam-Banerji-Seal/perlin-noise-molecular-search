"""Periodic n-dimensional Perlin gradient noise with analytic gradients.

The lattice is periodic with an integer period per axis, so the field is a
smooth function on a torus. That matters for torsion angles, which live on
T^k rather than R^k.

Cost per evaluated point is O(2^n * n): every corner of the enclosing
hypercube contributes. That is cheap for n <= ~8 and impractical beyond
~15, which is the main practical limit on putting a Perlin landscape
directly over all torsions of a molecule (see docs/feasibility.md).
"""

from __future__ import annotations

import itertools

import numpy as np


def _fade(t: np.ndarray) -> np.ndarray:
    """Quintic smoothstep 6t^5 - 15t^4 + 10t^3 (C2 continuous)."""
    return t * t * t * (t * (t * 6.0 - 15.0) + 10.0)


def _dfade(t: np.ndarray) -> np.ndarray:
    return 30.0 * t * t * (t * (t - 2.0) + 1.0)


def _leave_one_out_prod(a: np.ndarray) -> np.ndarray:
    """out[..., k] = prod_{j != k} a[..., j], without dividing."""
    pre = np.ones_like(a)
    suf = np.ones_like(a)
    pre[..., 1:] = np.cumprod(a[..., :-1], axis=-1)
    suf[..., :-1] = np.cumprod(a[..., ::-1], axis=-1)[..., ::-1][..., 1:]
    return pre * suf


class PerlinNoise:
    """Single-octave gradient noise on an n-torus.

    Parameters
    ----------
    ndim:
        Dimension of the input space.
    period:
        Integer lattice period (scalar or one per axis). N(x + period*e_d) == N(x).
    seed:
        Seed for the gradient table and permutation.
    table_size:
        Number of distinct random gradients / permutation entries.
    """

    def __init__(
        self,
        ndim: int,
        period: int | np.ndarray = 256,
        seed: int | None = None,
        table_size: int = 256,
    ) -> None:
        if ndim < 1:
            raise ValueError("ndim must be >= 1")
        rng = np.random.default_rng(seed)
        self.ndim = ndim
        self.period = np.broadcast_to(np.asarray(period, dtype=np.int64), (ndim,)).copy()
        if np.any(self.period < 1):
            raise ValueError("period must be >= 1")
        g = rng.standard_normal((table_size, ndim))
        if ndim == 1:
            # Unit "vectors" in 1-D are just +-1; use magnitudes too so the
            # field is not trivially symmetric.
            g = rng.uniform(-1.0, 1.0, (table_size, 1))
        else:
            g /= np.linalg.norm(g, axis=1, keepdims=True)
        self.gradients = g
        self.perm = rng.permutation(table_size).astype(np.int64)
        self.table_size = table_size
        self.corners = np.array(list(itertools.product((0, 1), repeat=ndim)), dtype=np.int64)

    def _hash(self, idx: np.ndarray) -> np.ndarray:
        h = np.zeros(idx.shape[:-1], dtype=np.int64)
        for d in range(self.ndim):
            h = self.perm[(h + idx[..., d]) % self.table_size]
        return h

    def evaluate(self, x: np.ndarray, grad: bool = False):
        """Evaluate noise (and optionally its gradient) at points ``x[..., ndim]``."""
        x = np.asarray(x, dtype=np.float64)
        if x.shape[-1] != self.ndim:
            raise ValueError(f"expected last axis of size {self.ndim}, got {x.shape}")
        shape = x.shape[:-1]
        x = x.reshape(-1, self.ndim)

        i0 = np.floor(x).astype(np.int64)
        f = x - i0
        u = _fade(f)
        du = _dfade(f) if grad else None

        val = np.zeros(x.shape[0])
        g_out = np.zeros_like(x) if grad else None
        for c in self.corners:
            gc = self.gradients[self._hash((i0 + c) % self.period)]
            dot = np.einsum("md,md->m", gc, f - c)
            w_d = np.where(c == 1, u, 1.0 - u)
            w = np.prod(w_d, axis=1)
            val += w * dot
            if grad:
                sign = np.where(c == 1, 1.0, -1.0)
                dw = _leave_one_out_prod(w_d) * sign * du
                g_out += dw * dot[:, None] + w[:, None] * gc

        val = val.reshape(shape)
        if grad:
            return val, g_out.reshape(shape + (self.ndim,))
        return val

    __call__ = evaluate


class FractalNoise:
    """Fractional Brownian motion: a sum of Perlin octaves.

    Octave ``o`` samples ``N_o(x * lacunarity**o)`` with amplitude
    ``persistence**o``. Each octave's lattice period is scaled so the sum
    stays periodic with the base ``period``.
    """

    def __init__(
        self,
        ndim: int,
        octaves: int = 3,
        persistence: float = 0.5,
        lacunarity: int = 2,
        period: int = 64,
        seed: int | None = None,
    ) -> None:
        if octaves < 1:
            raise ValueError("octaves must be >= 1")
        if int(lacunarity) != lacunarity or lacunarity < 1:
            raise ValueError("lacunarity must be a positive integer to keep periodicity")
        seeds = np.random.SeedSequence(seed).spawn(octaves)
        self.ndim = ndim
        self.period = int(period)
        self.octaves = [
            PerlinNoise(ndim, period=self.period * lacunarity**o, seed=s.generate_state(1)[0])
            for o, s in enumerate(seeds)
        ]
        self.amplitudes = np.array([persistence**o for o in range(octaves)])
        self.frequencies = np.array([float(lacunarity**o) for o in range(octaves)])

    def evaluate(self, x: np.ndarray, grad: bool = False):
        x = np.asarray(x, dtype=np.float64)
        val = np.zeros(x.shape[:-1])
        g_out = np.zeros(x.shape) if grad else None
        for noise, a, fr in zip(self.octaves, self.amplitudes, self.frequencies, strict=True):
            if grad:
                v, g = noise.evaluate(x * fr, grad=True)
                val += a * v
                g_out += a * fr * g
            else:
                val += a * noise.evaluate(x * fr)
        return (val, g_out) if grad else val

    __call__ = evaluate
