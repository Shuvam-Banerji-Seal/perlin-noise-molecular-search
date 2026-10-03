import numpy as np
import pytest

from perlin_mole.clusters import (
    KNOWN_GM,
    GaussianFieldDisplacement,
    IIDDisplacement,
    PerlinDisplacement,
    basin_hopping,
    lj_energy,
    lj_energy_grad,
    perlin_correlation_function,
    relax,
)


def _icosahedron13() -> np.ndarray:
    phi = (1 + 5**0.5) / 2
    v = []
    for a in (-1, 1):
        for b in (-phi, phi):
            v += [(0, a, b), (a, b, 0), (b, 0, a)]
    v = np.array(v, float)
    v /= np.linalg.norm(v[0])
    return np.vstack([np.zeros(3), 1.1 * v])


def test_gradient_matches_finite_difference():
    x = _icosahedron13().ravel() + np.random.default_rng(0).normal(0, 0.05, 39)
    _, g = lj_energy_grad(x)
    h = 1e-6
    fd = np.array([(lj_energy(x + h * e) - lj_energy(x - h * e)) / (2 * h) for e in np.eye(x.size)])
    assert np.allclose(g, fd, atol=1e-5)


def test_dimer_minimum():
    x, e = relax(np.array([[0, 0, 0], [1.3, 0, 0]], float))
    assert abs(e + 1.0) < 1e-8
    assert abs(np.linalg.norm(x[0] - x[1]) - 2 ** (1 / 6)) < 1e-5


def test_lj13_icosahedron_reproduces_known_gm():
    x = _icosahedron13() + np.random.default_rng(1).normal(0, 0.02, (13, 3))
    _, e = relax(x)
    assert abs(e - KNOWN_GM[13]["energy"]) < 1e-5


@pytest.mark.parametrize(
    "prop", [IIDDisplacement(0.3), PerlinDisplacement(0.3, 2.0), GaussianFieldDisplacement(0.3, 2.0)]
)
def test_proposals_have_marginal_std(prop):
    rng = np.random.default_rng(0)
    x = _icosahedron13()
    d = np.concatenate([prop.sample(x, rng) for _ in range(600)])
    assert abs(d.std() - 0.3) < 0.03
    assert abs(d.mean()) < 0.02


def test_matched_field_has_perlin_correlation():
    rng = np.random.default_rng(1)
    x = _icosahedron13()
    L = 2.0
    p, g = PerlinDisplacement(0.3, L, seed=3), GaussianFieldDisplacement(0.3, L)
    dp = np.array([p.sample(x, rng)[:, 0] for _ in range(4000)])
    dg = np.array([g.sample(x, rng)[:, 0] for _ in range(4000)])
    cp, cg = np.corrcoef(dp.T), np.corrcoef(dg.T)
    # centre atom vs a shell atom (distance 1.1 sigma) and two shell atoms
    assert cp[0, 1] > 0.1
    assert abs(cp[0, 1] - cg[0, 1]) < 0.12
    assert np.abs(cp - cg).mean() < 0.08


def test_correlation_function_shape():
    r, k = perlin_correlation_function(n=3000)
    assert abs(k[0] - 1) < 1e-9
    assert np.interp(0.25, r, k) > np.interp(0.75, r, k)
    assert abs(np.interp(3.0, r, k)) < 0.1


def test_basin_hopping_finds_lj13_quickly():
    out = basin_hopping(13, "gauss-iid", seed=0, budget=60, s=0.3)
    assert out["t_hit"] is not None and out["gap"] < 1e-4
