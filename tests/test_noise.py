import numpy as np
import pytest

from perlin_mole.noise import FractalNoise, PerlinNoise


@pytest.mark.parametrize("ndim", [1, 2, 3, 5])
def test_zero_on_lattice(ndim):
    n = PerlinNoise(ndim, period=16, seed=1)
    pts = np.random.default_rng(0).integers(0, 16, size=(50, ndim)).astype(float)
    assert np.allclose(n(pts), 0.0, atol=1e-12)


@pytest.mark.parametrize("ndim", [1, 2, 3, 4])
def test_periodic(ndim):
    period = 8
    n = PerlinNoise(ndim, period=period, seed=2)
    x = np.random.default_rng(1).uniform(0, period, size=(200, ndim))
    for d in range(ndim):
        shift = np.zeros(ndim)
        shift[d] = period
        assert np.allclose(n(x), n(x + shift), atol=1e-12)


@pytest.mark.parametrize("ndim", [1, 2, 3, 4, 6])
def test_gradient_matches_finite_difference(ndim):
    n = FractalNoise(ndim, octaves=3, period=8, seed=3)
    x = np.random.default_rng(2).uniform(0, 8, size=(40, ndim))
    _, g = n(x, grad=True)
    h = 1e-6
    fd = np.empty_like(g)
    for d in range(ndim):
        e = np.zeros(ndim)
        e[d] = h
        fd[:, d] = (n(x + e) - n(x - e)) / (2 * h)
    assert np.allclose(g, fd, atol=1e-5)


def test_fractal_periodic_and_bounded():
    n = FractalNoise(3, octaves=4, period=4, seed=4)
    x = np.random.default_rng(3).uniform(0, 4, size=(2000, 3))
    v = n(x)
    assert np.allclose(v, n(x + np.array([4.0, 0, 0])), atol=1e-12)
    assert np.all(np.abs(v) < 2.0)
    assert v.std() > 0.05


def test_smoothness_correlation_decays_with_distance():
    n = PerlinNoise(2, period=64, seed=5)
    x = np.random.default_rng(4).uniform(0, 64, size=(20000, 2))
    v0 = n(x)
    near = np.corrcoef(v0, n(x + [0.1, 0.0]))[0, 1]
    far = np.corrcoef(v0, n(x + [3.0, 0.0]))[0, 1]
    assert near > 0.9
    assert abs(far) < 0.1
