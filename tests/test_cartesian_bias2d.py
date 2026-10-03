import numpy as np

from perlin_mole.benchmark import _system
from perlin_mole.bias2d import MatchedGRFBias, PerlinBias, Surface, exact_dF, run_mc
from perlin_mole.cartesian import _heavy_graph, intact
from perlin_mole.isomers import cip_labels


def test_intact_detects_bond_change_and_inversion():
    smi = "CC(C)NC[C@@H](O)COc1cccc2ccccc12"
    s = _system(smi)
    x = s.coords()
    g, cip = _heavy_graph(s.mol), cip_labels(smi)
    assert intact(s, x, g, cip)
    assert not intact(s, x * np.array([-1.0, 1.0, 1.0]), g, cip)  # mirror image inverts the stereocentre
    y = x.copy()
    y[0] = y[12] + 0.05  # crush a methyl carbon onto a ring atom -> new "bond"
    assert not intact(s, y, g, cip)


def test_bias_fields_periodic_and_normalised():
    pb = PerlinBias(1.0, cells=4, seed=1)
    p = np.random.default_rng(0).uniform(-180, 180, 500)
    q = np.random.default_rng(1).uniform(-180, 180, 500)
    assert np.allclose(pb(p, q), pb(p + 360.0, q - 360.0), atol=1e-9)
    v = pb(np.random.default_rng(2).uniform(-180, 180, 20000), np.random.default_rng(3).uniform(-180, 180, 20000))
    assert 0.8 < v.std() < 1.2
    g = MatchedGRFBias(pb, n=48, seed=0, n_spec=8)
    w = g(np.random.default_rng(4).uniform(-180, 180, 5000), np.random.default_rng(5).uniform(-180, 180, 5000))
    assert 0.7 < w.std() < 1.3


def test_mc_on_flat_surface_is_unbiased():
    flat = Surface(np.zeros((36, 36)), 10.0)
    assert abs(exact_dF(flat)) < 1e-9
    r = run_mc(flat, None, n_steps=4000, n_walkers=16, step_deg=30.0, seed=0)
    assert r["transitions"].sum() > 0
    assert abs(np.nanmean(r["dF"])) < 0.25  # equal populations of phi>0 and phi<0
