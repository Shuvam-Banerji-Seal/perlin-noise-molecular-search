import numpy as np
import pytest

from perlin_mole.benchmark import Task, run_task
from perlin_mole.chem import TorsionSystem, classical_mds
from perlin_mole.landscape import TorsionPerlinBias, gradient_correlation
from perlin_mole.proposals import (
    CorrelatedGaussianProposal,
    PerlinFieldProposal,
    RandomResetProposal,
    build_proposal,
)
from perlin_mole.search import ConformerArchive, SymmetryRMSD, basin_hopping


@pytest.fixture(scope="module")
def pentane():
    return TorsionSystem.from_smiles("CCCCCC", seed=0)


def test_rotatable_torsions_exclude_terminal_and_amide():
    assert TorsionSystem.from_smiles("CCCCCC").n_dof == 3  # hexane: 3 internal C-C
    assert TorsionSystem.from_smiles("CC(=O)NC").n_dof == 0  # N-methylacetamide
    assert TorsionSystem.from_smiles("c1ccccc1").n_dof == 0


def test_apply_delta_is_exact(pentane):
    x = pentane.coords()
    d = np.array([0.3, -1.2, 2.0])
    before = pentane.torsion_values(x)
    after = pentane.torsion_values(pentane.apply_delta(x, d))
    assert np.allclose(np.angle(np.exp(1j * (after - before - d))), 0.0, atol=1e-10)


def test_minimize_lowers_energy(pentane):
    x = pentane.apply_delta(pentane.coords(), np.array([2.0, 2.0, 2.0]))
    e0 = pentane.energy(x)
    _, e1, conv = pentane.minimize(x)
    assert conv and e1 <= e0


def test_mds_recovers_chain():
    d = np.abs(np.subtract.outer(np.arange(5.0), np.arange(5.0)))
    x = classical_mds(d, 3)
    rec = np.linalg.norm(x[:, None] - x[None], axis=-1)
    assert np.allclose(rec, d, atol=1e-8)


def test_symmetry_rmsd_end_to_end_flip(pentane):
    rmsd = SymmetryRMSD(pentane.heavy_mol())
    xyz = pentane.coords()[pentane.heavy_idx]
    flipped = xyz[::-1]  # hexane is symmetric under reversing atom order
    assert rmsd(xyz, flipped) < 1e-6
    archive = ConformerArchive(rmsd)
    assert archive.add(xyz, 1.0, 0)
    assert not archive.add(flipped, 1.02, 1)
    assert archive.add(xyz + 0.0, 2.0, 2)  # same geometry, different energy => new entry


def _symmetric(samples):
    return np.abs(samples.mean(0)).max() < 0.05 and np.allclose(np.median(samples, 0), 0.0, atol=0.08)


def test_proposals_are_symmetric_and_scaled():
    coords = np.c_[np.arange(6.0), np.zeros(6), np.zeros(6)]
    rng = np.random.default_rng(0)
    perlin = PerlinFieldProposal(coords, sigma=1.0, corr_length=4.0)
    gauss = CorrelatedGaussianProposal.matched_to(perlin)
    for prop in (perlin, gauss):
        s = np.array([prop.sample(rng) for _ in range(4000)])
        assert _symmetric(s)
        assert np.allclose(s.std(0), 1.0, atol=0.08)
    # neighbouring torsions move together, distant ones do not
    c = perlin.correlation()
    assert c[0, 1] > 0.5 and abs(c[0, 5]) < 0.2
    s = np.array([gauss.sample(rng) for _ in range(4000)])
    assert np.allclose(np.corrcoef(s.T), c, atol=0.06)


def test_reset_changes_few_torsions():
    rng = np.random.default_rng(0)
    p = RandomResetProposal(10, max_changed=3)
    counts = [np.count_nonzero(p.sample(rng)) for _ in range(200)]
    assert min(counts) >= 1 and max(counts) <= 3


def test_full_dimensional_gradient_is_not_collective():
    c = gradient_correlation(4, n_points=8000)
    assert np.abs(c[~np.eye(4, dtype=bool)]).max() < 0.1


def test_torsion_bias_periodic():
    b = TorsionPerlinBias(3, strength=2.0, seed=1)
    phi = np.random.default_rng(0).uniform(-np.pi, np.pi, (50, 3))
    assert np.allclose(b(phi), b(phi + 2 * np.pi), atol=1e-10)


@pytest.mark.parametrize("method", ["reset", "gauss-iid", "gauss-corr", "perlin"])
def test_basin_hopping_runs(pentane, method):
    prop = build_proposal(method, pentane.dof_coordinates(3), np.deg2rad(60), 4.0)
    archive = ConformerArchive(SymmetryRMSD(pentane.heavy_mol()))
    res = basin_hopping(pentane, prop, 20, 1.0, np.random.default_rng(0), pentane.coords(), archive)
    assert len(res.best_trace) == 21
    assert np.all(np.diff(res.best_trace) <= 1e-12)
    assert 1 <= len(res.minima) <= 21


def test_run_task_etkdg_is_deterministic():
    t = Task("hexane", "CCCCCC", "etkdg", seed=1, n_steps=5)
    a, b = run_task(t), run_task(t)
    assert a["best_trace"] == b["best_trace"]


def test_parallel_deduplicate_matches_sequential_archive(pentane):
    from perlin_mole.search import deduplicate

    rng = np.random.default_rng(0)
    base = pentane.coords()[pentane.heavy_idx]
    xyz, es = [], []
    for k in range(60):  # clusters of near-identical geometries with near-identical energies
        center = base + rng.normal(0, 0.6, base.shape)
        for _ in range(3):
            xyz.append(center + rng.normal(0, 0.03, base.shape))
            es.append(k * 0.05 + rng.normal(0, 0.01))
    order = np.argsort(es, kind="stable")
    es, xyz = np.asarray(es)[order], np.asarray(xyz)[order]
    rmsd = SymmetryRMSD(pentane.heavy_mol())

    archive = ConformerArchive(rmsd)
    seq_rep = []
    for e, x in zip(es, xyz, strict=True):
        i = archive.find(x, e)
        if i < 0:
            archive.add(x, e, 0)
            i = len(archive.minima) - 1
        seq_rep.append(float(archive.minima[i].energy))

    for jobs in (1, 3):
        rep = deduplicate(es, xyz, rmsd, jobs=jobs)
        assert np.allclose(es[rep], seq_rep)
    assert 40 <= int((rep == np.arange(len(rep))).sum()) < len(es)
