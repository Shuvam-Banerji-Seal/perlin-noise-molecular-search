import pickle

import numpy as np
from rdkit import Chem

from perlin_mole.benchmark import Task, run_task
from perlin_mole.export import analyze
from perlin_mole.validate import classify_ala, systematic_search


def test_classify_alanine_basins():
    assert classify_ala(-83, 73) == "C7eq"
    assert classify_ala(-157, 165) == "C5"
    assert classify_ala(74, -60) == "C7ax"
    assert classify_ala(0, 180) == "other"


def test_systematic_search_butane_finds_anti_and_two_gauche():
    e, heavy, full, n = systematic_search("CCCC", [0.0, 60.0, 120.0, 180.0, 240.0, 300.0], jobs=2)
    assert n == 6
    assert len(e) == 3  # anti + gauche(+) + gauche(-), mirror images kept distinct
    assert np.all(np.diff(e) >= 0)
    assert abs(e[1] - e[2]) < 1e-3 and e[1] - e[0] > 0.3  # degenerate gauche pair above anti


def test_analyze_writes_sdf_and_pools(tmp_path):
    runs = [run_task(Task("hexane", "CCCCCC", m, seed=0, n_steps=15)) for m in ("perlin", "gauss-iid")]
    (tmp_path / "raw").mkdir()
    with open(tmp_path / "raw" / "raw.pkl", "wb") as fh:
        pickle.dump(runs, fh)
    rep = analyze([tmp_path / "raw"], tmp_path / "out", systematic_dir=None, jobs=1)
    r = rep["molecules"]["hexane"]
    assert r["n_low"] >= 1 and set(r["coverage_common"]) == {"perlin(L=4) @60", "gauss-iid @60"}
    mols = [m for m in Chem.SDMolSupplier(str(tmp_path / "out" / "sdf" / "hexane_low_energy.sdf"), removeHs=False)]
    assert len(mols) == r["n_low"]
    assert all(m.GetNumAtoms() == 20 for m in mols)  # hydrogens rebuilt
    assert r["h_rebuild_max_abs_dE"] < 0.05


def test_mirror_and_stereo_utilities():
    from perlin_mole.chem import TorsionSystem
    from perlin_mole.isomers import cip_labels, mirror_match, mirror_smiles, stereo_fidelity
    from perlin_mole.search import SymmetryRMSD

    s = "CC(C)NC[C@@H](O)COc1cccc2ccccc12"
    m = mirror_smiles(s)
    assert cip_labels(s) == {5: "R"} and cip_labels(m) == {5: "S"}
    t = TorsionSystem.from_smiles(s)
    x = t.coords()[t.heavy_idx]
    assert stereo_fidelity(s, [x]) == 1.0 and stereo_fidelity(s, [x * [-1, 1, 1]]) == 0.0
    rm = SymmetryRMSD(t.heavy_mol())
    assert mirror_match([(1.0, x)], [(1.0, x * [-1, 1, 1])], rm) == 1.0
    assert mirror_match([(1.0, x)], [(1.0, x)], rm) == 0.0  # a chiral conformer is not its own mirror image


def test_cov_mat_metrics():
    from perlin_mole.ablation import _cov_mat
    from perlin_mole.chem import TorsionSystem
    from perlin_mole.search import SymmetryRMSD

    t = TorsionSystem.from_smiles("CCCCCC")
    x = t.coords()[t.heavy_idx]
    rm = SymmetryRMSD(t.heavy_mol())
    far = x + np.random.default_rng(0).normal(0, 1.5, x.shape)
    m = _cov_mat(np.stack([x, far]), np.stack([x]), rm)
    assert m["cov_r_05"] == 0.5 and m["cov_p_05"] == 1.0 and m["mat_p"] < 1e-6 and m["mat_r"] > 0.2
    empty = _cov_mat(np.stack([x]), np.empty((0, *x.shape)), rm)
    assert empty["cov_r_05"] == 0.0 and empty["mat_p"] != empty["mat_p"]  # NaN when nothing generated


def test_mirror_match_survives_atom_reordering():
    from perlin_mole.chem import TorsionSystem
    from perlin_mole.isomers import atom_order_map, mirror_match, mirror_smiles
    from perlin_mole.search import SymmetryRMSD

    a = "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"  # ibuprofen: mirrored SMILES reorders atoms
    b = mirror_smiles(a)
    perm = atom_order_map(a, b)
    assert not np.array_equal(perm, np.arange(len(perm)))
    ta = TorsionSystem.from_smiles(a)
    xa = ta.coords()[ta.heavy_idx]
    xb = np.empty_like(xa)
    xb[perm] = xa * [-1, 1, 1]  # exact mirror image, written in b's atom order
    rm = SymmetryRMSD(ta.heavy_mol())
    assert mirror_match([(1.0, xa)], [(1.0, xb[perm])], rm) == 1.0
