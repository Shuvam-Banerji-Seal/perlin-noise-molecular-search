import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation

from perlin_mole.experimental import (
    attach_coords,
    best_alignment,
    canonical_with_coords,
    load_structure,
    relax_experimental,
    run_metrics,
    write_structure,
)
from perlin_mole.search import SymmetryRMSD

CHIRAL = "C[C@@H](C(=O)O)c1ccc(CC(C)C)cc1"  # (S)-ibuprofen


def _embedded_heavy(smiles: str, seed: int = 7) -> Chem.Mol:
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=seed)
    return Chem.RemoveHs(m)


def _as_instance(heavy: Chem.Mol, rng: np.random.Generator) -> Chem.Mol:
    """Shuffle atom order like a PDB/COD file would (bond orders kept, coordinates moved along)."""
    perm = [int(i) for i in rng.permutation(heavy.GetNumAtoms())]
    return Chem.RenumberAtoms(heavy, perm)


def test_mapping_round_trip_gives_zero_rmsd(tmp_path):
    heavy = _embedded_heavy(CHIRAL)
    inst = _as_instance(heavy, np.random.default_rng(0))
    smiles, mol = canonical_with_coords(Chem.MolFromSmiles(CHIRAL), inst)
    path = write_structure(tmp_path / "x.sdf", mol, "x")
    xyz = load_structure(path, smiles)
    # same molecule, original coordinates: RMSD must be ~0 after mapping into SMILES order
    ref = Chem.MolFromSmiles(smiles)
    order = heavy.GetSubstructMatch(ref, useChirality=True)
    original = heavy.GetConformer().GetPositions()[list(order)]
    # SDF stores 4 decimals, so "zero" means below 1e-4 Angstrom
    assert SymmetryRMSD(ref)(original, xyz) < 1e-4
    assert np.abs(original - xyz).max() < 1e-4  # identical frame, not just superposable


def test_wrong_enantiomer_is_rejected_and_unspecified_stereo_is_read_from_3d():
    heavy = _embedded_heavy(CHIRAL)
    mirrored = attach_coords(heavy, heavy.GetConformer().GetPositions() * np.array([1.0, 1.0, -1.0]))
    template = Chem.MolFromSmiles(CHIRAL)
    assert canonical_with_coords(template, heavy) is not None
    assert canonical_with_coords(template, mirrored) is None  # (R) structure vs (S) template
    flat = Chem.MolFromSmiles("CC(C(=O)O)c1ccc(CC(C)C)cc1")  # racemic template
    s1, _ = canonical_with_coords(flat, heavy)
    s2, _ = canonical_with_coords(flat, mirrored)
    assert "@" in s1 and "@" in s2 and s1 != s2


def test_load_structure_rejects_atom_order_mismatch(tmp_path):
    heavy = _embedded_heavy("CCOC(=O)c1ccccc1")
    path = write_structure(tmp_path / "y.sdf", heavy, "y")
    with pytest.raises(ValueError):
        load_structure(path, "O=C(OCC)c1ccccc1")


def test_run_metrics_on_toy_run():
    heavy = _embedded_heavy(CHIRAL)
    exp = heavy.GetConformer().GetPositions()
    rng = np.random.default_rng(1)
    rot = Rotation.random(random_state=2).as_matrix()
    minima = [
        (5.0, exp + rng.normal(0, 1.2, exp.shape), 0),  # far
        (2.0, exp @ rot.T + 3.0, 7),  # exact match, rotated + translated, found at step 7
        (1.0, exp + rng.normal(0, 0.9, exp.shape), 9),  # lowest energy, not a match
    ]
    m = run_metrics({"minima": minima}, exp, SymmetryRMSD(heavy), n_steps=20)
    assert m["min_rmsd"] < 1e-6 and m["success_1"] and m["success_05"]
    assert m["t_1"] == 7 and m["match_rank"] == 2 and abs(m["match_dE"] - 1.0) < 1e-12
    assert m["rmsd_best"] > 0.5


def test_best_alignment_recovers_permuted_rotated_copy():
    heavy = _embedded_heavy("CCCCc1ccccc1")
    exp = heavy.GetConformer().GetPositions()
    perms = SymmetryRMSD(heavy).perms
    rot = Rotation.random(random_state=3).as_matrix()
    pred = (exp @ rot.T)[perms[-1]]  # rotated and relabelled by an automorphism
    aligned, r = best_alignment(exp, pred, perms)
    assert r < 1e-6 and np.abs(aligned - exp).max() < 1e-5


def test_relax_experimental_is_close_for_an_mmff_minimum():
    m = Chem.AddHs(Chem.MolFromSmiles("CCCCO"))
    AllChem.EmbedMolecule(m, randomSeed=1)
    AllChem.MMFFOptimizeMolecule(m, maxIters=5000)
    heavy = Chem.RemoveHs(m).GetConformer().GetPositions()
    out = relax_experimental("CCCCO", heavy)
    assert out["rmsd_exp_to_own_minimum"] < 0.1
    assert out["e_exp_frozen"] >= out["e_exp_relaxed"] - 1e-6
