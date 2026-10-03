"""Ground-truth targets must keep the input's stereocentres and trans secondary amides."""

import numpy as np

from perlin_mole.chem import TorsionSystem
from perlin_mole.validate import ALANINE_DIPEPTIDE, _dihedral, keeps_fixed_features, reference_targets


def _heavy(smiles: str) -> np.ndarray:
    s = TorsionSystem.from_smiles(smiles)
    xyz, _, _ = s.minimize(s.initial_xyz)
    return xyz[s.heavy_idx]


def test_input_geometry_is_kept_and_mirror_image_is_rejected():
    x = _heavy(ALANINE_DIPEPTIDE)
    assert keeps_fixed_features(ALANINE_DIPEPTIDE, x)
    mirror = x * np.array([-1.0, 1.0, 1.0])  # inverts the stereocentre, keeps the amides trans
    assert not keeps_fixed_features(ALANINE_DIPEPTIDE, mirror)


def test_cis_amide_is_rejected():
    from rdkit import Chem

    from perlin_mole.validate import _AMIDE

    x = _heavy(ALANINE_DIPEPTIDE)
    mol = Chem.MolFromSmiles(ALANINE_DIPEPTIDE)
    o, c, n, r = mol.GetSubstructMatches(Chem.MolFromSmarts(_AMIDE))[0]
    assert abs(_dihedral(*x[[o, c, n, r]])) < 90.0
    # rotate the N-substituent side of the amide by 180 degrees about the C-N bond: trans -> cis
    axis = (x[n] - x[c]) / np.linalg.norm(x[n] - x[c])
    side = [r]
    th = np.pi
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    rot = np.eye(3) + np.sin(th) * k + (1 - np.cos(th)) * k @ k
    y = x.copy()
    y[side] = (x[side] - x[n]) @ rot.T + x[n]
    assert abs(_dihedral(*y[[o, c, n, r]])) > 90.0
    assert not keeps_fixed_features(ALANINE_DIPEPTIDE, y)


def test_reference_targets_skips_excluded_minima_and_rebases_window():
    x = _heavy(ALANINE_DIPEPTIDE)
    mirror = x * np.array([-1.0, 1.0, 1.0])
    energies = np.array([0.0, 1.0, 5.0])
    heavy = np.stack([mirror, x, x])  # the lowest minimum is the wrong stereoisomer
    idx, e0, excluded = reference_targets(energies, heavy, ALANINE_DIPEPTIDE, window=3.0)
    assert list(idx) == [1]
    assert e0 == 1.0
    assert excluded == 1
