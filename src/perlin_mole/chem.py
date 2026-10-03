"""Molecule handling: rotatable torsions as degrees of freedom, MMFF94 as the judge.

The molecular graph is fixed here (Level 1: conformers). Proposals act on
torsion angles; every proposed geometry is relaxed with MMFF94 so the search
compares local minima, never raw perturbed geometries.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from rdkit.Chem import rdForceFieldHelpers as ffh

# Acyclic single bonds between two non-terminal heavy atoms, no triple bonds.
_ROTATABLE = Chem.MolFromSmarts("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]")
# Amide C-N bonds are near-planar with a high barrier; flipping them changes
# cis/trans isomer class, so they are excluded as search DOFs.
_AMIDE = Chem.MolFromSmarts("[CX3](=O)-&!@[NX3]")


def rotatable_torsions(heavy: Chem.Mol) -> list[tuple[int, int, int, int]]:
    """Return one (i, j, k, l) dihedral per rotatable bond j-k of a heavy-atom mol."""
    amide = {frozenset((m[0], m[2])) for m in heavy.GetSubstructMatches(_AMIDE)}
    out = []
    for j, k in heavy.GetSubstructMatches(_ROTATABLE):
        if frozenset((j, k)) in amide:
            continue
        aj, ak = heavy.GetAtomWithIdx(j), heavy.GetAtomWithIdx(k)
        i = min(n.GetIdx() for n in aj.GetNeighbors() if n.GetIdx() != k)
        l = min(n.GetIdx() for n in ak.GetNeighbors() if n.GetIdx() != j)
        out.append((i, j, k, l))
    return sorted(out, key=lambda t: (t[1], t[2]))


def bond_graph_distance(heavy: Chem.Mol, torsions: list[tuple[int, int, int, int]]) -> np.ndarray:
    """Topological distance between torsion bonds, in bonds (adjacent bonds -> 1)."""
    dm = Chem.GetDistanceMatrix(heavy)
    n = len(torsions)
    out = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):
            _, j1, k1, _ = torsions[a]
            _, j2, k2, _ = torsions[b]
            d = min(dm[j1, j2], dm[j1, k2], dm[k1, j2], dm[k1, k2]) + 1.0
            out[a, b] = out[b, a] = d
    return out


def classical_mds(dist: np.ndarray, dim: int) -> np.ndarray:
    """Embed a distance matrix in ``dim`` Euclidean dimensions (zero-padded if rank < dim)."""
    n = dist.shape[0]
    if n == 1:
        return np.zeros((1, dim))
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ (dist**2) @ j
    w, v = np.linalg.eigh(b)
    order = np.argsort(w)[::-1][:dim]
    w, v = np.clip(w[order], 0.0, None), v[:, order]
    x = v * np.sqrt(w)
    if x.shape[1] < dim:
        x = np.hstack([x, np.zeros((n, dim - x.shape[1]))])
    return x


@dataclass
class TorsionSystem:
    """A molecule whose search space is its rotatable torsions."""

    smiles: str
    mol: Chem.Mol  # with explicit H, exactly one conformer (the working buffer)
    heavy_idx: np.ndarray
    torsions: list[tuple[int, int, int, int]]
    props: object = field(repr=False)
    max_iters: int = 2000
    # Immutable copy of the construction-time embedding. ``coords()`` reads the
    # mutable working buffer, which every minimize()/apply_delta() overwrites.
    initial_xyz: np.ndarray | None = field(default=None, repr=False)

    @classmethod
    def from_smiles(cls, smiles: str, seed: int = 0, max_iters: int = 2000) -> TorsionSystem:
        heavy = Chem.MolFromSmiles(smiles)
        if heavy is None:
            raise ValueError(f"could not parse SMILES {smiles!r}")
        torsions = rotatable_torsions(heavy)
        mol = Chem.AddHs(heavy)  # AddHs keeps heavy-atom indices unchanged
        params = AllChem.ETKDGv3()
        params.randomSeed = int(seed)
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise RuntimeError(f"ETKDG embedding failed for {smiles!r}")
        props = ffh.MMFFGetMoleculeProperties(mol)
        if props is None:
            raise ValueError(f"MMFF94 has no parameters for {smiles!r}")
        heavy_idx = np.array([a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1])
        initial = mol.GetConformer().GetPositions().copy()
        return cls(smiles, mol, heavy_idx, torsions, props, max_iters, initial)

    @property
    def n_dof(self) -> int:
        return len(self.torsions)

    @property
    def n_atoms(self) -> int:
        return self.mol.GetNumAtoms()

    def heavy_mol(self) -> Chem.Mol:
        return Chem.RemoveHs(self.mol)

    def coords(self) -> np.ndarray:
        """Current contents of the working buffer (changes after every minimize/apply_delta)."""
        return self.mol.GetConformer().GetPositions()

    def _load(self, xyz: np.ndarray) -> Chem.Conformer:
        conf = self.mol.GetConformer()
        for a, p in enumerate(np.asarray(xyz, dtype=np.float64)):
            conf.SetAtomPosition(a, p.tolist())
        return conf

    def torsion_values(self, xyz: np.ndarray) -> np.ndarray:
        conf = self._load(xyz)
        return np.array([rdMolTransforms.GetDihedralRad(conf, *t) for t in self.torsions])

    def apply_delta(self, xyz: np.ndarray, delta: np.ndarray) -> np.ndarray:
        """Rotate every torsion by ``delta`` (radians). Rotations commute exactly for acyclic bonds."""
        conf = self._load(xyz)
        for t, d in zip(self.torsions, delta, strict=True):
            if d != 0.0:
                rdMolTransforms.SetDihedralRad(conf, *t, rdMolTransforms.GetDihedralRad(conf, *t) + float(d))
        return conf.GetPositions()

    def energy(self, xyz: np.ndarray) -> float:
        self._load(xyz)
        return ffh.MMFFGetMoleculeForceField(self.mol, self.props).CalcEnergy()

    def minimize(self, xyz: np.ndarray) -> tuple[np.ndarray, float, bool]:
        """MMFF94 local relaxation. Returns (coords, energy kcal/mol, converged)."""
        self._load(xyz)
        ff = ffh.MMFFGetMoleculeForceField(self.mol, self.props)
        status = ff.Minimize(maxIts=self.max_iters)
        return self.coords(), ff.CalcEnergy(), status == 0

    def embed_random(self, seed: int) -> np.ndarray:
        """A fresh ETKDG geometry (the standard independent-sampling baseline)."""
        work = Chem.Mol(self.mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = int(seed)
        if AllChem.EmbedMolecule(work, params) != 0:
            params.useRandomCoords = True
            if AllChem.EmbedMolecule(work, params) != 0:
                raise RuntimeError("embedding failed")
        return work.GetConformer().GetPositions()

    def dof_coordinates(self, dim: int = 3) -> np.ndarray:
        """Place each torsion in R^dim so Euclidean distance ~ topological distance (bonds)."""
        return classical_mds(bond_graph_distance(self.heavy_mol(), self.torsions), dim)
