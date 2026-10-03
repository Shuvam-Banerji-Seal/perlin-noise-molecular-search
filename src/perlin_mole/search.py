"""Basin hopping in torsion space and a symmetry-aware archive of minima."""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem

from perlin_mole.chem import TorsionSystem


def _kabsch_rmsd_batch(p: np.ndarray, qs: np.ndarray) -> np.ndarray:
    """RMSD of p (n,3) to each qs[m] (n,3) after optimal proper rotation."""
    p = p - p.mean(axis=0)
    qs = qs - qs.mean(axis=1, keepdims=True)
    h = np.einsum("ni,mnj->mij", p, qs)
    u, s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(np.einsum("mij,mjk->mik", u, vt)))
    s[:, -1] *= d
    e = (p**2).sum() + (qs**2).sum(axis=(1, 2)) - 2.0 * s.sum(axis=1)
    return np.sqrt(np.clip(e, 0.0, None) / p.shape[0])


class SymmetryRMSD:
    """Heavy-atom RMSD minimised over the graph automorphisms of the molecule."""

    def __init__(self, heavy: Chem.Mol, max_matches: int = 1000) -> None:
        matches = heavy.GetSubstructMatches(heavy, uniquify=False, maxMatches=max_matches)
        self.perms = np.array(matches, dtype=np.int64)

    def __call__(self, a: np.ndarray, b: np.ndarray) -> float:
        return float(_kabsch_rmsd_batch(a, b[self.perms]).min())

    def many(self, a: np.ndarray, bs: np.ndarray) -> np.ndarray:
        """Symmetry-minimised RMSD of ``a`` to each ``bs[c]``, in one batched SVD."""
        c, m = bs.shape[0], self.perms.shape[0]
        qs = bs[:, self.perms].reshape(c * m, *a.shape)
        return _kabsch_rmsd_batch(a, qs).reshape(c, m).min(axis=1)


@dataclass
class Minimum:
    energy: float
    heavy_xyz: np.ndarray
    step: int


@dataclass
class ConformerArchive:
    """Unique local minima. Two minima are the same if |dE| < e_tol and RMSD < rmsd_tol."""

    rmsd: SymmetryRMSD
    e_tol: float = 0.1
    rmsd_tol: float = 0.25
    minima: list[Minimum] = field(default_factory=list)
    _e: np.ndarray = field(default_factory=lambda: np.empty(64), repr=False)
    _xyz: np.ndarray | None = field(default=None, repr=False)

    def find(self, heavy_xyz: np.ndarray, energy: float) -> int:
        n = len(self.minima)
        if n == 0:
            return -1
        cand = np.flatnonzero(np.abs(self._e[:n] - energy) < self.e_tol)
        if cand.size == 0:
            return -1
        hits = np.flatnonzero(self.rmsd.many(heavy_xyz, self._xyz[cand]) < self.rmsd_tol)
        return int(cand[hits[0]]) if hits.size else -1

    def add(self, heavy_xyz: np.ndarray, energy: float, step: int) -> bool:
        if self.find(heavy_xyz, energy) >= 0:
            return False
        n = len(self.minima)
        if self._xyz is None:
            self._xyz = np.empty((64, *np.shape(heavy_xyz)))
        if n == self._e.shape[0]:
            self._e = np.concatenate([self._e, np.empty(n)])
            self._xyz = np.concatenate([self._xyz, np.empty_like(self._xyz)])
        self._e[n] = energy
        self._xyz[n] = heavy_xyz
        self.minima.append(Minimum(energy, np.array(heavy_xyz), step))
        return True

    @property
    def energies(self) -> np.ndarray:
        return self._e[: len(self.minima)].copy()


@dataclass
class RunResult:
    method: str
    seed: int
    energies: np.ndarray  # relaxed energy at every step (index 0 = start)
    best_trace: np.ndarray
    minima: list[Minimum]
    acceptance: float
    unconverged: int
    wall_time: float
    rejected_chem: int = 0


def basin_hopping(
    system: TorsionSystem,
    proposal,
    n_steps: int,
    kT: float,
    rng: np.random.Generator,
    start_xyz: np.ndarray,
    archive: ConformerArchive,
    guard=None,
) -> RunResult:
    """Metropolis basin hopping: kick torsions -> MMFF relax -> accept on relaxed energy.

    ``guard(xyz) -> bool`` (optional) rejects relaxed structures that fail a chemistry check, e.g. an
    inverted stereocentre: relaxing a badly clashed geometry can push a substituent through the centre.
    """
    t0 = time.perf_counter()
    heavy = system.heavy_idx
    xyz, e, conv = system.minimize(start_xyz)
    unconverged = int(not conv)
    archive.add(xyz[heavy], e, 0)
    energies, best = [e], [e]
    accepted = 0
    rejected_chem = 0
    for step in range(1, n_steps + 1):
        trial = system.apply_delta(xyz, proposal.sample(rng))
        trial, e_t, conv = system.minimize(trial)
        unconverged += int(not conv)
        if guard is not None and not guard(trial):
            rejected_chem += 1
            energies.append(e_t)
            best.append(best[-1])
            continue
        archive.add(trial[heavy], e_t, step)
        if e_t <= e or rng.random() < np.exp(-(e_t - e) / kT):
            xyz, e = trial, e_t
            accepted += 1
        energies.append(e_t)
        best.append(min(best[-1], e_t))
    return RunResult(
        getattr(proposal, "name", type(proposal).__name__),
        -1,
        np.array(energies),
        np.array(best),
        archive.minima,
        accepted / max(n_steps, 1),
        unconverged,
        time.perf_counter() - t0,
        rejected_chem,
    )


def random_embedding(
    system: TorsionSystem,
    n_steps: int,
    seed: int,
    start_xyz: np.ndarray,
    archive: ConformerArchive,
) -> RunResult:
    """Independent ETKDG embeddings + MMFF relax: the standard conformer-generation baseline."""
    t0 = time.perf_counter()
    heavy = system.heavy_idx
    xyz, e, conv = system.minimize(start_xyz)
    unconverged = int(not conv)
    archive.add(xyz[heavy], e, 0)
    energies, best = [e], [e]
    for step in range(1, n_steps + 1):
        trial, e_t, conv = system.minimize(system.embed_random(seed * 1_000_003 + step))
        unconverged += int(not conv)
        archive.add(trial[heavy], e_t, step)
        energies.append(e_t)
        best.append(min(best[-1], e_t))
    return RunResult(
        "etkdg",
        -1,
        np.array(energies),
        np.array(best),
        archive.minima,
        float("nan"),
        unconverged,
        time.perf_counter() - t0,
    )


# --- parallel offline de-duplication -------------------------------------------------

_DEDUP: dict = {}


def _dedup_init(energies, xyz, perms, e_tol, rmsd_tol) -> None:
    _DEDUP.update(e=energies, x=xyz, perms=perms, e_tol=e_tol, rmsd_tol=rmsd_tol)


def _dedup_chunk(idx: np.ndarray) -> list[np.ndarray]:
    e, x, perms = _DEDUP["e"], _DEDUP["x"], _DEDUP["perms"]
    out = []
    for i in idx:
        lo = int(np.searchsorted(e, e[i] - _DEDUP["e_tol"], side="right"))
        if lo >= i:
            out.append(np.empty(0, dtype=np.int64))
            continue
        c, m = i - lo, perms.shape[0]
        qs = x[lo:i][:, perms].reshape(c * m, *x.shape[1:])
        r = _kabsch_rmsd_batch(x[i], qs).reshape(c, m).min(axis=1)
        out.append(lo + np.flatnonzero(r < _DEDUP["rmsd_tol"]))
    return out


def deduplicate(
    energies: np.ndarray,
    xyz: np.ndarray,
    rmsd: SymmetryRMSD,
    e_tol: float = 0.1,
    rmsd_tol: float = 0.25,
    jobs: int = 1,
) -> np.ndarray:
    """Cluster minima; returns ``rep[i]`` = index of the unique minimum that ``i`` duplicates.

    ``energies`` must be sorted ascending. Same rule as ``ConformerArchive``
    (|dE| < e_tol and symmetry RMSD < rmsd_tol), applied greedily in energy
    order. The expensive part, all pairwise RMSDs inside each energy window,
    runs on ``jobs`` processes. The greedy resolution is a cheap linear pass.
    """
    n = len(energies)
    if np.any(np.diff(energies) < 0):
        raise ValueError("energies must be sorted ascending")
    chunks = [np.arange(k, n, max(1, jobs * 16)) for k in range(min(n, max(1, jobs * 16)))]
    args = (energies, xyz, rmsd.perms, e_tol, rmsd_tol)
    matches: list[np.ndarray] = [None] * n  # type: ignore[list-item]
    if jobs <= 1:
        _dedup_init(*args)
        results = map(_dedup_chunk, chunks)
    else:
        from concurrent.futures import ProcessPoolExecutor

        pool = ProcessPoolExecutor(max_workers=jobs, initializer=_dedup_init, initargs=args)
        results = pool.map(_dedup_chunk, chunks)
    for idx, res in zip(chunks, results, strict=True):
        for i, mt in zip(idx, res, strict=True):
            matches[i] = mt
    if jobs > 1:
        pool.shutdown()

    rep = np.arange(n)
    for i in range(n):
        for j in matches[i]:
            if rep[j] == j:
                rep[i] = j
                break
    return rep
