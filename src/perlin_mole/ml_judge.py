"""A second, more accurate judge on the GPU: the ANI-2x neural-network potential (TorchANI).

MMFF94 decides what "stable" means everywhere else in this package. Here every
discovered low-energy minimum is re-judged with ANI-2x: a single-point energy at
the MMFF geometry, then a batched FIRE relaxation on the GPU. We ask whether the
energy ordering, and in particular the identity of the global minimum, survives a
change of judge. ANI-2x covers H, C, N, O, F, S and Cl; molecules with other
elements (e.g. Br) are skipped.
"""

from __future__ import annotations

import numpy as np

HARTREE_TO_KCAL = 627.509474
ANI2X_ELEMENTS = {1, 6, 7, 8, 9, 16, 17}

_MODEL = {}


def load_model(device: str | None = None):
    import torch
    import torchani

    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    if device not in _MODEL:
        _MODEL[device] = torchani.models.ANI2x(periodic_table_index=True).to(device).double()
    return _MODEL[device], device


def supported(atomic_numbers) -> bool:
    return set(int(z) for z in atomic_numbers) <= ANI2X_ELEMENTS


def energies(z: np.ndarray, xyz: np.ndarray, device: str | None = None, batch: int = 256) -> np.ndarray:
    """ANI-2x energies (kcal/mol) for a batch of conformers of one molecule. xyz: (n_conf, n_atoms, 3) in Angstrom."""
    import torch

    model, device = load_model(device)
    out = []
    for i in range(0, len(xyz), batch):
        c = torch.tensor(xyz[i : i + batch], dtype=torch.float64, device=device)
        s = torch.tensor(np.broadcast_to(z, c.shape[:2]).copy(), dtype=torch.long, device=device)
        with torch.no_grad():
            out.append(model((s, c)).energies.cpu().numpy())
    return np.concatenate(out) * HARTREE_TO_KCAL


def fire_relax(
    z: np.ndarray,
    xyz: np.ndarray,
    device: str | None = None,
    fmax: float = 0.05,
    max_steps: int = 2000,
    dt0: float = 0.05,
    dt_max: float = 0.5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Batched FIRE minimisation on the GPU. fmax in kcal/mol/Angstrom.

    Returns (relaxed xyz, energies in kcal/mol, converged flags).
    """
    import torch

    model, device = load_model(device)
    x = torch.tensor(xyz, dtype=torch.float64, device=device).requires_grad_(True)
    s = torch.tensor(np.broadcast_to(z, xyz.shape[:2]).copy(), dtype=torch.long, device=device)
    v = torch.zeros_like(x)
    nb = x.shape[0]
    dt = torch.full((nb, 1, 1), dt0, dtype=torch.float64, device=device)
    alpha = torch.full((nb, 1, 1), 0.1, dtype=torch.float64, device=device)
    n_pos = torch.zeros(nb, dtype=torch.long, device=device)
    done = torch.zeros(nb, dtype=torch.bool, device=device)
    for _ in range(max_steps):
        e = model((s, x)).energies * HARTREE_TO_KCAL
        (g,) = torch.autograd.grad(e.sum(), x)
        f = -g
        fmag = f.norm(dim=-1).amax(dim=-1)
        done = fmag < fmax
        if bool(done.all()):
            break
        with torch.no_grad():
            p = (f * v).sum(dim=(1, 2))
            fn = f / (f.norm(dim=(1, 2), keepdim=True) + 1e-12)
            vn = v.norm(dim=(1, 2), keepdim=True)
            v = (1 - alpha) * v + alpha * vn * fn
            up = p > 0
            n_pos = torch.where(up, n_pos + 1, torch.zeros_like(n_pos))
            grow = up & (n_pos > 5)
            dt = torch.where(grow[:, None, None], torch.clamp(dt * 1.1, max=dt_max), dt)
            alpha = torch.where(grow[:, None, None], alpha * 0.99, alpha)
            down = ~up
            dt = torch.where(down[:, None, None], dt * 0.5, dt)
            alpha = torch.where(down[:, None, None], torch.full_like(alpha, 0.1), alpha)
            v = torch.where(down[:, None, None], torch.zeros_like(v), v)
            v = v + dt * f
            step = dt * v
            step_len = step.norm(dim=-1, keepdim=True)
            step = torch.where(step_len > 0.2, step * 0.2 / step_len, step)  # cap atomic step at 0.2 A
            step = torch.where(done[:, None, None], torch.zeros_like(step), step)
            x += step
    with torch.no_grad():
        e = model((s, x)).energies * HARTREE_TO_KCAL
    return x.detach().cpu().numpy(), e.cpu().numpy(), done.cpu().numpy()
