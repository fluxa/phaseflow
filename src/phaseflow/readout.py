"""Forwarding readout shared by every model.

P[i, j] = 2 x_i^T K_ij x_j is the rate at which neighbour j feeds |x_i|^2.
For Hermitian models it is the probability current into i from j.
"""
from __future__ import annotations

import numpy as np


def power_inflow(K: np.ndarray, x: np.ndarray) -> np.ndarray:
    return 2.0 * np.einsum("ia,ijab,jb->ij", x, K, x)


def weights(P: np.ndarray) -> np.ndarray:
    return np.maximum(P, 0.0)


def cancellation(K: np.ndarray, x: np.ndarray, is_complex: bool) -> np.ndarray:
    """1 - |sum of neighbour contributions| / sum of |contributions|, per node.

    A complex state is one channel (modulus over both components); a real
    state is one channel per component. NaN where nothing is received.
    """
    contrib = np.einsum("ijab,jb->ija", K, x)
    if is_complex:
        total = np.linalg.norm(contrib.sum(axis=1), axis=1)
        parts = np.linalg.norm(contrib, axis=2).sum(axis=1)
    else:
        total = np.abs(contrib.sum(axis=1)).sum(axis=1)
        parts = np.abs(contrib).sum(axis=(1, 2))
    out = np.full(len(x), np.nan)
    np.divide(total, parts, out=out, where=parts > 0)
    return 1.0 - out


def next_hops(W: np.ndarray, dest: int, rtol: float = 1e-9) -> tuple[np.ndarray, np.ndarray]:
    """Deterministic next hop per node (-1 if no positive weight) and a tie flag."""
    n = W.shape[0]
    nxt = np.full(n, -1, dtype=int)
    ties = np.zeros(n, dtype=bool)
    nxt[dest] = dest
    for i in range(n):
        if i == dest:
            continue
        best = W[i].max()
        if best <= 0:
            continue
        top = np.flatnonzero(W[i] >= best * (1 - rtol))
        nxt[i] = top[0]
        ties[i] = len(top) > 1
    return nxt, ties


def policy(W: np.ndarray, dest: int) -> np.ndarray:
    """Stochastic forwarding p(i -> j) proportional to W[i, j]."""
    Pm = W.copy()
    Pm[dest] = 0.0
    total = Pm.sum(axis=1, keepdims=True)
    return np.divide(Pm, total, out=np.zeros_like(Pm), where=total > 0)
