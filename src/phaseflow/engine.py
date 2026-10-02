"""Local synchronous iteration.

One round: every node sends its k reals to each neighbour, then updates from
what it received. Node i uses only its own state, its own blocks K[i, :] and
its neighbours' states. Nothing from the oracle enters here.
"""
from __future__ import annotations

import numpy as np


def step_sizes(K: np.ndarray, gamma: float) -> np.ndarray:
    """alpha_i = gamma / (gamma^2 + s_i^2), s_i = sum_j ||K_ij||_2 (node-local)."""
    s = np.linalg.norm(K, ord=2, axis=(2, 3)).sum(axis=1)
    return gamma / (gamma**2 + s**2)


def step(x: np.ndarray, K: np.ndarray, b: np.ndarray, gamma: float, alpha: np.ndarray) -> np.ndarray:
    received = np.einsum("ijab,jb->ia", K, x)
    return x + alpha[:, None] * (-gamma * x + received + b)


def run(
    x: np.ndarray, K: np.ndarray, b: np.ndarray, gamma: float, rounds: int, diverge_at: float = 1e12
) -> tuple[np.ndarray, int, bool]:
    """Iterate; stop early and report divergence if the state leaves [-diverge_at, diverge_at]."""
    alpha = step_sizes(K, gamma)
    for t in range(rounds):
        x = step(x, K, b, gamma, alpha)
        if not np.isfinite(x).all() or np.abs(x).max() > diverge_at:
            return x, t + 1, True
    return x, rounds, False
