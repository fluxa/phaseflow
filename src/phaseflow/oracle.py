"""Exact solutions used only to validate the local dynamics.

Nothing here may be passed to the engine as an input.
"""
from __future__ import annotations

import numpy as np


def flat(K: np.ndarray) -> np.ndarray:
    n, _, k, _ = K.shape
    return K.transpose(0, 2, 1, 3).reshape(n * k, n * k)


def fixed_point(K: np.ndarray, b: np.ndarray, gamma: float) -> np.ndarray:
    n, _, k, _ = K.shape
    return np.linalg.solve(gamma * np.eye(n * k) - flat(K), b.reshape(n * k)).reshape(n, k)


def abscissa(K: np.ndarray, gamma: float) -> float:
    """Largest real part of the spectrum of -gamma*I + K. Negative means stable."""
    return float(np.linalg.eigvals(flat(K)).real.max() - gamma)


def iteration_radius(K: np.ndarray, gamma: float, alpha: np.ndarray) -> float:
    """Spectral radius of the local iteration matrix. Below 1 means convergent."""
    n, _, k, _ = K.shape
    M = flat(K) - gamma * np.eye(n * k)
    G = np.eye(n * k) + np.repeat(alpha, k)[:, None] * M
    return float(np.abs(np.linalg.eigvals(G)).max())
