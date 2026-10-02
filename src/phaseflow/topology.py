"""Topologies for the smoke test. All costs are 1 unless changed afterwards."""
from __future__ import annotations

import numpy as np

from .graph import Graph


def _unit(n: int, pairs: list[tuple[int, int]], names: str = "") -> Graph:
    return Graph(n, {(min(i, j), max(i, j)): 1.0 for i, j in pairs}, tuple(names))


def path(n: int) -> Graph:
    return _unit(n, [(i, i + 1) for i in range(n - 1)])


def ring(n: int) -> Graph:
    return _unit(n, [(i, (i + 1) % n) for i in range(n)])


def binary_tree(depth: int) -> Graph:
    n = 2 ** (depth + 1) - 1
    return _unit(n, [((child - 1) // 2, child) for child in range(1, n)])


def diamond() -> Graph:
    """S-A, then arms A-B-D-T (upper) and A-C-E-T (lower). One 6-cycle."""
    return _unit(7, [(0, 1), (1, 2), (2, 4), (4, 6), (1, 3), (3, 5), (5, 6)], "SABCDET")


def lopsided_diamond() -> Graph:
    """S-A, then arms A-B-T (upper) and A-C-E-T (lower). One 5-cycle."""
    return _unit(6, [(0, 1), (1, 2), (2, 5), (1, 3), (3, 4), (4, 5)], "SABCET")


def random_costs(g: Graph, seed: int, lo: float = 0.5, hi: float = 2.0) -> Graph:
    rng = np.random.default_rng(seed)
    return Graph(g.n, {e: float(rng.uniform(lo, hi)) for e in g.edges}, g.names)
