"""The five models. Each is a rule for the real k-by-k block K_ij on an edge.

K[i, j] is what node i receives from node j. A complex number p + iq is the
block [[p, -q], [q, p]], so complex models are two reals per node.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import networkx as nx
import numpy as np

from .graph import Graph

KINDS = ("M0", "M1", "M2", "M3H", "M3S")


def rot(theta: float) -> np.ndarray:
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, -s], [s, c]])


def orientation_signs(g: Graph, convention: str, dest: int) -> dict[tuple[int, int], int]:
    """sigma_ij for every edge key (i, j) with i < j. sigma_ji = -sigma_ij."""
    edges = g.edges
    if convention == "id":
        return {e: 1 for e in edges}
    if convention == "reverse":
        return {e: -1 for e in edges}
    if convention.startswith("random:"):
        rng = np.random.default_rng(int(convention.split(":")[1]))
        return {e: int(s) for e, s in zip(edges, rng.choice([-1, 1], size=len(edges)))}
    if convention == "oracle_hops":
        # Diagnostic only: uses global knowledge of the destination.
        hops = nx.single_source_shortest_path_length(g.to_networkx(), dest)
        return {(i, j): -1 if hops[i] < hops[j] else 1 for i, j in edges}
    raise ValueError(f"unknown orientation convention {convention!r}")


@dataclass(frozen=True)
class Model:
    kind: str
    kappa: float = 0.0  # M0, M1: attenuation scale
    kappa2: float = 0.0  # M1: second attenuation scale
    beta: float = 0.0  # M3H, M3S: phase per unit cost
    theta0: float = -math.pi / 2  # M3S: per-hop phase offset
    orientation: str = "id"  # M3H: edge sign convention
    a: float = 1.0

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise ValueError(f"unknown model kind {self.kind!r}")

    @property
    def k(self) -> int:
        return 1 if self.kind == "M0" else 2

    @property
    def is_complex(self) -> bool:
        return self.kind in ("M2", "M3H", "M3S")

    @property
    def variant(self) -> str:
        if self.kind == "M3H":
            return f"M3H[{self.orientation}]"
        if self.kind == "M3S":
            return f"M3S[theta0={'0' if self.theta0 == 0 else '-pi/2'}]"
        return self.kind

    @property
    def param(self) -> float:
        if self.kind in ("M0", "M1"):
            return self.kappa
        return self.beta if self.kind in ("M3H", "M3S") else 0.0

    def block(self, cost: float, sigma: int) -> np.ndarray:
        """K_ij for one direction of an edge. sigma is the M3-H sign for that direction."""
        if self.kind == "M0":
            return np.array([[self.a * math.exp(-self.kappa * cost)]])
        if self.kind == "M1":
            return self.a * np.diag([math.exp(-self.kappa * cost), math.exp(-self.kappa2 * cost)])
        if self.kind == "M2":
            return self.a * rot(-math.pi / 2)
        if self.kind == "M3H":
            return self.a * rot(-math.pi / 2 + sigma * self.beta * cost)
        return self.a * rot(self.theta0 + self.beta * cost)

    def coupling(self, g: Graph, dest: int) -> np.ndarray:
        K = np.zeros((g.n, g.n, self.k, self.k))
        signs = orientation_signs(g, self.orientation, dest) if self.kind == "M3H" else {}
        for i, j in g.edges:
            sigma = signs.get((i, j), 1)
            K[i, j] = self.block(g.cost(i, j), sigma)
            K[j, i] = self.block(g.cost(i, j), -sigma)
        return K

    def source(self, n: int, dest: int) -> np.ndarray:
        b = np.zeros((n, self.k))
        b[dest] = (1.0, 0.0) if self.is_complex else 1.0
        return b


def to_complex(x: np.ndarray) -> np.ndarray:
    """(n, 2) real state -> complex vector."""
    return x[:, 0] + 1j * x[:, 1]


def complex_matrix(K: np.ndarray) -> np.ndarray:
    """(n, n, 2, 2) rotation-scaling blocks -> complex (n, n) matrix."""
    return K[:, :, 0, 0] + 1j * K[:, :, 1, 0]
