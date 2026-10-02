"""Undirected weighted graph with positive symmetric link costs."""
from __future__ import annotations

from dataclasses import dataclass

import networkx as nx


@dataclass(frozen=True)
class Graph:
    n: int
    costs: dict[tuple[int, int], float]  # key (i, j) with i < j
    names: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for (i, j), c in self.costs.items():
            if not 0 <= i < j < self.n:
                raise ValueError(f"edge key must satisfy 0 <= i < j < n, got {(i, j)}")
            if not c > 0:
                raise ValueError(f"cost must be positive, got {c} on {(i, j)}")
        if self.names and len(self.names) != self.n:
            raise ValueError("names must have one entry per node")

    @property
    def edges(self) -> list[tuple[int, int]]:
        return sorted(self.costs)

    def cost(self, i: int, j: int) -> float:
        return self.costs[(min(i, j), max(i, j))]

    def neighbors(self, i: int) -> list[int]:
        return sorted(b if a == i else a for a, b in self.costs if i in (a, b))

    def with_cost(self, i: int, j: int, cost: float) -> Graph:
        key = (min(i, j), max(i, j))
        if key not in self.costs:
            raise KeyError(f"no edge {key}")
        return Graph(self.n, {**self.costs, key: cost}, self.names)

    def name(self, i: int) -> str:
        return self.names[i] if self.names else str(i)

    def to_networkx(self) -> nx.Graph:
        g = nx.Graph()
        g.add_nodes_from(range(self.n))
        g.add_weighted_edges_from((i, j, c) for (i, j), c in self.costs.items())
        return g
