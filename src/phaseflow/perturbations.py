"""Deterministic link-cost change schedule."""
from __future__ import annotations

from dataclasses import dataclass

from .graph import Graph


@dataclass(frozen=True)
class CostChange:
    round: int
    edge: tuple[int, int]
    cost: float


def phases(g: Graph, changes, total_rounds: int) -> list[tuple[int, int, Graph]]:
    """Split [0, total_rounds) at each change. Returns (start, end, graph in force)."""
    out, start, current = [], 0, g
    for change in sorted(changes, key=lambda c: c.round):
        if not start < change.round < total_rounds:
            raise ValueError(f"change round {change.round} must lie strictly inside ({start}, {total_rounds})")
        out.append((start, change.round, current))
        current = current.with_cost(*change.edge, change.cost)
        start = change.round
    out.append((start, total_rounds, current))
    return out
