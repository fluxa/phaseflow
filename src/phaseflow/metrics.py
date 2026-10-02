"""Route quality metrics. Shortest paths come from networkx as an independent reference."""
from __future__ import annotations

import math

import networkx as nx
import numpy as np

from .graph import Graph
from .readout import next_hops, policy

DARK = 1 - 1e-6


def shortest(g: Graph, dest: int) -> tuple[np.ndarray, list[set[int]]]:
    """Cost-to-destination per node, and the set of optimal next hops per node."""
    lengths = nx.single_source_dijkstra_path_length(g.to_networkx(), dest)
    dist = np.array([lengths[i] for i in range(g.n)])
    opt = [
        {j for j in g.neighbors(i) if math.isclose(g.cost(i, j) + dist[j], dist[i], rel_tol=1e-9)}
        for i in range(g.n)
    ]
    opt[dest] = set()
    return dist, opt


def follow(nxt: np.ndarray, src: int, dest: int) -> tuple[str, list[int]]:
    """Follow deterministic next hops from src. Status: delivered, loop or dead."""
    route, seen = [src], {src}
    while route[-1] != dest:
        j = int(nxt[route[-1]])
        if j < 0:
            return "dead", route
        if j in seen:
            return "loop", route + [j]
        route.append(j)
        seen.add(j)
    return "delivered", route


def route_metrics(g: Graph, dest: int, W: np.ndarray, canc: np.ndarray) -> dict:
    dist, opt = shortest(g, dest)
    nxt, ties = next_hops(W, dest)
    nodes = [i for i in range(g.n) if i != dest]
    statuses, stretches, margins = [], [], []
    for i in nodes:
        status, route = follow(nxt, i, dest)
        statuses.append(status)
        if status == "delivered":
            stretches.append(sum(g.cost(u, v) for u, v in zip(route, route[1:])) / dist[i])
        if len(g.neighbors(i)) >= 2:
            total = W[i].sum()
            top = np.sort(W[i])[::-1]
            margins.append(float((top[0] - top[1]) / total) if total > 0 else 0.0)
    c = canc[nodes]
    seen = c[~np.isnan(c)]
    return {
        "nexthop_acc": float(np.mean([int(nxt[i]) in opt[i] for i in nodes])),
        "delivered_frac": statuses.count("delivered") / len(nodes),
        "n_loop_nodes": statuses.count("loop"),
        "n_dead_nodes": statuses.count("dead"),
        "n_ties": int(ties[nodes].sum()),
        "n_dark_nodes": int((seen > DARK).sum()),
        "mean_stretch": float(np.mean(stretches)) if stretches else math.nan,
        "max_stretch": float(np.max(stretches)) if stretches else math.nan,
        "mean_cancellation": float(seen.mean()) if len(seen) else math.nan,
        "max_cancellation": float(seen.max()) if len(seen) else math.nan,
        "min_margin": min(margins) if margins else math.nan,
    }


def stochastic_metrics(g: Graph, dest: int, W: np.ndarray) -> dict:
    """Delivery probability, expected-cost stretch and entropy of p(i -> j) ∝ W[i, j]."""
    dist, _ = shortest(g, dest)
    Pm = policy(W, dest)
    nodes = [i for i in range(g.n) if i != dest]
    Q, r = Pm[np.ix_(nodes, nodes)], Pm[nodes, dest]
    eye = np.eye(len(nodes))
    if np.abs(np.linalg.eigvals(Q)).max() < 1 - 1e-12:
        h = np.linalg.solve(eye - Q, r)
    else:  # a closed class that excludes the destination: take the minimal solution
        h = np.zeros(len(nodes))
        for _ in range(20000):
            h = r + Q @ h
    if h.min() > 1 - 1e-9:
        step_cost = np.array([sum(Pm[i, j] * g.cost(i, j) for j in g.neighbors(i)) for i in nodes])
        exp_stretch = float(np.mean(np.linalg.solve(eye - Q, step_cost) / dist[nodes]))
    else:
        exp_stretch = math.inf
    logs = np.log(Pm, out=np.zeros_like(Pm), where=Pm > 0)
    entropy = -(Pm * logs).sum(axis=1)
    return {
        "p_delivery_min": float(h.min()),
        "exp_stretch_mean": exp_stretch,
        "entropy_mean": float(entropy[nodes].mean()),
    }


def circulation(P: np.ndarray, cycle: tuple[int, ...]) -> float:
    """Mean antisymmetric flow along the cycle in the order given."""
    m = len(cycle)
    flows = [(P[cycle[(t + 1) % m], cycle[t]] - P[cycle[t], cycle[(t + 1) % m]]) / 2 for t in range(m)]
    return float(np.mean(flows))


def split_share(W: np.ndarray, node: int, upper: int, lower: int) -> float:
    """Share of node's two-arm forwarding weight placed on the upper arm."""
    total = W[node, upper] + W[node, lower]
    return float(W[node, upper] / total) if total > 0 else math.nan
