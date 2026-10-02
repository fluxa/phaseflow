import networkx as nx
import pytest

from phaseflow import topology
from phaseflow.graph import Graph


def test_graph_rejects_bad_edges_and_costs():
    with pytest.raises(ValueError):
        Graph(2, {(1, 0): 1.0})
    with pytest.raises(ValueError):
        Graph(2, {(0, 1): 0.0})
    with pytest.raises(ValueError):
        Graph(2, {(0, 1): float("nan")})


def test_cost_is_symmetric_and_with_cost_does_not_mutate():
    g = topology.path(3)
    h = g.with_cost(2, 1, 2.5)
    assert h.cost(1, 2) == h.cost(2, 1) == 2.5
    assert g.cost(1, 2) == 1.0
    with pytest.raises(KeyError):
        g.with_cost(0, 2, 1.0)


def test_topology_shapes():
    assert len(topology.path(6).edges) == 5
    tree = topology.binary_tree(2)
    assert tree.n == 7 and len(tree.edges) == 6
    assert len(topology.ring(8).edges) == 8
    d = topology.diamond()
    assert d.n == 7 and len(d.edges) == 7
    assert d.neighbors(1) == [0, 2, 3] and d.name(6) == "T"
    lop = topology.lopsided_diamond()
    assert lop.n == 6 and len(lop.edges) == 6 and lop.neighbors(5) == [2, 4]


def test_cycle_parity():
    assert nx.is_bipartite(topology.diamond().to_networkx())
    assert not nx.is_bipartite(topology.lopsided_diamond().to_networkx())
    assert nx.is_tree(topology.binary_tree(2).to_networkx())


def test_random_costs_are_seeded_and_bounded():
    a = topology.random_costs(topology.path(6), 1)
    b = topology.random_costs(topology.path(6), 1)
    assert a.costs == b.costs
    assert all(0.5 <= c <= 2.0 for c in a.costs.values())
