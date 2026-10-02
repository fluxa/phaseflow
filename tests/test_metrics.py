import math

import numpy as np
import pytest

from phaseflow import metrics, topology


def test_shortest_reference_on_lopsided_diamond():
    g = topology.lopsided_diamond()  # upper A-B-T costs 2, lower A-C-E-T costs 3
    dist, opt = metrics.shortest(g, 5)
    assert list(dist) == [3, 2, 1, 2, 1, 0]
    assert opt[1] == {2} and opt[3] == {4} and opt[5] == set()
    cheap = g.with_cost(1, 3, 0.5).with_cost(3, 4, 0.5).with_cost(4, 5, 0.5)
    assert metrics.shortest(cheap, 5)[1][1] == {3}
    assert metrics.shortest(topology.diamond(), 6)[1][1] == {2, 3}


def test_follow_statuses():
    assert metrics.follow(np.array([0, 0, 1]), 2, 0) == ("delivered", [2, 1, 0])
    assert metrics.follow(np.array([0, 2, 1]), 1, 0)[0] == "loop"
    assert metrics.follow(np.array([0, -1, 1]), 2, 0)[0] == "dead"


def test_route_metrics_on_hand_built_weights():
    g = topology.lopsided_diamond()
    W = np.zeros((6, 6))
    W[0, 1] = 1
    W[1, 3], W[1, 2] = 2, 1  # A prefers the three-hop lower arm: wrong
    W[2, 5] = W[3, 4] = W[4, 5] = 1
    m = metrics.route_metrics(g, 5, W, np.zeros(6))
    assert m["nexthop_acc"] == pytest.approx(4 / 5)
    assert m["delivered_frac"] == 1.0
    assert m["max_stretch"] == pytest.approx(3 / 2)
    assert m["mean_stretch"] == pytest.approx((4 / 3 + 3 / 2 + 1 + 1 + 1) / 5)
    assert m["n_loop_nodes"] == 0 and m["n_dead_nodes"] == 0 and m["n_ties"] == 0
    assert m["min_margin"] == pytest.approx(1 / 3)


def test_loop_dead_and_dark_nodes_are_counted():
    g = topology.path(4)
    W = np.zeros((4, 4))
    W[1, 2] = W[2, 1] = 1  # 1 <-> 2 loop; node 3 has no positive weight
    canc = np.array([0.0, 0.2, 1.0, np.nan])
    m = metrics.route_metrics(g, 0, W, canc)
    assert m["n_loop_nodes"] == 2 and m["n_dead_nodes"] == 1
    assert m["delivered_frac"] == 0.0 and math.isnan(m["mean_stretch"])
    assert m["n_dark_nodes"] == 1
    assert m["max_cancellation"] == 1.0 and m["mean_cancellation"] == pytest.approx(0.6)


def test_tie_counts_as_correct_when_the_chosen_hop_is_optimal():
    g = topology.diamond()
    W = np.zeros((7, 7))
    W[0, 1] = 1
    W[1, 2] = W[1, 3] = 1  # exact tie at A; both arms are optimal
    W[2, 4] = W[4, 6] = W[3, 5] = W[5, 6] = 1
    m = metrics.route_metrics(g, 6, W, np.zeros(7))
    assert m["n_ties"] == 1 and m["nexthop_acc"] == 1.0 and m["min_margin"] == 0.0


def test_stochastic_metrics_simple_chain():
    g = topology.path(3)
    W = np.zeros((3, 3))
    W[1, 0] = W[1, 2] = 1  # node 1: to the destination or back to 2, each 1/2
    W[2, 1] = 1
    s = metrics.stochastic_metrics(g, 0, W)
    assert s["p_delivery_min"] == pytest.approx(1.0)
    # E1 = 1 + E2 / 2 and E2 = 1 + E1 give E1 = 3, E2 = 4; stretch = (3/1 + 4/2) / 2
    assert s["exp_stretch_mean"] == pytest.approx(2.5)
    assert s["entropy_mean"] == pytest.approx(math.log(2) / 2)


def test_stochastic_metrics_closed_loop_without_destination():
    g = topology.path(3)
    W = np.zeros((3, 3))
    W[1, 2] = W[2, 1] = 1
    s = metrics.stochastic_metrics(g, 0, W)
    assert s["p_delivery_min"] == 0.0
    assert s["exp_stretch_mean"] == math.inf


def test_circulation_and_split_share():
    P = np.zeros((3, 3))
    for a, b in ((0, 1), (1, 2), (2, 0)):  # unit current 0 -> 1 -> 2 -> 0
        P[b, a], P[a, b] = 1, -1
    assert metrics.circulation(P, (0, 1, 2)) == pytest.approx(1.0)
    assert metrics.circulation(np.abs(P), (0, 1, 2)) == pytest.approx(0.0)
    W = np.zeros((4, 4))
    W[1, 2], W[1, 3] = 3, 1
    assert metrics.split_share(W, 1, 2, 3) == pytest.approx(0.75)
    assert math.isnan(metrics.split_share(W, 0, 2, 3))
