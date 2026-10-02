import numpy as np
import pytest

from phaseflow import oracle, readout, topology
from phaseflow.models import Model

ALL = (
    Model("M0", kappa=1.0),
    Model("M1", kappa=1.0, kappa2=2.0),
    Model("M2"),
    Model("M3H", beta=0.7),
    Model("M3S", beta=0.7),
    Model("M3S", beta=0.7, theta0=0.0),
)
LOP = topology.random_costs(topology.lopsided_diamond(), 3)


def solve(m, g, dest, gamma):
    K, b = m.coupling(g, dest), m.source(g.n, dest)
    x = oracle.fixed_point(K, b, gamma)
    return K, b, x, readout.power_inflow(K, x)


def test_a1_two_node_sign_and_forwarding_direction():
    g = topology.path(2)
    for m in ALL:
        for gamma in (1.5, 3.0):  # above every model's threshold on two nodes
            _, _, x, P = solve(m, g, 0, gamma)
            assert P[1, 0] > 0
            assert P[1, 0] == pytest.approx(2 * gamma * (x[1] @ x[1]))
            nxt, _ = readout.next_hops(readout.weights(P), 0)
            assert nxt[1] == 0
    for m in (Model("M2"), Model("M3H", beta=0.7)):
        for gamma in (0.2, 1.0, 3.0):
            _, _, _, P = solve(m, g, 0, gamma)
            assert P[1, 0] > 0 and P[0, 1] == pytest.approx(-P[1, 0])


def test_a2_balance_identity():
    gamma = 3.5
    for m in ALL:
        _, b, x, P = solve(m, LOP, 5, gamma)
        lhs = P.sum(axis=1) + 2 * np.einsum("ia,ia->i", x, b)
        assert np.allclose(lhs, 2 * gamma * np.einsum("ia,ia->i", x, x), atol=1e-12)


def test_inflow_is_a_current_only_for_hermitian_models():
    for m in (Model("M2"), Model("M3H", beta=0.7, orientation="random:1")):
        _, _, _, P = solve(m, LOP, 5, 0.3)
        assert np.allclose(P, -P.T, atol=1e-12)
    _, _, _, P = solve(Model("M3S", beta=0.7), LOP, 5, 3.5)
    assert not np.allclose(P, -P.T, atol=1e-6)
    _, _, _, P = solve(Model("M0", kappa=1.0), LOP, 5, 3.5)
    assert (P >= 0).all()


def test_a3_tree_gauge_invariance_on_observables_not_raw_state():
    g = topology.random_costs(topology.binary_tree(2), 2)
    for dest in (0, 6):
        _, _, x0, P0 = solve(Model("M2"), g, dest, 0.3)
        for orientation in ("id", "reverse", "random:0", "random:1", "oracle_hops"):
            for beta in (0.5, 1.7):
                _, _, x, P = solve(Model("M3H", beta=beta, orientation=orientation), g, dest, 0.3)
                assert np.allclose(np.linalg.norm(x, axis=1), np.linalg.norm(x0, axis=1), atol=1e-12)
                assert np.allclose(P, P0, atol=1e-12)
    _, _, x, _ = solve(Model("M3H", beta=1.7, orientation="id"), g, 0, 0.3)
    _, _, x0, _ = solve(Model("M2"), g, 0, 0.3)
    assert not np.allclose(x, x0, atol=1e-3)  # raw states differ; only observables agree


def test_a4_bipartite_flux_sign_blindness_and_odd_cycle_counterexample():
    def flip_difference(g):
        dest = g.n - 1
        _, _, _, Pa = solve(Model("M3H", beta=0.8, orientation="id"), g, dest, 0.3)
        _, _, _, Pb = solve(Model("M3H", beta=0.8, orientation="reverse"), g, dest, 0.3)
        return np.abs(Pa - Pb).max()

    for g in (topology.ring(8), topology.diamond()):
        assert flip_difference(topology.random_costs(g, 4)) < 1e-12
    for g in (topology.ring(7), topology.lopsided_diamond()):
        assert flip_difference(topology.random_costs(g, 4)) > 1e-3


def test_a5_diamond_split_is_exactly_half_for_m3h():
    for seed in range(3):
        g = topology.random_costs(topology.diamond(), seed)
        for orientation in ("id", "random:0", "oracle_hops"):
            for beta in (0.5, 1.0, 2.0):
                for gamma in (0.1, 1.0):
                    _, _, _, P = solve(Model("M3H", beta=beta, orientation=orientation), g, 6, gamma)
                    W = readout.weights(P)
                    assert W[1, 2] > 0
                    assert W[1, 2] == pytest.approx(W[1, 3], rel=1e-9)


def test_cancellation_is_zero_for_real_models_and_positive_for_coherent_ones():
    for m in (Model("M0", kappa=1.0), Model("M1", kappa=1.0, kappa2=2.0)):
        K, _, x, _ = solve(m, LOP, 5, 3.5)
        assert np.allclose(readout.cancellation(K, x, m.is_complex), 0.0, atol=1e-12)
    m = Model("M3S", beta=1.0)
    K, _, x, _ = solve(m, LOP, 5, 3.5)
    assert readout.cancellation(K, x, m.is_complex)[1] > 1e-3  # node A has three neighbours


def test_next_hops_handles_ties_dead_nodes_and_destination():
    W = np.array(
        [
            [0.0, 0.0, 0.0, 0.0],
            [2.0, 0.0, 2.0, 0.0],  # tie between 0 and 2 -> lowest id
            [0.0, 0.0, 0.0, 0.0],  # dead
            [0.0, 1.0, 3.0, 0.0],
        ]
    )
    nxt, ties = readout.next_hops(W, 0)
    assert list(nxt) == [0, 0, -1, 2]
    assert list(ties) == [False, True, False, False]


def test_policy_rows_sum_to_one_or_zero():
    W = np.array([[0.0, 5.0, 0.0, 0.0], [2.0, 0.0, 2.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 1.0, 3.0, 0.0]])
    Pm = readout.policy(W, 0)
    assert np.allclose(Pm.sum(axis=1), [0.0, 1.0, 0.0, 1.0])
    assert Pm[3, 2] == pytest.approx(0.75)


def test_zero_field_does_not_crash_the_readout():
    m = Model("M3S", beta=1.0)
    K = m.coupling(topology.path(3), 0)
    x = np.zeros((3, 2))
    assert np.isnan(readout.cancellation(K, x, True)).all()
    nxt, ties = readout.next_hops(readout.weights(readout.power_inflow(K, x)), 0)
    assert list(nxt) == [0, -1, -1] and not ties.any()


def test_m1_with_equal_scales_forwards_like_m0():
    _, _, _, P0 = solve(Model("M0", kappa=1.0), LOP, 5, 3.5)
    _, _, _, P1 = solve(Model("M1", kappa=1.0, kappa2=1.0), LOP, 5, 3.5)
    assert np.allclose(P1, 2 * P0)
