import numpy as np

from phaseflow import engine, oracle, topology
from phaseflow.models import Model

G = topology.random_costs(topology.lopsided_diamond(), 3)
DEST = 5


def test_converges_to_oracle_fixed_point():
    cases = (
        (Model("M3H", beta=1.0), 0.5),
        (Model("M3S", beta=1.0), 3.5),
        (Model("M0", kappa=2.0), 1.0),
        (Model("M1", kappa=2.0, kappa2=4.0), 1.0),
    )
    for m, gamma in cases:
        K, b = m.coupling(G, DEST), m.source(G.n, DEST)
        x, rounds, diverged = engine.run(np.zeros((G.n, m.k)), K, b, gamma, 4000)
        assert not diverged and rounds == 4000
        assert np.allclose(x, oracle.fixed_point(K, b, gamma), atol=1e-9)


def test_step_size_uses_only_local_coupling_strength():
    K = Model("M2").coupling(G, DEST)
    degree = np.array([len(G.neighbors(i)) for i in range(G.n)])
    assert np.allclose(engine.step_sizes(K, 0.5), 0.5 / (0.25 + degree**2))


def test_update_is_local():
    m = Model("M3S", beta=1.0)
    K, b = m.coupling(G, DEST), m.source(G.n, DEST)
    alpha = engine.step_sizes(K, 3.5)
    x = np.random.default_rng(0).normal(size=(G.n, 2))
    y = x.copy()
    y[4] += 10.0  # node 4 (E) has neighbours 3 (C) and 5 (T) only
    a, c = engine.step(x, K, b, 3.5, alpha), engine.step(y, K, b, 3.5, alpha)
    changed = {i for i in range(G.n) if not np.allclose(a[i], c[i])}
    assert changed == {3, 4, 5}


def test_hermitian_iteration_converges_at_small_damping_on_irregular_graph():
    for beta in (0.0, 1.0):
        K = Model("M3H", beta=beta, orientation="random:1").coupling(G, DEST)
        for gamma in (0.05, 0.3, 2.0):
            assert oracle.iteration_radius(K, gamma, engine.step_sizes(K, gamma)) < 1.0


def test_divergence_below_threshold_is_predicted_and_detected():
    m = Model("M3S", beta=1.0)
    K, b = m.coupling(G, DEST), m.source(G.n, DEST)
    assert oracle.abscissa(K, 0.1) > 0
    assert oracle.iteration_radius(K, 0.1, engine.step_sizes(K, 0.1)) > 1.0
    x, rounds, diverged = engine.run(np.zeros((G.n, 2)), K, b, 0.1, 20000)
    assert diverged and rounds < 20000
    assert np.isfinite(x).all()
