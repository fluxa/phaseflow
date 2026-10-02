import numpy as np
import pytest

from phaseflow import oracle, topology
from phaseflow.models import Model, complex_matrix, orientation_signs, to_complex

G = topology.random_costs(topology.lopsided_diamond(), 3)
DEST = 5


def test_state_sizes():
    assert Model("M0", kappa=1.0).k == 1
    assert all(Model(kind).k == 2 for kind in ("M1", "M2", "M3H", "M3S"))
    assert not Model("M1").is_complex and Model("M3S").is_complex


def test_unknown_kind_is_rejected():
    with pytest.raises(ValueError):
        Model("M9")


def test_coupling_is_zero_off_edges():
    for m in (Model("M0", kappa=1.0), Model("M3H", beta=1.0), Model("M3S", beta=1.0)):
        K = m.coupling(G, DEST)
        for i in range(G.n):
            for j in range(G.n):
                if j not in G.neighbors(i):
                    assert not K[i, j].any()


def test_m3h_is_hermitian_and_m3s_is_symmetric_not_hermitian():
    H = 1j * complex_matrix(Model("M3H", beta=0.9, orientation="random:1").coupling(G, DEST))
    assert np.allclose(H, H.conj().T)
    S = complex_matrix(Model("M3S", beta=0.9, theta0=0.0).coupling(G, DEST))
    assert np.allclose(S, S.T)
    assert not np.allclose(S, S.conj().T)


def test_m3s_kickoff_form_is_minus_i_times_s():
    S = complex_matrix(Model("M3S", beta=0.9, theta0=0.0).coupling(G, DEST))
    K = complex_matrix(Model("M3S", beta=0.9).coupling(G, DEST))
    assert np.allclose(K, -1j * S)


def test_m2_equals_both_m3_models_at_beta_zero():
    K2 = Model("M2").coupling(G, DEST)
    assert np.allclose(K2, Model("M3H", beta=0.0).coupling(G, DEST))
    assert np.allclose(K2, Model("M3S", beta=0.0).coupling(G, DEST))


def test_orientation_conventions():
    assert set(orientation_signs(G, "id", DEST).values()) == {1}
    assert set(orientation_signs(G, "reverse", DEST).values()) == {-1}
    r = orientation_signs(G, "random:0", DEST)
    assert r == orientation_signs(G, "random:0", DEST)
    assert set(r.values()) <= {-1, 1}
    # Path from destination 0: for every key (i, j) with i < j, i is nearer the
    # destination, so "i receives from j" points toward it and gets -1.
    assert set(orientation_signs(topology.path(4), "oracle_hops", 0).values()) == {-1}
    with pytest.raises(ValueError):
        orientation_signs(G, "nope", DEST)


def test_variant_and_param_labels():
    assert Model("M0", kappa=2.0).variant == "M0" and Model("M0", kappa=2.0).param == 2.0
    assert Model("M2").param == 0.0
    assert Model("M3H", beta=0.5, orientation="random:3").variant == "M3H[random:3]"
    assert Model("M3S", beta=0.5).variant == "M3S[theta0=-pi/2]"
    assert Model("M3S", beta=0.5, theta0=0.0).variant == "M3S[theta0=0]"


def test_real_block_form_matches_direct_complex_solve():
    gamma = 0.4
    for m in (Model("M3H", beta=1.3, orientation="random:2"), Model("M3S", beta=1.3), Model("M2")):
        K, b = m.coupling(G, DEST), m.source(G.n, DEST)
        x = to_complex(oracle.fixed_point(K, b, gamma))
        direct = np.linalg.solve(gamma * np.eye(G.n) - complex_matrix(K), to_complex(b))
        assert np.allclose(x, direct)


def test_stability_regions():
    Kh = Model("M3H", beta=1.0).coupling(G, DEST)
    assert oracle.abscissa(Kh, 0.01) == pytest.approx(-0.01)
    Ks = Model("M3S", beta=1.0).coupling(G, DEST)
    assert oracle.abscissa(Ks, 0.1) > 0
    assert oracle.abscissa(Ks, 3.5) < 0  # 3.5 exceeds a * deg_max = 3


def test_loop_flux_is_odd_under_reversal_and_wrapped():
    from phaseflow.models import loop_flux

    g = topology.random_costs(topology.ring(7), 5)
    cycle = tuple(range(7))
    forward = loop_flux(g, orientation_signs(g, "id", 0), 0.3, cycle)
    backward = loop_flux(g, orientation_signs(g, "reverse", 0), 0.3, cycle)
    assert forward == pytest.approx(-backward) and abs(forward) > 0.1
    assert abs(loop_flux(g, orientation_signs(g, "id", 0), 50.0, cycle)) <= np.pi


def test_only_the_loop_flux_of_an_orientation_matters_on_a_single_cycle():
    from phaseflow import readout
    from phaseflow.models import loop_flux

    g, cycle, gamma = topology.ring(8), tuple(range(8)), 0.3  # unit costs: flux = beta * sum of signs

    def observe(convention):
        m = Model("M3H", beta=0.4, orientation=convention)
        K = m.coupling(g, 0)
        flux = loop_flux(g, orientation_signs(g, convention, 0), 0.4, cycle)
        return tuple(orientation_signs(g, convention, 0).values()), flux, readout.power_inflow(
            K, oracle.fixed_point(K, m.source(g.n, 0), gamma)
        )

    seen = [observe(f"random:{s}") for s in range(40)]
    same = [(a, b) for a in seen for b in seen if a[0] != b[0] and abs(a[1] - b[1]) < 1e-12]
    different = [(a, b) for a in seen for b in seen if abs(abs(a[1]) - abs(b[1])) > 1e-6]
    assert same and different
    assert all(np.allclose(a[2], b[2], atol=1e-12) for a, b in same)
    assert not np.allclose(different[0][0][2], different[0][1][2], atol=1e-6)
