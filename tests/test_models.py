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
