import pytest

from phaseflow import metrics, scenarios, topology
from phaseflow.perturbations import CostChange, phases


def test_phases_split_the_schedule():
    g = topology.diamond()
    changes = [CostChange(10, (2, 4), 3.0), CostChange(20, (2, 4), 1.0)]
    ph = phases(g, changes, 30)
    assert [(start, end) for start, end, _ in ph] == [(0, 10), (10, 20), (20, 30)]
    assert [graph.cost(2, 4) for _, _, graph in ph] == [1.0, 3.0, 1.0]
    with pytest.raises(ValueError):
        phases(g, [CostChange(30, (2, 4), 3.0)], 30)
    with pytest.raises(ValueError):
        phases(g, [CostChange(10, (2, 4), 3.0), CostChange(10, (2, 4), 1.0)], 30)


def test_static_scenarios_have_the_intended_optimal_arms():
    scn = {s.name: s for s in scenarios.static_scenarios()}
    assert len(scn) == 15
    assert {s.family for s in scn.values()} == {"tree", "diamond", "lopsided", "ring"}

    def optimal_at_a(name):
        s = scn[name]
        return metrics.shortest(s.graph, s.dest)[1][1]

    assert optimal_at_a("diamond_equal") == {2, 3}
    assert optimal_at_a("diamond_lower+0.2") == {2}
    assert optimal_at_a("diamond_lower+1.0") == {2}
    assert optimal_at_a("diamond_upper+1.0") == {3}
    assert optimal_at_a("lop_unit") == {2}
    assert optimal_at_a("lop_lower_cheap") == {3}
    assert optimal_at_a("lop_near_tie") == {2}
    assert metrics.shortest(scn["ring8_heavy"].graph, 0)[1][3] == {4}  # detour the long way
    assert metrics.shortest(scn["ring7_heavy"].graph, 0)[1][3] == {4}


def test_model_grid_matches_the_spec():
    grid = scenarios.model_grid()
    assert len(grid) == 5 + 5 + 1 + 8 * 5 + 2 * 5
    assert len({(m.variant, m.param) for m in grid}) == len(grid)


def test_full_config_matches_the_spec():
    cfg = scenarios.full_config()
    assert cfg.gammas == (0.1, 0.3, 1.0, 2.0, 3.5)
    assert cfg.curve_gammas == (0.3, 3.5) and cfg.dynamic_gammas == (0.3, 1.0, 3.5)
    assert [c.name for c in cfg.curves] == ["diamond", "lopsided"]
    assert all(case.total_rounds == 7500 and len(case.changes) == 2 for case in cfg.dynamic_cases)
    assert [m.variant for m in cfg.dynamic_models] == [
        "M0", "M1", "M2", "M3H[id]", "M3S[theta0=-pi/2]", "M3S[theta0=0]",
    ]


def test_quick_config_is_a_small_subset():
    cfg = scenarios.quick_config()
    assert len(cfg.scenarios) == 4 and len(cfg.gammas) == 2
    assert all(case.total_rounds == 600 for case in cfg.dynamic_cases)
