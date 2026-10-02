import math

import numpy as np

from phaseflow import experiment0, scenarios
from phaseflow.models import Model
from phaseflow.perturbations import CostChange

LOP = {s.name: s for s in scenarios.static_scenarios()}["lop_unit"]


def test_settle_round():
    assert experiment0.settle_round(np.array([False, False, True, True])) == 3
    assert experiment0.settle_round(np.array([True, True])) == 1
    assert math.isnan(experiment0.settle_round(np.array([True, False])))
    assert math.isnan(experiment0.settle_round(np.array([], dtype=bool)))


def test_reconverges_after_a_cost_change_and_after_restoring_it():
    changes = (CostChange(1500, (2, 5), 3.0), CostChange(3000, (2, 5), 1.0))
    res = experiment0.run_dynamic(LOP, changes, Model("M0", kappa=2.0), 1.0, 4500)
    ph = res.phases
    assert list(ph["phase"]) == [0, 1, 2]
    assert not ph["diverged"].any() and ph["stable"].all()
    assert ph["settle_rounds"].notna().all() and ph["route_settle_rounds"].notna().all()
    assert list(ph["final_nexthop_acc"]) == [1.0, 1.0, 1.0]
    split = res.series.set_index("round")["split_upper"]
    assert split[1499] > 0.5 > split[2999]  # B-T raised to 3 makes the lower arm cheaper
    assert abs(split[4499] - split[1499]) < 1e-6  # no hysteresis


def test_logging_is_subsampled_but_keeps_phase_ends():
    changes = (CostChange(200, (2, 5), 3.0),)
    res = experiment0.run_dynamic(LOP, changes, Model("M2"), 3.5, 400)
    rounds = sorted(res.series["round"].unique())
    assert rounds[:60] == list(range(60)) and 199 in rounds and 399 in rounds and 61 not in rounds
    assert len(res.states) == len(rounds) * LOP.graph.n
    assert len(res.currents) == len(rounds) * 2 * len(LOP.graph.edges)


def test_divergence_is_flagged_and_stops_the_run():
    changes = (CostChange(3000, (2, 5), 3.0),)
    res = experiment0.run_dynamic(LOP, changes, Model("M3S", beta=1.0), 0.1, 6000)
    assert len(res.phases) == 1
    row = res.phases.iloc[0]
    assert row["diverged"] and not row["stable"]
    assert math.isnan(row["settle_rounds"]) and math.isnan(row["final_nexthop_acc"])
    assert np.isfinite(res.states["magnitude"]).all()


def test_dynamic_run_is_deterministic():
    changes = (CostChange(100, (2, 5), 3.0),)
    a = experiment0.run_dynamic(LOP, changes, Model("M3H", beta=1.0), 3.5, 200)
    b = experiment0.run_dynamic(LOP, changes, Model("M3H", beta=1.0), 3.5, 200)
    assert a.series.equals(b.series) and a.states.equals(b.states)
