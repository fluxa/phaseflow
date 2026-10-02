import math

import numpy as np
import pandas as pd
import pytest

from phaseflow import experiment0, scenarios
from phaseflow.models import Model

BY_NAME = {s.name: s for s in scenarios.static_scenarios()}


def test_evaluate_stable_cell_has_metrics():
    row, P = experiment0.evaluate(BY_NAME["lop_unit"], Model("M0", kappa=2.0), 1.0)
    assert row["stable"] and row["converges"]
    assert row["nexthop_acc"] == 1.0 and row["split_upper"] > 0.5
    assert P.shape == (6, 6)


def test_unstable_cell_reports_no_readout():
    row, P = experiment0.evaluate(BY_NAME["lop_unit"], Model("M3S", beta=1.0), 0.1)
    assert not row["stable"] and P is None
    assert math.isnan(row["nexthop_acc"]) and math.isnan(row["split_upper"])


def test_tie_at_equal_diamond_is_recorded_and_correct():
    row, _ = experiment0.evaluate(BY_NAME["diamond_equal"], Model("M0", kappa=2.0), 1.0)
    assert row["n_ties"] == 1 and row["nexthop_acc"] == 1.0
    assert row["split_upper"] == pytest.approx(0.5)


def test_static_sweep_is_complete_and_deterministic():
    scns = [BY_NAME["path6_rand"], BY_NAME["diamond_lower+1.0"]]
    models = [Model("M2"), Model("M3H", beta=1.0), Model("M3S", beta=1.0), Model("M0", kappa=1.0)]
    a = experiment0.run_static(scns, models, (0.3, 3.5))
    b = experiment0.run_static(scns, models, (0.3, 3.5))
    assert len(a.metrics) == 2 * 4 * 2
    pd.testing.assert_frame_equal(a.metrics, b.metrics)
    pd.testing.assert_frame_equal(a.readouts, b.readouts)
    assert set(a.readouts.columns) >= {"scenario", "variant", "param", "gamma", "i", "j", "inflow", "weight"}


def test_control_checks_on_a_small_sweep():
    scns = [BY_NAME["path6_rand"], BY_NAME["btree7_leaf"], BY_NAME["diamond_lower+1.0"]]
    models = [Model("M2")] + [
        Model("M3H", beta=beta, orientation=o) for beta in (0.5, 2.0) for o in ("id", "random:0")
    ]
    checks = experiment0.control_checks(experiment0.run_static(scns, models, (0.3, 1.0)))
    assert checks["tree_gauge_max_rel_diff"] < 1e-10
    assert checks["diamond_m3h_split_max_dev"] < 1e-9
    assert checks["tree_cells_all_correct"] == {"M2": 1.0, "M3H": 1.0}


def test_transfer_table_marks_unstable_cells():
    scns = [BY_NAME["diamond_lower+1.0"], BY_NAME["lop_unit"], BY_NAME["path6"]]
    res = experiment0.run_static(scns, [Model("M3S", beta=1.0), Model("M0", kappa=2.0)], (0.1, 3.5))
    t = experiment0.transfer_table(res.metrics).set_index(["variant", "gamma"])
    assert (t["n_scenarios"] == 2).all()  # the tree is excluded
    unstable = t.loc[("M3S[theta0=-pi/2]", 0.1)]
    assert not unstable["stable_all"] and np.isnan(unstable["mean_acc"])
    stable = t.loc[("M0", 3.5)]
    assert stable["stable_all"] and stable["min_acc"] == 1.0


def test_summarize_groups_keeps_the_oracle_diagnostic_apart():
    scns = [BY_NAME["lop_unit"], BY_NAME["lop_lower_cheap"]]
    models = [Model("M3H", beta=1.0, orientation=o) for o in ("id", "reverse", "oracle_hops")]
    res = experiment0.run_static(scns, models, (0.3,))
    groups = experiment0.summarize_groups(experiment0.transfer_table(res.metrics))
    assert groups["M3H"]["n_cells"] == 2
    assert groups["M3H-oracle (diagnostic)"]["n_cells"] == 1
    assert "M3H worst arbitrary orientation" in groups


def test_discrimination_curves_cover_both_graphs():
    disc = experiment0.run_discrimination(scenarios.quick_config())
    assert set(disc["curve"]) == {"diamond", "lopsided"}
    m0 = disc[(disc["curve"] == "diamond") & (disc["variant"] == "M0")].sort_values("swept_cost")
    assert m0["split_upper"].is_monotonic_increasing  # dearer lower arm -> more weight on the upper arm


def test_m3h_rows_carry_the_loop_flux_and_other_models_do_not():
    lop = BY_NAME["lop_unit"]
    forward, _ = experiment0.evaluate(lop, Model("M3H", beta=1.0, orientation="id"), 0.3)
    backward, _ = experiment0.evaluate(lop, Model("M3H", beta=1.0, orientation="reverse"), 0.3)
    assert forward["loop_flux"] == pytest.approx(-backward["loop_flux"]) and forward["loop_flux"] != 0
    assert math.isnan(experiment0.evaluate(lop, Model("M0", kappa=2.0), 1.0)[0]["loop_flux"])
    assert math.isnan(experiment0.evaluate(BY_NAME["path6"], Model("M3H", beta=1.0), 0.3)[0]["loop_flux"])


def test_tie_between_an_optimal_and_a_wrong_hop_scores_half_under_fair_accuracy():
    row, _ = experiment0.evaluate(BY_NAME["diamond_lower+1.0"], Model("M2"), 1.0)
    assert row["n_ties"] == 1 and row["nexthop_acc"] == 1.0  # right only because B has the lower id
    assert row["nexthop_acc_fair"] == pytest.approx(5.5 / 6)
