import json
import math

import pandas as pd
import pytest

from phaseflow import analysis, experiment0, scenarios
from phaseflow.models import Model

BY_NAME = {s.name: s for s in scenarios.static_scenarios()}
COLUMNS = ["scenario", "family", "kind", "variant", "param", "gamma", "stable", "nexthop_acc_fair"]
HAND = pd.DataFrame(
    [
        ("easy", "ring", "M2", "M2", 0.0, 1.0, True, 1.0),
        ("hard", "ring", "M2", "M2", 0.0, 1.0, True, 0.5),
        ("leaf", "tree", "M2", "M2", 0.0, 1.0, True, 0.2),
        ("easy", "ring", "M0", "M0", 1.0, 1.0, True, 1.0),
        ("hard", "ring", "M0", "M0", 1.0, 1.0, True, 1.0),
        ("hard", "ring", "M0", "M0", 2.0, 1.0, False, math.nan),
        ("easy", "ring", "M3H", "M3H[oracle_hops]", 1.0, 1.0, True, 0.5),
        ("hard", "ring", "M3H", "M3H[id]", 1.0, 1.0, True, 0.5),
    ],
    columns=COLUMNS,
)


def counts(table, group, columns):
    return table.set_index("group").loc[group, columns].tolist()


def test_cost_sensitive_scenarios_are_the_cyclic_ones_the_cost_blind_model_gets_wrong():
    assert analysis.cost_sensitive(HAND) == ["hard"]


def test_versus_m2_counts_stable_cells_only_and_keeps_the_oracle_apart():
    t = analysis.versus_m2(HAND)
    columns = ["cells", "better", "equal", "worse"]
    assert counts(t, "M0", columns) == [2, 1, 1, 0]
    assert counts(t, "M3H", columns) == [1, 0, 1, 0]
    assert counts(t, "M3H-oracle", columns) == [1, 0, 0, 1]
    assert "M2" not in set(t["group"])


def test_perfect_tally_on_named_scenarios():
    t = analysis.perfect_tally(HAND, ["hard"])
    assert counts(t, "M0", ["stable", "perfect"]) == [1, 1]
    assert counts(t, "M3H", ["stable", "perfect"]) == [1, 0]


def test_cell_summary_requires_stability_on_every_cyclic_scenario():
    t = analysis.cell_summary(HAND)
    assert counts(t, "M0", ["cells", "stable_all", "perfect"]) == [2, 1, 1]  # kappa=2 is unstable on "hard"


def test_flux_sign_table_groups_m3h_cells_by_chirality():
    models = [Model("M3H", beta=0.5, orientation=o) for o in ("id", "reverse", "random:4")]
    res = experiment0.run_static([BY_NAME["ring7_heavy"]], models, (0.3,))
    t = analysis.flux_sign_table(res.metrics, "ring7_heavy")
    assert sorted(t["sign"]) == [-1, 1]
    assert t["cells"].sum() == 3 and t["distinct_flux"].sum() == 2  # random:4 repeats id on seven edges


def test_duplicate_orientations_are_reported():
    names = ("ring7_unit", "ring8_unit", "lop_unit")
    dup = analysis.duplicate_orientations([BY_NAME[n] for n in names], scenarios.ORIENTATIONS)
    assert ["id", "random:4"] in dup["ring7_unit"]
    assert dup["ring8_unit"] == []
    assert any({"id", "oracle_hops"} <= set(group) for group in dup["lop_unit"])


def test_mirroring_the_ring_swaps_which_chirality_is_right():
    t = analysis.mirrored_ring_check(betas=(0.5,), orientations=("id", "reverse"), gammas=(0.3,))
    acc = t.set_index(["scenario", "variant"])["nexthop_acc_fair"]
    assert set(t["scenario"]) == {"ring7_heavy", "ring7_heavy_mirrored"}
    assert acc["ring7_heavy", "M3H[id]"] == pytest.approx(acc["ring7_heavy_mirrored", "M3H[reverse]"])
    assert acc["ring7_heavy", "M3H[reverse]"] == pytest.approx(acc["ring7_heavy_mirrored", "M3H[id]"])
    assert acc["ring7_heavy", "M3H[id]"] != pytest.approx(acc["ring7_heavy", "M3H[reverse]"])


def test_lopsided_flux_scan_reports_the_minimum_and_where():
    scan = analysis.lopsided_flux_scan(n_flux=13, costs=(0.3, 1.0), gammas=(0.1, 1.0))
    assert 0.0 <= scan["min_upper_share"] <= 1.0
    assert set(scan["at"]) == {"flux", "gamma", "lower_cost"}


def test_main_writes_posthoc_json(tmp_path):
    experiment0.main(tmp_path, quick=True, make_plots=False)
    analysis.main(tmp_path)
    loaded = json.loads((tmp_path / "posthoc.json").read_text())
    assert set(loaded) >= {
        "note", "cost_sensitive_scenarios", "cell_summary", "versus_m2", "perfect_on_cost_sensitive",
        "flux_sign", "mirrored_ring", "matched_gamma", "stochastic_policy", "duplicate_orientations",
        "distinct_m3h_cells", "lopsided_flux_scan", "dynamic_outcomes",
    }
    assert "post-hoc" in loaded["note"].lower()
