import json

from phaseflow import experiment0

FILES = (
    "config.json", "summary.json", "route_metrics.csv", "static_readouts.csv", "transfer.csv",
    "discrimination.csv", "dynamic_metrics.csv", "dynamic_series.csv", "states.csv", "currents.csv",
)


def test_quick_run_writes_every_output_and_is_deterministic(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    summary = experiment0.main(a, quick=True)
    experiment0.main(b, quick=True, make_plots=False)
    for name in FILES:
        assert (a / name).stat().st_size > 0
        assert (a / name).read_bytes() == (b / name).read_bytes()
    assert len(list((a / "plots").glob("*.png"))) >= 5
    loaded = json.loads((a / "summary.json").read_text())
    assert loaded["control_checks"]["tree_gauge_max_rel_diff"] < 1e-10
    assert set(loaded["transfer_by_group"]) >= {"M0", "M1", "M2", "M3H", "M3S"}
    assert summary["n_static_rows"] == loaded["n_static_rows"]
    config = json.loads((a / "config.json").read_text())
    assert config["gammas"] == [0.3, 3.5] and "numpy" in config["versions"]
