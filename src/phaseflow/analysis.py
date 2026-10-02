"""Post-hoc analyses of the Experiment 0 sweep.

Everything in this module was defined after looking at the results. It lives
in the repository so that every number quoted in docs/experiment-0.md can be
regenerated; none of it was planned before the run.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from . import metrics, oracle, readout
from .experiment0 import _clean, run_static
from .models import Model, orientation_signs, rot
from .scenarios import ORIENTATIONS, SCALES, Scenario, _lopsided, static_scenarios
from .topology import ring

TOL = 1e-9
FAIR = "nexthop_acc_fair"  # a tie scores the share of tied hops that are optimal
NOTE = (
    "Post-hoc: every table here was defined after the sweep was inspected. "
    "Accuracy is tie-fair unless a key says otherwise."
)


def group_of(m: pd.DataFrame) -> np.ndarray:
    """Reporting group: the oracle orientation and the two M3-S rules are kept apart."""
    return np.where(
        m["variant"] == "M3H[oracle_hops]", "M3H-oracle", np.where(m["kind"] == "M3S", m["variant"], m["kind"])
    )


def cost_sensitive(m: pd.DataFrame, column: str = FAIR) -> list[str]:
    """Cyclic scenarios where the cost-blind model M2 gets at least one node wrong."""
    m2 = m[(m["kind"] == "M2") & (m["family"] != "tree")]
    return sorted(m2.loc[m2[column] < 1 - TOL, "scenario"].unique())


def versus_m2(m: pd.DataFrame, column: str = FAIR) -> pd.DataFrame:
    """Each stable cyclic scenario-cell against M2 at the same scenario and gamma."""
    cyc = m[(m["family"] != "tree") & m["stable"]]
    ref = cyc[cyc["kind"] == "M2"].set_index(["scenario", "gamma"])[column]
    rest = cyc[cyc["kind"] != "M2"]
    diff = rest[column].to_numpy() - ref.reindex(pd.MultiIndex.from_frame(rest[["scenario", "gamma"]])).to_numpy()
    table = pd.DataFrame({"group": group_of(rest), "diff": diff})
    return table.groupby("group")["diff"].agg(
        cells="size",
        better=lambda d: int((d > TOL).sum()),
        equal=lambda d: int((d.abs() <= TOL).sum()),
        worse=lambda d: int((d < -TOL).sum()),
    ).reset_index()


def perfect_tally(m: pd.DataFrame, names: list[str], column: str = FAIR) -> pd.DataFrame:
    """Stable cells on the named scenarios, and how many forward every node correctly."""
    sub = m[m["scenario"].isin(names) & m["stable"]]
    table = pd.DataFrame({"group": group_of(sub), "perfect": sub[column].to_numpy() >= 1 - TOL})
    return table.groupby("group")["perfect"].agg(stable="size", perfect="sum").astype(int).reset_index()


def cell_summary(m: pd.DataFrame, column: str = FAIR) -> pd.DataFrame:
    """Per parameter cell over all cyclic scenarios: usable, perfect, and the accuracy range."""
    cyc = m[m["family"] != "tree"].assign(group=group_of(m[m["family"] != "tree"]))
    rows = []
    for (group, _, _, _), grp in cyc.groupby(["group", "variant", "param", "gamma"]):
        ok = bool(grp["stable"].all())
        rows.append({
            "group": group, "stable_all": ok,
            "perfect": ok and grp[column].min() >= 1 - TOL,
            "mean_acc": grp[column].mean() if ok else math.nan,
        })
    cells = pd.DataFrame(rows)
    return cells.groupby("group").agg(
        cells=("stable_all", "size"), stable_all=("stable_all", "sum"), perfect=("perfect", "sum"),
        mean_acc_min=("mean_acc", "min"), mean_acc_max=("mean_acc", "max"),
    ).reset_index()


def flux_sign_table(m: pd.DataFrame, scenario: str, column: str = FAIR) -> pd.DataFrame:
    """M3-H accuracy on one single-cycle scenario, by gamma and the sign of the wrapped loop flux."""
    sub = m[(m["scenario"] == scenario) & (m["kind"] == "M3H")]
    flux = sub["loop_flux"].round(9)
    table = pd.DataFrame({
        "gamma": sub["gamma"].to_numpy(), "sign": np.sign(flux).astype(int).to_numpy(),
        "flux": flux.to_numpy(), "acc": sub[column].to_numpy(),
    })
    return table.groupby(["gamma", "sign"]).agg(
        cells=("acc", "size"), distinct_flux=("flux", "nunique"), min_acc=("acc", "min"), max_acc=("acc", "max"),
    ).reset_index()


def matched_gamma(m: pd.DataFrame, column: str = FAIR, gamma_min: float = 2.0) -> pd.DataFrame:
    """Mean accuracy over cyclic scenario-cells at strong damping, where few cells are unstable."""
    sub = m[(m["family"] != "tree") & (m["gamma"] >= gamma_min)]
    table = pd.DataFrame({"group": group_of(sub), "stable": sub["stable"].to_numpy(), "acc": sub[column].to_numpy()})
    return table.groupby("group").agg(
        cells=("stable", "size"), unstable=("stable", lambda s: int((~s).sum())), mean_acc=("acc", "mean"),
    ).reset_index()


def stochastic_table(m: pd.DataFrame) -> pd.DataFrame:
    """Proportional forwarding p(i -> j) ∝ max(P, 0): expected-cost stretch against M2."""
    cyc = m[(m["family"] != "tree") & m["stable"]]
    ref = cyc[cyc["kind"] == "M2"].set_index(["scenario", "gamma"])["exp_stretch_mean"]
    diff = cyc["exp_stretch_mean"].to_numpy() - ref.reindex(pd.MultiIndex.from_frame(cyc[["scenario", "gamma"]])).to_numpy()
    group = group_of(cyc)
    damping = np.where(cyc["gamma"] >= 1, " (gamma >= 1)", " (gamma < 1)")
    group = np.where(cyc["kind"] == "M3H", np.char.add(group.astype(str), damping), group)
    table = pd.DataFrame({
        "group": group, "stretch": cyc["exp_stretch_mean"].to_numpy(), "diff": diff,
        "delivery": cyc["p_delivery_min"].to_numpy(), "entropy": cyc["entropy_mean"].to_numpy(),
    })
    return table.groupby("group").agg(
        cells=("stretch", "size"), exp_stretch=("stretch", "mean"), min_delivery=("delivery", "min"),
        entropy=("entropy", "mean"),
        lower_than_m2=("diff", lambda d: int((d < -TOL).sum())), higher_than_m2=("diff", lambda d: int((d > TOL).sum())),
    ).reset_index()


def duplicate_orientations(scns, conventions) -> dict[str, list[list[str]]]:
    """Per scenario, the sets of orientation conventions that are edge-for-edge identical."""
    out = {}
    for scn in scns:
        seen: dict[tuple, list[str]] = {}
        for convention in conventions:
            signs = tuple(orientation_signs(scn.graph, convention, scn.dest).values())
            seen.setdefault(signs, []).append(convention)
        out[scn.name] = [names for names in seen.values() if len(names) > 1]
    return out


def distinct_m3h_cells(m: pd.DataFrame) -> dict:
    """On a single-cycle graph only (flux, gamma) matters, so count distinct physical cells."""
    sub = m[(m["kind"] == "M3H") & (m["family"] != "tree") & (m["variant"] != "M3H[oracle_hops]")]
    keys = pd.DataFrame({"s": sub["scenario"], "f": sub["loop_flux"].round(9), "g": sub["gamma"]})
    return {"scenario_cells": len(sub), "distinct": len(keys.drop_duplicates())}


def mirrored_ring_check(betas=SCALES, orientations=ORIENTATIONS, gammas=(0.1, 0.3)) -> pd.DataFrame:
    """M3-H on ring7 with the heavy edge at (1, 2) and at its mirror image (5, 6).

    The orientation, and so the flux, is the same in both; only which side is
    dear differs. If the flux carried cost information the same cells would be
    right in both scenarios.
    """
    cycle = tuple(range(7))
    scns = [
        Scenario("ring7_heavy", "ring", ring(7).with_cost(1, 2, 4.0), 0, cycle),
        Scenario("ring7_heavy_mirrored", "ring", ring(7).with_cost(5, 6, 4.0), 0, cycle),
    ]
    models = [Model("M3H", beta=b, orientation=o) for o in orientations for b in betas]
    m = run_static(scns, models, gammas).metrics
    return m[["scenario", "variant", "param", "gamma", "loop_flux", FAIR]].assign(
        sign=np.sign(m["loop_flux"].round(9)).astype(int)
    )


def lopsided_flux_scan(n_flux: int = 73, costs=tuple(np.linspace(0.3, 1.5, 13)), gammas=(0.02, 0.05, 0.1, 0.3, 1.0, 3.5)) -> dict:
    """Smallest share of A's weight on the two-hop arm over every loop flux, cost and gamma.

    By gauge invariance the whole flux can sit on one edge, so this covers
    every orientation and every beta at once.
    """
    base = Model("M2")
    best = {"min_upper_share": math.inf, "at": {}}
    for cost in costs:
        scn = _lopsided("scan", float(cost))
        b = base.source(scn.graph.n, scn.dest)
        for flux in np.linspace(-math.pi, math.pi, n_flux):
            K = base.coupling(scn.graph, scn.dest)
            K[4, 5], K[5, 4] = rot(-math.pi / 2 + flux), rot(-math.pi / 2 - flux)
            for gamma in gammas:
                W = readout.weights(readout.power_inflow(K, oracle.fixed_point(K, b, gamma)))
                share = metrics.split_share(W, *scn.split)
                if share < best["min_upper_share"]:
                    best = {"min_upper_share": share, "at": {"flux": float(flux), "gamma": gamma, "lower_cost": float(cost)}}
    return best


def main(out_dir="results/experiment0-smoke") -> dict:
    out = Path(out_dir)
    m = pd.read_csv(out / "route_metrics.csv")
    dyn = pd.read_csv(out / "dynamic_metrics.csv")
    hard = cost_sensitive(m)
    cyclic = [s for s in static_scenarios() if s.family != "tree" and s.name in set(m["scenario"])]
    odd = [s.name for s in cyclic if s.cycle and len(s.cycle) % 2 == 1]
    mirror = mirrored_ring_check()
    dynamic_columns = [
        "scenario", "variant", "gamma", "phase", "stable", "diverged",
        "settle_rounds", "route_settle_rounds", "loop_rounds", "final_nexthop_acc",
    ]

    def records(table: pd.DataFrame) -> list[dict]:
        return table.to_dict("records")

    result = _clean({
        "note": NOTE,
        "cost_sensitive_scenarios": hard,
        "cost_sensitive_scenarios_lowest_id_ties": cost_sensitive(m, "nexthop_acc"),
        "cell_summary": records(cell_summary(m)),
        "versus_m2": records(versus_m2(m)),
        "versus_m2_lowest_id_ties": records(versus_m2(m, "nexthop_acc")),
        "perfect_on_cost_sensitive": records(perfect_tally(m, hard)),
        "perfect_per_cost_sensitive_scenario": {name: records(perfect_tally(m, [name])) for name in hard},
        "flux_sign": {name: records(flux_sign_table(m, name)) for name in odd},
        "mirrored_ring": {
            name: records(flux_sign_table(mirror.assign(kind="M3H"), name))
            for name in ("ring7_heavy", "ring7_heavy_mirrored")
        },
        "matched_gamma": records(matched_gamma(m)),
        "stochastic_policy": records(stochastic_table(m)),
        "duplicate_orientations": duplicate_orientations(cyclic, ORIENTATIONS),
        "distinct_m3h_cells": distinct_m3h_cells(m),
        "lopsided_flux_scan": lopsided_flux_scan(),
        "dynamic_outcomes": records(dyn[dynamic_columns]),
    })
    (out / "posthoc.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Post-hoc analyses of the Experiment 0 sweep.")
    parser.add_argument("--out", default="results/experiment0-smoke")
    main(parser.parse_args().out)
