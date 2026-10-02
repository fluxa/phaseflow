"""Experiment 0 smoke test: static sweep, discrimination curves, dynamic run."""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import platform
import subprocess
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd

from . import engine, metrics, oracle, readout
from .models import Model, loop_flux, orientation_signs
from .perturbations import phases
from .scenarios import Config, Scenario, full_config, quick_config

STABLE_TOL = 1e-9
SETTLE_TOL = 1e-6
DIVERGE_AT = 1e12
LOG_DENSE = 60
LOG_EVERY = 50
ID_COLS = ["scenario", "family", "kind", "variant", "param", "gamma"]
METRIC_COLS = [
    "nexthop_acc", "nexthop_acc_fair", "delivered_frac", "n_loop_nodes", "n_dead_nodes", "n_ties", "n_dark_nodes",
    "mean_stretch", "max_stretch", "mean_cancellation", "max_cancellation", "min_margin",
    "p_delivery_min", "exp_stretch_mean", "entropy_mean", "circulation", "split_upper",
]
ORACLE_GROUP = "M3H-oracle (diagnostic)"


def _ids(scn: Scenario, model: Model, gamma: float) -> dict:
    return {
        "scenario": scn.name, "family": scn.family, "kind": model.kind,
        "variant": model.variant, "param": model.param, "gamma": gamma,
    }


def evaluate(scn: Scenario, model: Model, gamma: float) -> tuple[dict, np.ndarray | None]:
    """One cell at the oracle fixed point. Unstable cells get NaN metrics and no inflow."""
    g, dest = scn.graph, scn.dest
    K, b = model.coupling(g, dest), model.source(g.n, dest)
    absc = oracle.abscissa(K, gamma)
    rho = oracle.iteration_radius(K, gamma, engine.step_sizes(K, gamma))
    row = {
        **_ids(scn, model, gamma),
        "abscissa": absc, "rho_discrete": rho,
        "stable": absc < -STABLE_TOL, "converges": rho < 1 - STABLE_TOL,
        "loop_flux": math.nan,
        **dict.fromkeys(METRIC_COLS, math.nan),
    }
    if model.kind == "M3H" and scn.cycle:
        signs = orientation_signs(g, model.orientation, dest)
        row["loop_flux"] = loop_flux(g, signs, model.beta, scn.cycle)
    if not row["stable"]:
        return row, None
    x = oracle.fixed_point(K, b, gamma)
    P = readout.power_inflow(K, x)
    W = readout.weights(P)
    row.update(metrics.route_metrics(g, dest, W, readout.cancellation(K, x, model.is_complex)))
    row.update(metrics.stochastic_metrics(g, dest, W))
    if scn.cycle:
        row["circulation"] = metrics.circulation(P, scn.cycle)
    if scn.split:
        row["split_upper"] = metrics.split_share(W, *scn.split)
    return row, P


@dataclass
class StaticResult:
    metrics: pd.DataFrame
    readouts: pd.DataFrame
    inflows: dict  # (scenario, variant, param, gamma) -> P


def run_static(scns, models, gammas) -> StaticResult:
    rows, edges, inflows = [], [], {}
    for scn in scns:
        for model in models:
            for gamma in gammas:
                row, P = evaluate(scn, model, gamma)
                rows.append(row)
                if P is None:
                    continue
                inflows[(scn.name, model.variant, model.param, gamma)] = P
                ids = {key: row[key] for key in ID_COLS}
                for i in range(scn.graph.n):
                    for j in scn.graph.neighbors(i):
                        edges.append({**ids, "i": i, "j": j, "inflow": P[i, j], "weight": max(P[i, j], 0.0)})
    return StaticResult(pd.DataFrame(rows), pd.DataFrame(edges), inflows)


def control_checks(res: StaticResult) -> dict:
    """Residuals of the controls derived in the spec (A3 tree gauge, A5 diamond split)."""
    m = res.metrics
    family = dict(zip(m["scenario"], m["family"]))
    worst = 0.0
    for (name, variant, _, gamma), P in res.inflows.items():
        ref = res.inflows.get((name, "M2", 0.0, gamma))
        if variant.startswith("M3H") and family[name] == "tree" and ref is not None:
            worst = max(worst, float(np.abs(P - ref).max() / np.abs(ref).max()))
    dia = m[(m["family"] == "diamond") & (m["kind"] == "M3H") & m["stable"]]
    tree = m[(m["family"] == "tree") & m["stable"]]
    correct = (tree["nexthop_acc"] == 1.0) & (tree["delivered_frac"] == 1.0)
    return {
        "tree_gauge_max_rel_diff": worst,
        "diamond_m3h_split_max_dev": float((dia["split_upper"] - 0.5).abs().max()) if len(dia) else math.nan,
        "tree_cells_all_correct": {kind: float(grp.mean()) for kind, grp in correct.groupby(tree["kind"])},
    }


def transfer_table(m: pd.DataFrame) -> pd.DataFrame:
    """Per parameter cell, aggregate over the cyclic scenarios. NaN if unstable on any."""
    rows = []
    cyclic = m[m["family"] != "tree"]
    for (kind, variant, param, gamma), grp in cyclic.groupby(["kind", "variant", "param", "gamma"]):
        ok = bool(grp["stable"].all())

        def agg(column: str, how: str) -> float:
            return float(getattr(grp[column], how)()) if ok else math.nan

        rows.append({
            "kind": kind, "variant": variant, "param": param, "gamma": gamma,
            "n_scenarios": len(grp), "stable_all": ok,
            "min_acc": agg("nexthop_acc", "min"), "mean_acc": agg("nexthop_acc", "mean"),
            "min_delivered": agg("delivered_frac", "min"), "max_stretch": agg("max_stretch", "max"),
            "mean_cancellation": agg("mean_cancellation", "mean"),
        })
    return pd.DataFrame(rows)


def summarize_groups(t: pd.DataFrame) -> dict:
    """Describe the whole distribution of cells per model. The best cell is listed, never alone."""
    group = np.where(t["variant"] == "M3H[oracle_hops]", ORACLE_GROUP, t["kind"])
    out = {}
    for name, grp in t.groupby(group):
        ok = grp[grp["stable_all"]]
        perfect = ok[(ok["min_acc"] == 1.0) & (ok["min_delivered"] == 1.0)]
        best = ok.sort_values(["mean_acc", "min_acc"], ascending=False).head(1)
        columns = ["variant", "param", "gamma", "min_acc", "mean_acc", "min_delivered"]
        out[name] = {
            "n_cells": len(grp), "n_stable_all": len(ok), "n_perfect": len(perfect),
            "mean_acc_quartiles": [float(q) for q in ok["mean_acc"].quantile([0, 0.25, 0.5, 0.75, 1])] if len(ok) else [],
            "best_cell": best[columns].iloc[0].to_dict() if len(best) else None,
        }
    arbitrary = t[(t["kind"] == "M3H") & (t["variant"] != "M3H[oracle_hops]")]
    if len(arbitrary):
        worst = arbitrary.groupby(["param", "gamma"])[["min_acc", "mean_acc"]].min().reset_index()
        worst = worst.sort_values(["mean_acc", "min_acc"], ascending=False)
        out["M3H worst arbitrary orientation"] = {
            "n_cells": len(worst), "n_perfect": int((worst["min_acc"] == 1.0).sum()),
            "best_cell": worst.iloc[0].to_dict(),
        }
    return out


def run_discrimination(cfg: Config) -> pd.DataFrame:
    """Share of A's forwarding weight on the upper arm as one cost is swept."""
    keep = ("kind", "variant", "param", "gamma", "stable", "split_upper", "nexthop_acc", "mean_cancellation")
    rows = []
    for curve in cfg.curves:
        for cost, scn in curve.points:
            for model in cfg.curve_models:
                for gamma in cfg.curve_gammas:
                    row, _ = evaluate(scn, model, gamma)
                    rows.append({
                        "curve": curve.name, "swept_cost": cost, "crossover": curve.crossover,
                        **{key: row[key] for key in keep},
                    })
    return pd.DataFrame(rows)


@dataclass
class DynamicResult:
    phases: pd.DataFrame
    series: pd.DataFrame
    states: pd.DataFrame
    currents: pd.DataFrame


def settle_round(ok: np.ndarray) -> float:
    """Rounds needed before ok holds for good (1-based). NaN if it does not hold at the end."""
    if len(ok) == 0 or not ok[-1]:
        return math.nan
    bad = np.flatnonzero(~ok)
    return float(bad[-1] + 2) if len(bad) else 1.0


def run_dynamic(scn: Scenario, changes, model: Model, gamma: float, total_rounds: int) -> DynamicResult:
    """Run the local engine through a cost-change schedule, starting from x = 0."""
    dest, ids = scn.dest, _ids(scn, model, gamma)
    nodes = [i for i in range(scn.graph.n) if i != dest]
    x = np.zeros((scn.graph.n, model.k))
    phase_rows, series, states, currents = [], [], [], []
    for p, (start, end, g) in enumerate(phases(scn.graph, changes, total_rounds)):
        K, b = model.coupling(g, dest), model.source(g.n, dest)
        alpha = engine.step_sizes(K, gamma)  # recomputed locally when costs change
        stable = oracle.abscissa(K, gamma) < -STABLE_TOL
        target = oracle.fixed_point(K, b, gamma) if stable else None
        _, opt = metrics.shortest(g, dest)
        length = end - start
        settled = np.zeros(length, dtype=bool)
        hops, loop_rounds, diverged = [], 0, False
        for r in range(length):
            x = engine.step(x, K, b, gamma, alpha)
            if not np.isfinite(x).all() or np.abs(x).max() > DIVERGE_AT:
                diverged = True
                break
            P = readout.power_inflow(K, x)
            W = readout.weights(P)
            nxt, _ = readout.next_hops(W, dest)
            hops.append(tuple(int(v) for v in nxt))
            loop_rounds += any(metrics.follow(nxt, i, dest)[0] == "loop" for i in nodes)
            residual = float(np.linalg.norm(x - target) / np.linalg.norm(target)) if stable else math.nan
            settled[r] = residual < SETTLE_TOL
            if r < LOG_DENSE or r % LOG_EVERY == 0 or r == length - 1:
                t = start + r
                series.append({
                    **ids, "phase": p, "round": t, "residual": residual,
                    "split_upper": metrics.split_share(W, *scn.split) if scn.split else math.nan,
                    "n_correct": sum(hops[-1][i] in opt[i] for i in nodes),
                })
                for i in range(g.n):
                    states.append({
                        **ids, "round": t, "node": g.name(i), "x0": x[i, 0],
                        "x1": x[i, 1] if model.k == 2 else math.nan,
                        "magnitude": float(np.linalg.norm(x[i])),
                    })
                    for j in g.neighbors(i):
                        currents.append({**ids, "round": t, "i": g.name(i), "j": g.name(j), "inflow": P[i, j]})
        done = bool(hops) and not diverged
        phase_rows.append({
            **ids, "phase": p, "start": start, "rounds": length, "stable": stable, "diverged": diverged,
            "settle_rounds": settle_round(settled) if done else math.nan,
            "route_settle_rounds": settle_round(np.array([h == hops[-1] for h in hops])) if done else math.nan,
            "loop_rounds": loop_rounds,
            "final_nexthop_acc": float(np.mean([hops[-1][i] in opt[i] for i in nodes])) if done else math.nan,
            "alpha_mean": float(alpha.mean()),
        })
        if diverged:
            break
    return DynamicResult(pd.DataFrame(phase_rows), pd.DataFrame(series), pd.DataFrame(states), pd.DataFrame(currents))


def run_all_dynamic(cfg: Config) -> DynamicResult:
    runs = [
        run_dynamic(case.scenario, case.changes, model, gamma, case.total_rounds)
        for case in cfg.dynamic_cases
        for model in cfg.dynamic_models
        for gamma in cfg.dynamic_gammas
    ]
    return DynamicResult(*(
        pd.concat([getattr(run, name) for run in runs], ignore_index=True)
        for name in ("phases", "series", "states", "currents")
    ))


def _clean(obj):
    """Make an object JSON-safe: numpy scalars to Python, non-finite floats to null."""
    if isinstance(obj, dict):
        return {str(key): _clean(value) for key, value in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(value) for value in obj]
    if isinstance(obj, np.generic):
        obj = obj.item()
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj


def _git(*args: str) -> str | None:
    try:
        done = subprocess.run(
            ["git", *args], capture_output=True, text=True, check=True, cwd=Path(__file__).parent
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return done.stdout.strip()


def _scenario(scn: Scenario) -> dict:
    return {
        "name": scn.name, "family": scn.family, "dest": scn.dest, "n": scn.graph.n,
        "costs": [[i, j, scn.graph.cost(i, j)] for i, j in scn.graph.edges],
    }


def describe(cfg: Config) -> dict:
    """Everything needed to reproduce the run. No wall-clock values."""
    return {
        "gammas": cfg.gammas,
        "models": [dataclasses.asdict(m) for m in cfg.models],
        "scenarios": [_scenario(s) for s in cfg.scenarios],
        "curves": [
            {"name": c.name, "crossover": c.crossover, "swept_costs": [cost for cost, _ in c.points]}
            for c in cfg.curves
        ],
        "curve_models": [dataclasses.asdict(m) for m in cfg.curve_models],
        "curve_gammas": cfg.curve_gammas,
        "dynamic_cases": [
            {
                "scenario": _scenario(case.scenario), "total_rounds": case.total_rounds,
                "changes": [dataclasses.asdict(change) for change in case.changes],
            }
            for case in cfg.dynamic_cases
        ],
        "dynamic_models": [dataclasses.asdict(m) for m in cfg.dynamic_models],
        "dynamic_gammas": cfg.dynamic_gammas,
        "tolerances": {
            "stable": STABLE_TOL, "settle": SETTLE_TOL, "diverge_at": DIVERGE_AT,
            "log_dense": LOG_DENSE, "log_every": LOG_EVERY,
        },
        "seeds": {"path6_rand": 1, "btree7": 2, "random_orientations": [0, 1, 2, 3, 4]},
        "versions": {
            "python": platform.python_version(),
            **{name: version(name) for name in ("numpy", "networkx", "pandas", "matplotlib")},
        },
        "git_commit": _git("rev-parse", "HEAD"),
    }


def main(out_dir="results/experiment0-smoke", quick: bool = False, make_plots: bool = True) -> dict:
    cfg = quick_config() if quick else full_config()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    static = run_static(cfg.scenarios, cfg.models, cfg.gammas)
    disc = run_discrimination(cfg)
    dyn = run_all_dynamic(cfg)
    transfer = transfer_table(static.metrics)
    summary = _clean({
        "n_static_rows": len(static.metrics),
        "n_unstable_rows": int((~static.metrics["stable"]).sum()),
        "control_checks": control_checks(static),
        "transfer_by_group": summarize_groups(transfer),
    })
    tables = {
        "route_metrics": static.metrics, "static_readouts": static.readouts, "transfer": transfer,
        "discrimination": disc, "dynamic_metrics": dyn.phases, "dynamic_series": dyn.series,
        "states": dyn.states, "currents": dyn.currents,
    }
    for name, table in tables.items():
        table.to_csv(out / f"{name}.csv", index=False, float_format="%.8g")
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "config.json").write_text(json.dumps(_clean(describe(cfg)), indent=2) + "\n")
    if make_plots:
        from . import plots

        plots.make_all(out)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the Experiment 0 smoke test.")
    parser.add_argument("--out", default="results/experiment0-smoke")
    parser.add_argument("--quick", action="store_true", help="small configuration used by the tests")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    result = main(args.out, quick=args.quick, make_plots=not args.no_plots)
    print(json.dumps(result["control_checks"], indent=2))
