"""Figures for the Experiment 0 smoke test, drawn from the CSVs the run wrote."""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

SURFACE, TEXT, MUTED, GRID, MASK = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1", "#d9d8d3"
# One fixed colour and line style per model variant, in every figure.
SERIES = {
    "M3H[id]": ("#2a78d6", "-"),
    "M3S[theta0=-pi/2]": ("#eb6834", "-"),
    "M3S[theta0=0]": ("#1baf7a", "-"),
    "M0": ("#eda100", "--"),
    "M1": ("#e87ba4", "--"),
    "M2": ("#008300", ":"),
}
ORDINAL = ("#86b6ef", "#2a78d6", "#0d366b")  # ordered scale parameter, light to dark
SEQUENTIAL = LinearSegmentedColormap.from_list("phaseflow_blue", ["#cde2fb", "#3987e5", "#0d366b"]).with_extremes(bad=MASK)
ORACLE = "M3H[oracle_hops]"

STYLE = {
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": TEXT, "axes.labelcolor": MUTED, "axes.titlecolor": TEXT,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.edgecolor": GRID,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
    "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True,
    "lines.linewidth": 2, "lines.solid_capstyle": "round", "font.size": 9, "axes.titlesize": 10,
    "legend.frameon": False, "savefig.dpi": 150, "savefig.bbox": "tight",
}


def _save(fig, path: Path) -> Path:
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_discrimination(disc: pd.DataFrame, plots: Path) -> list[Path]:
    """Does node A put more weight on the cheaper arm as one cost is swept?"""
    made = []
    for curve, df in disc.groupby("curve"):
        variants = list(dict.fromkeys(df["variant"]))
        params = sorted(df.loc[df["kind"] != "M2", "param"].unique())
        gammas = sorted(df["gamma"].unique())
        cols = min(4, len(variants))
        rows = math.ceil(len(variants) / cols)
        fig, axes = plt.subplots(
            rows, cols, figsize=(3.3 * cols, 2.6 * rows + 0.9), sharex=True, sharey=True,
            squeeze=False, layout="constrained",
        )
        for ax, variant in zip(axes.flat, variants):
            sub = df[(df["variant"] == variant) & df["stable"]]
            for (param, gamma), line in sub.groupby(["param", "gamma"]):
                line = line.sort_values("swept_cost")
                colour = ORDINAL[min(params.index(param), 2)] if param in params else MUTED
                ax.plot(line["swept_cost"], line["split_upper"], color=colour, linestyle="-" if gamma == gammas[0] else "--")
            ax.axhline(0.5, color=MUTED, linewidth=0.8)
            ax.axvline(df["crossover"].iloc[0], color=MUTED, linewidth=0.8)
            ax.set_title(variant)
            ax.set_ylim(-0.03, 1.03)
        for ax in axes.flat[len(variants):]:
            ax.set_visible(False)
        handles = [Line2D([], [], color=ORDINAL[min(i, 2)], label=f"β or κ = {p:g}") for i, p in enumerate(params)]
        handles += [Line2D([], [], color=MUTED, linestyle="-" if i == 0 else "--", label=f"γ = {g:g}") for i, g in enumerate(gammas)]
        fig.legend(handles=handles, loc="outside right upper")
        fig.supxlabel("swept cost on the lower arm (the arms cost the same at the vertical line)", color=MUTED, fontsize=9)
        fig.supylabel("share of A's forwarding weight on the upper arm", color=MUTED, fontsize=9)
        fig.suptitle(f"{curve}: a correct model is below 0.5 left of the line and above 0.5 right of it")
        made.append(_save(fig, plots / f"discrimination_{curve}.png"))
    return made


def plot_accuracy_heatmaps(transfer: pd.DataFrame, plots: Path) -> list[Path]:
    """Mean next-hop accuracy over the cyclic scenarios, per parameter cell."""
    panels = {}
    arbitrary = transfer[(transfer["kind"] == "M3H") & (transfer["variant"] != ORACLE)]
    for variant in ("M0", "M1", "M2", "M3H[id]"):
        sub = transfer[transfer["variant"] == variant]
        if len(sub):
            panels[variant] = sub.pivot(index="gamma", columns="param", values="mean_acc")
    if len(arbitrary):
        panels["M3H worst arbitrary orientation"] = arbitrary.groupby(["gamma", "param"])["mean_acc"].min().unstack()
    for variant in (ORACLE, "M3S[theta0=-pi/2]", "M3S[theta0=0]"):
        sub = transfer[transfer["variant"] == variant]
        if len(sub):
            panels[variant] = sub.pivot(index="gamma", columns="param", values="mean_acc")
    if not panels:
        return []
    cols = min(4, len(panels))
    rows = math.ceil(len(panels) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(3.2 * cols + 0.8, 2.8 * rows + 0.4), squeeze=False, layout="constrained")
    image = None
    for ax, (title, grid) in zip(axes.flat, panels.items()):
        image = ax.imshow(
            np.ma.masked_invalid(grid.to_numpy(dtype=float)), cmap=SEQUENTIAL, vmin=0.7, vmax=1,
            origin="lower", aspect="auto",
        )
        ax.set_xticks(range(len(grid.columns)), [f"{c:g}" for c in grid.columns])
        ax.set_yticks(range(len(grid.index)), [f"{g:g}" for g in grid.index])
        ax.grid(False)
        ax.set_title(title)
        ax.set_xlabel("β or κ")
        ax.set_ylabel("γ")
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    fig.suptitle("Mean next-hop accuracy over cyclic scenarios (grey: unstable on at least one)")
    fig.colorbar(image, ax=axes, shrink=0.8, label="share of nodes forwarding along a shortest path\n(colour scale starts at 0.7)")
    return [_save(fig, plots / "accuracy_heatmaps.png")]


def plot_orientation_spread(metrics: pd.DataFrame, plots: Path) -> list[Path]:
    """M3-H next-hop accuracy per cyclic scenario, across orientations and parameters."""
    df = metrics[(metrics["kind"] == "M3H") & (metrics["family"] != "tree") & metrics["stable"]]
    if df.empty:
        return []
    names = list(dict.fromkeys(df["scenario"]))
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(1.1 * len(names) + 2, 3.8), layout="constrained")
    groups = (
        ("arbitrary orientations", df[df["variant"] != ORACLE], "#2a78d6"),
        ("oracle_hops (diagnostic)", df[df["variant"] == ORACLE], "#eb6834"),
    )
    for offset, (label, sub, colour) in zip((-0.18, 0.18), groups):
        if sub.empty:
            continue
        x = sub["scenario"].map(names.index) + offset + rng.uniform(-0.12, 0.12, len(sub))
        ax.scatter(x, sub["nexthop_acc"], s=28, color=colour, alpha=0.55, edgecolors=SURFACE, linewidths=0.8, label=label)
    ax.set_xticks(range(len(names)), names, rotation=30, ha="right")
    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel("next-hop accuracy")
    ax.set_title("M3-H: one dot per orientation × β × γ")
    fig.legend(loc="outside lower center", ncol=2)
    return [_save(fig, plots / "m3h_orientation_spread.png")]


def plot_cancellation(metrics: pd.DataFrame, plots: Path) -> list[Path]:
    """Does more destructive interference go with better or worse forwarding?"""
    df = metrics[(metrics["family"] != "tree") & metrics["stable"]]
    fig, ax = plt.subplots(figsize=(5.6, 4.2), layout="constrained")
    drawn = 0
    for variant in ("M3H[id]", "M3S[theta0=-pi/2]", "M3S[theta0=0]"):
        sub = df[df["variant"] == variant]
        if sub.empty:
            continue
        ax.scatter(
            sub["mean_cancellation"], sub["nexthop_acc"], s=28, color=SERIES[variant][0], alpha=0.6,
            edgecolors=SURFACE, linewidths=0.8, label=variant,
        )
        drawn += 1
    if not drawn:
        plt.close(fig)
        return []
    ax.set_xlabel("mean cancellation at a node (0: contributions add, 1: they cancel)")
    ax.set_ylabel("next-hop accuracy")
    ax.set_ylim(-0.03, 1.03)
    ax.set_title("Interference against forwarding quality, one dot per cyclic cell")
    fig.legend(loc="outside lower center", ncol=3)
    return [_save(fig, plots / "cancellation_vs_accuracy.png")]


def plot_dynamic(series: pd.DataFrame, phases: pd.DataFrame, plots: Path) -> list[Path]:
    """Distance to the fixed point and A's split over time, through the cost changes."""
    made = []
    for (scenario, gamma), df in series.groupby(["scenario", "gamma"]):
        fig, (top, bottom) = plt.subplots(2, 1, figsize=(7.5, 5.8), sharex=True, layout="constrained")
        for variant, (colour, style) in SERIES.items():
            sub = df[df["variant"] == variant].sort_values("round")
            if sub.empty:
                continue
            if sub["residual"].notna().any():
                top.semilogy(sub["round"], sub["residual"].clip(lower=1e-16), color=colour, linestyle=style)
            bottom.plot(sub["round"], sub["split_upper"], color=colour, linestyle=style, label=variant)
        starts = phases[(phases["scenario"] == scenario) & (phases["start"] > 0)]["start"].unique()
        for ax in (top, bottom):
            for start in starts:
                ax.axvline(start, color=MUTED, linewidth=0.8)
        top.set_ylabel("relative distance to the\noracle fixed point")
        bottom.set_ylabel("share of A's weight\non the upper arm")
        bottom.set_xlabel("round (vertical lines: cost raised, then restored)")
        bottom.set_ylim(-0.03, 1.03)
        fig.legend(loc="outside lower center", ncol=3)
        fig.suptitle(f"{scenario}, γ = {gamma:g}: unstable models have no distance curve")
        made.append(_save(fig, plots / f"dynamic_{scenario}_gamma{gamma:g}.png"))
    return made


def make_all(out: Path) -> list[Path]:
    out = Path(out)
    plots = out / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(STYLE):
        made = plot_discrimination(pd.read_csv(out / "discrimination.csv"), plots)
        made += plot_accuracy_heatmaps(pd.read_csv(out / "transfer.csv"), plots)
        metrics = pd.read_csv(out / "route_metrics.csv")
        made += plot_orientation_spread(metrics, plots)
        made += plot_cancellation(metrics, plots)
        made += plot_dynamic(pd.read_csv(out / "dynamic_series.csv"), pd.read_csv(out / "dynamic_metrics.csv"), plots)
    return made
