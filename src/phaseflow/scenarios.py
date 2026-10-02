"""Scenarios and the parameter grid, declared before the run (spec sections 7-9)."""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .graph import Graph
from .models import Model
from .perturbations import CostChange
from .topology import binary_tree, diamond, lopsided_diamond, path, random_costs, ring

GAMMAS = (0.1, 0.3, 1.0, 2.0, 3.5)
SCALES = (0.25, 0.5, 1.0, 2.0, 4.0)  # kappa for M0/M1, beta for M3-H/M3-S
ARBITRARY_ORIENTATIONS = ("id", "reverse", "random:0", "random:1", "random:2", "random:3", "random:4")
ORIENTATIONS = ARBITRARY_ORIENTATIONS + ("oracle_hops",)
THETA0S = (-math.pi / 2, 0.0)
ROUNDS_PER_PHASE = 2500

DIAMOND_CYCLE, DIAMOND_SPLIT = (1, 2, 4, 6, 5, 3), (1, 2, 3)
LOP_CYCLE, LOP_SPLIT = (1, 2, 5, 4, 3), (1, 2, 3)


@dataclass(frozen=True)
class Scenario:
    name: str
    family: str  # "tree" | "diamond" | "lopsided" | "ring"
    graph: Graph
    dest: int
    cycle: tuple[int, ...] = ()  # node order around the single cycle, if any
    split: tuple[int, ...] = ()  # (node, upper neighbour, lower neighbour)


@dataclass(frozen=True)
class Curve:
    name: str
    crossover: float  # swept cost at which the two arms cost the same
    points: tuple[tuple[float, Scenario], ...]


@dataclass(frozen=True)
class DynamicCase:
    scenario: Scenario
    changes: tuple[CostChange, ...]
    total_rounds: int


@dataclass(frozen=True)
class Config:
    scenarios: tuple[Scenario, ...]
    models: tuple[Model, ...]
    gammas: tuple[float, ...]
    curves: tuple[Curve, ...]
    curve_models: tuple[Model, ...]
    curve_gammas: tuple[float, ...]
    dynamic_cases: tuple[DynamicCase, ...]
    dynamic_models: tuple[Model, ...]
    dynamic_gammas: tuple[float, ...]


def _diamond(name: str, g: Graph) -> Scenario:
    return Scenario(name, "diamond", g, 6, DIAMOND_CYCLE, DIAMOND_SPLIT)


def _lopsided(name: str, lower_cost: float) -> Scenario:
    g = lopsided_diamond()
    for edge in ((1, 3), (3, 4), (4, 5)):
        g = g.with_cost(*edge, lower_cost)
    return Scenario(name, "lopsided", g, 5, LOP_CYCLE, LOP_SPLIT)


def _ring(name: str, n: int, heavy: bool) -> Scenario:
    g = ring(n).with_cost(1, 2, 4.0) if heavy else ring(n)
    return Scenario(name, "ring", g, 0, tuple(range(n)))


def static_scenarios() -> tuple[Scenario, ...]:
    d = diamond()
    tree = random_costs(binary_tree(2), 2)
    return (
        Scenario("path6", "tree", path(6), 0),
        Scenario("path6_rand", "tree", random_costs(path(6), 1), 0),
        Scenario("btree7_root", "tree", tree, 0),
        Scenario("btree7_leaf", "tree", tree, 6),
        _diamond("diamond_equal", d),
        _diamond("diamond_lower+0.2", d.with_cost(3, 5, 1.2)),
        _diamond("diamond_lower+1.0", d.with_cost(3, 5, 2.0)),
        _diamond("diamond_upper+1.0", d.with_cost(2, 4, 2.0)),
        _lopsided("lop_unit", 1.0),
        _lopsided("lop_lower_cheap", 0.5),
        _lopsided("lop_near_tie", 0.7),
        _ring("ring8_unit", 8, False),
        _ring("ring8_heavy", 8, True),
        _ring("ring7_unit", 7, False),
        _ring("ring7_heavy", 7, True),
    )


def model_grid(scales=SCALES, orientations=ORIENTATIONS, theta0s=THETA0S) -> tuple[Model, ...]:
    out = [Model("M0", kappa=s) for s in scales]
    out += [Model("M1", kappa=s, kappa2=2 * s) for s in scales]
    out += [Model("M2")]
    out += [Model("M3H", beta=s, orientation=o) for o in orientations for s in scales]
    out += [Model("M3S", beta=s, theta0=t) for t in theta0s for s in scales]
    return tuple(out)


def curves(n_diamond: int = 21, n_lopsided: int = 25) -> tuple[Curve, ...]:
    d = diamond()
    dia = tuple(
        (float(c), _diamond(f"diamond_CE={c:.2f}", d.with_cost(3, 5, float(c))))
        for c in np.linspace(0.5, 2.5, n_diamond)
    )
    lop = tuple(
        (float(c), _lopsided(f"lop_lower={c:.2f}", float(c))) for c in np.linspace(0.3, 1.5, n_lopsided)
    )
    return (Curve("diamond", 1.0, dia), Curve("lopsided", 2 / 3, lop))


def dynamic_cases(rounds: int = ROUNDS_PER_PHASE) -> tuple[DynamicCase, ...]:
    d = _diamond("diamond", diamond())
    lop = _lopsided("lopsided", 1.0)
    return (
        DynamicCase(d, (CostChange(rounds, (2, 4), 3.0), CostChange(2 * rounds, (2, 4), 1.0)), 3 * rounds),
        DynamicCase(lop, (CostChange(rounds, (2, 5), 3.0), CostChange(2 * rounds, (2, 5), 1.0)), 3 * rounds),
    )


def dynamic_models() -> tuple[Model, ...]:
    return (
        Model("M0", kappa=1.0),
        Model("M1", kappa=1.0, kappa2=2.0),
        Model("M2"),
        Model("M3H", beta=1.0),
        Model("M3S", beta=1.0),
        Model("M3S", beta=1.0, theta0=0.0),
    )


def full_config() -> Config:
    return Config(
        scenarios=static_scenarios(),
        models=model_grid(),
        gammas=GAMMAS,
        curves=curves(),
        curve_models=model_grid((0.5, 1.0, 2.0), ("id", "random:0", "oracle_hops")),
        curve_gammas=(0.3, 3.5),
        dynamic_cases=dynamic_cases(),
        dynamic_models=dynamic_models(),
        dynamic_gammas=(0.3, 1.0, 3.5),
    )


def quick_config() -> Config:
    """A small configuration that exercises every code path; used by the tests."""
    by_name = {s.name: s for s in static_scenarios()}
    return Config(
        scenarios=tuple(by_name[n] for n in ("path6_rand", "diamond_lower+1.0", "lop_unit", "ring7_heavy")),
        models=model_grid((1.0,), ("id", "oracle_hops")),
        gammas=(0.3, 3.5),
        curves=curves(5, 5),
        curve_models=model_grid((1.0,), ("id",)),
        curve_gammas=(3.5,),
        dynamic_cases=dynamic_cases(200),
        dynamic_models=dynamic_models(),
        dynamic_gammas=(3.5,),
    )
