# Experiment 0 Smoke Test Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the Experiment 0 smoke test: five routing-field models on tree, diamond and ring graphs, with an exact oracle, a shared forwarding readout, one cost perturbation, and reproducible CSV/JSON/plot outputs.

**Architecture:** Every model is the same linear local system `dx_i/dt = -γ x_i + Σ_j K_ij x_j + b_i` over real `k`-vectors; a model is only a rule for the `k×k` block `K_ij`. One engine iterates it locally, one oracle solves it exactly, one readout (`P_i←j = 2 x_iᵀ K_ij x_j`) turns it into forwarding weights. The experiment runner sweeps a predeclared grid and writes flat files; plots are drawn from those files.

**Tech Stack:** Python ≥ 3.12, `uv`, `numpy`, `networkx`, `pandas`, `matplotlib`, `pytest`.

**Spec:** `docs/superpowers/specs/2026-10-02-experiment0-smoke-design.md`

## Global Constraints

- Exploratory only. No code path selects a "best" cell and reports it alone; the whole grid is written out.
- M3-H and M3-S are separate hypotheses. No output combines them into one figure of merit.
- The oracle (exact fixed point, spectra) is never an input to `engine.step` or `engine.step_sizes`.
- `a = 1`. Costs are symmetric and strictly positive.
- `K[i, j]` always means "what `i` receives from `j`". `P[i, j]` always means "power fed into `i` by `j`".
- Forwarding moves upstream: node `i` forwards to `argmax_j max(P[i, j], 0)`.
- Grid: `γ ∈ {0.1, 0.3, 1.0, 2.0, 3.5}`; scale `κ` or `β ∈ {0.25, 0.5, 1, 2, 4}`; M1 uses `(κ, 2κ)`; M3-H orientations `id, reverse, random:0..4, oracle_hops`; M3-S `θ₀ ∈ {-π/2, 0}`.
- `oracle_hops` is a diagnostic and is excluded from M3-H summaries.
- Unstable cells (spectral abscissa ≥ `-1e-9`) get NaN readout metrics.
- Runs are deterministic: no wall-clock values in any output, all randomness seeded.
- Commit messages end with the two attribution lines used in the first commit of this repository.

## Review Focus

1. **Unstable cell.** The linear solve still returns a "fixed point" below the stability threshold. Expected: `stable = False`, readout metrics NaN, no inflow matrix stored. Test in Task 7.
2. **Dark or zero-field node.** A node whose incoming contributions cancel, or an all-zero state at round 0. Expected: no exception, next hop `-1`, cancellation NaN when there is nothing to cancel. Test in Task 4.
3. **Exact ties.** `diamond_equal` gives node A two equal weights. Expected: lowest id chosen, tie recorded, counted correct. Test in Task 4 and Task 7.
4. **Stochastic policy with a closed loop that excludes the destination.** `I - Q` is singular. Expected: delivery probability 0, expected stretch `inf`, no `LinAlgError`. Test in Task 5.
5. **Divergence during the dynamic run.** Expected: the run stops at `1e12`, the phase is flagged `diverged`, later phases are not run, logged values are finite. Tests in Task 3 and Task 8.

---

## File map

| File | Responsibility |
|---|---|
| `pyproject.toml` | package and dependencies |
| `src/phaseflow/graph.py` | undirected weighted graph, immutable cost updates |
| `src/phaseflow/topology.py` | path, binary tree, ring, diamond, lopsided diamond, seeded costs |
| `src/phaseflow/models.py` | `Model`: block rule, orientation conventions, source; complex helpers |
| `src/phaseflow/oracle.py` | exact fixed point, spectral abscissa, iteration spectral radius |
| `src/phaseflow/engine.py` | local step rule and synchronous iteration |
| `src/phaseflow/readout.py` | power inflow, weights, cancellation, next hops, stochastic policy |
| `src/phaseflow/metrics.py` | shortest-path reference, routes, stretch, loops, stochastic metrics, circulation |
| `src/phaseflow/perturbations.py` | cost-change schedule |
| `src/phaseflow/scenarios.py` | scenarios, grids, curves, dynamic cases, `Config` |
| `src/phaseflow/experiment0.py` | static sweep, curves, dynamic run, summaries, writers, CLI |
| `src/phaseflow/plots.py` | figures from the written CSVs |
| `tests/test_*.py` | one test module per source module |

---

### Task 1: Project scaffold, graph and topologies

**Files:**
- Create: `pyproject.toml`, `src/phaseflow/__init__.py`, `src/phaseflow/graph.py`, `src/phaseflow/topology.py`
- Test: `tests/test_graph_topology.py`

**Interfaces:**
- Produces: `Graph(n: int, costs: dict[tuple[int, int], float], names: tuple[str, ...] = ())` with `.edges -> list[tuple[int, int]]`, `.cost(i, j) -> float`, `.neighbors(i) -> list[int]`, `.with_cost(i, j, cost) -> Graph`, `.name(i) -> str`, `.to_networkx() -> nx.Graph` (edge attribute `weight`).
- Produces: `topology.path(n)`, `topology.ring(n)`, `topology.binary_tree(depth)`, `topology.diamond()`, `topology.lopsided_diamond()`, `topology.random_costs(g, seed, lo=0.5, hi=2.0)`, all returning `Graph`.
- Node ids — diamond: `S0 A1 B2 C3 D4 E5 T6`, arms `A-B-D-T` and `A-C-E-T`. Lopsided: `S0 A1 B2 C3 E4 T5`, arms `A-B-T` and `A-C-E-T`.

- [ ] **Step 1: Write `pyproject.toml` and the empty package**

```toml
[project]
name = "phaseflow"
version = "0.1.0"
description = "Exploratory study of interference-based routing fields on classical networks"
requires-python = ">=3.12"
dependencies = ["numpy>=2.0", "networkx>=3.3", "pandas>=2.2", "matplotlib>=3.9"]

[dependency-groups]
dev = ["pytest>=8"]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.hatch.build.targets.wheel]
packages = ["src/phaseflow"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

`src/phaseflow/__init__.py`:

```python
"""PhaseFlow: exploratory interference-based routing fields."""
```

- [ ] **Step 2: Write the failing tests**

`tests/test_graph_topology.py`:

```python
import networkx as nx
import pytest

from phaseflow import topology
from phaseflow.graph import Graph


def test_graph_rejects_bad_edges_and_costs():
    with pytest.raises(ValueError):
        Graph(2, {(1, 0): 1.0})
    with pytest.raises(ValueError):
        Graph(2, {(0, 1): 0.0})
    with pytest.raises(ValueError):
        Graph(2, {(0, 1): float("nan")})


def test_cost_is_symmetric_and_with_cost_does_not_mutate():
    g = topology.path(3)
    h = g.with_cost(2, 1, 2.5)
    assert h.cost(1, 2) == h.cost(2, 1) == 2.5
    assert g.cost(1, 2) == 1.0
    with pytest.raises(KeyError):
        g.with_cost(0, 2, 1.0)


def test_topology_shapes():
    assert len(topology.path(6).edges) == 5
    tree = topology.binary_tree(2)
    assert tree.n == 7 and len(tree.edges) == 6
    assert len(topology.ring(8).edges) == 8
    d = topology.diamond()
    assert d.n == 7 and len(d.edges) == 7
    assert d.neighbors(1) == [0, 2, 3] and d.name(6) == "T"
    lop = topology.lopsided_diamond()
    assert lop.n == 6 and len(lop.edges) == 6 and lop.neighbors(5) == [2, 4]


def test_cycle_parity():
    assert nx.is_bipartite(topology.diamond().to_networkx())
    assert not nx.is_bipartite(topology.lopsided_diamond().to_networkx())
    assert nx.is_tree(topology.binary_tree(2).to_networkx())


def test_random_costs_are_seeded_and_bounded():
    a = topology.random_costs(topology.path(6), 1)
    b = topology.random_costs(topology.path(6), 1)
    assert a.costs == b.costs
    assert all(0.5 <= c <= 2.0 for c in a.costs.values())
```

- [ ] **Step 3: Run the tests and confirm they fail**

Run: `uv run pytest tests/test_graph_topology.py -q`
Expected: collection error, `ModuleNotFoundError: No module named 'phaseflow.graph'` (or `phaseflow.topology`).

- [ ] **Step 4: Implement `graph.py` and `topology.py`**

`src/phaseflow/graph.py`:

```python
"""Undirected weighted graph with positive symmetric link costs."""
from __future__ import annotations

from dataclasses import dataclass

import networkx as nx


@dataclass(frozen=True)
class Graph:
    n: int
    costs: dict[tuple[int, int], float]  # key (i, j) with i < j
    names: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for (i, j), c in self.costs.items():
            if not 0 <= i < j < self.n:
                raise ValueError(f"edge key must satisfy 0 <= i < j < n, got {(i, j)}")
            if not c > 0:
                raise ValueError(f"cost must be positive, got {c} on {(i, j)}")
        if self.names and len(self.names) != self.n:
            raise ValueError("names must have one entry per node")

    @property
    def edges(self) -> list[tuple[int, int]]:
        return sorted(self.costs)

    def cost(self, i: int, j: int) -> float:
        return self.costs[(min(i, j), max(i, j))]

    def neighbors(self, i: int) -> list[int]:
        return sorted(b if a == i else a for a, b in self.costs if i in (a, b))

    def with_cost(self, i: int, j: int, cost: float) -> Graph:
        key = (min(i, j), max(i, j))
        if key not in self.costs:
            raise KeyError(f"no edge {key}")
        return Graph(self.n, {**self.costs, key: cost}, self.names)

    def name(self, i: int) -> str:
        return self.names[i] if self.names else str(i)

    def to_networkx(self) -> nx.Graph:
        g = nx.Graph()
        g.add_nodes_from(range(self.n))
        g.add_weighted_edges_from((i, j, c) for (i, j), c in self.costs.items())
        return g
```

`src/phaseflow/topology.py`:

```python
"""Topologies for the smoke test. All costs are 1 unless changed afterwards."""
from __future__ import annotations

import numpy as np

from .graph import Graph


def _unit(n: int, pairs: list[tuple[int, int]], names: str = "") -> Graph:
    return Graph(n, {(min(i, j), max(i, j)): 1.0 for i, j in pairs}, tuple(names))


def path(n: int) -> Graph:
    return _unit(n, [(i, i + 1) for i in range(n - 1)])


def ring(n: int) -> Graph:
    return _unit(n, [(i, (i + 1) % n) for i in range(n)])


def binary_tree(depth: int) -> Graph:
    n = 2 ** (depth + 1) - 1
    return _unit(n, [((child - 1) // 2, child) for child in range(1, n)])


def diamond() -> Graph:
    """S-A, then arms A-B-D-T (upper) and A-C-E-T (lower). One 6-cycle."""
    return _unit(7, [(0, 1), (1, 2), (2, 4), (4, 6), (1, 3), (3, 5), (5, 6)], "SABCDET")


def lopsided_diamond() -> Graph:
    """S-A, then arms A-B-T (upper) and A-C-E-T (lower). One 5-cycle."""
    return _unit(6, [(0, 1), (1, 2), (2, 5), (1, 3), (3, 4), (4, 5)], "SABCET")


def random_costs(g: Graph, seed: int, lo: float = 0.5, hi: float = 2.0) -> Graph:
    rng = np.random.default_rng(seed)
    return Graph(g.n, {e: float(rng.uniform(lo, hi)) for e in g.edges}, g.names)
```

- [ ] **Step 5: Run the tests and confirm they pass**

Run: `uv run pytest tests/test_graph_topology.py -q`
Expected: `5 passed`.

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml uv.lock src/phaseflow/__init__.py src/phaseflow/graph.py src/phaseflow/topology.py tests/test_graph_topology.py
git commit -m "Add graph type and smoke-test topologies"
```

---

### Task 2: Models and oracle

**Files:**
- Create: `src/phaseflow/models.py`, `src/phaseflow/oracle.py`
- Test: `tests/test_models.py`

**Interfaces:**
- Consumes: `Graph`, `topology`.
- Produces: `Model(kind, kappa=0.0, kappa2=0.0, beta=0.0, theta0=-π/2, orientation="id", a=1.0)` with `.k -> int`, `.is_complex -> bool`, `.variant -> str`, `.param -> float`, `.coupling(g, dest) -> ndarray (n, n, k, k)`, `.source(n, dest) -> ndarray (n, k)`.
- `variant` strings: `"M0"`, `"M1"`, `"M2"`, `"M3H[<orientation>]"`, `"M3S[theta0=-pi/2]"`, `"M3S[theta0=0]"`. `param` is `kappa` for M0/M1, `beta` for M3H/M3S, `0.0` for M2.
- Produces: `orientation_signs(g, convention, dest) -> dict[tuple[int, int], int]`, `rot(theta)`, `to_complex(x)`, `complex_matrix(K)`.
- Produces: `oracle.flat(K) -> (nk, nk)`, `oracle.fixed_point(K, b, gamma) -> (n, k)`, `oracle.abscissa(K, gamma) -> float`, `oracle.iteration_radius(K, gamma, alpha) -> float`.

- [ ] **Step 1: Write the failing tests**

`tests/test_models.py`:

```python
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
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run pytest tests/test_models.py -q`
Expected: `ModuleNotFoundError: No module named 'phaseflow.oracle'` (or `phaseflow.models`).

- [ ] **Step 3: Implement `models.py` and `oracle.py`**

`src/phaseflow/models.py`:

```python
"""The five models. Each is a rule for the real k-by-k block K_ij on an edge.

K[i, j] is what node i receives from node j. A complex number p + iq is the
block [[p, -q], [q, p]], so complex models are two reals per node.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import networkx as nx
import numpy as np

from .graph import Graph

KINDS = ("M0", "M1", "M2", "M3H", "M3S")


def rot(theta: float) -> np.ndarray:
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, -s], [s, c]])


def orientation_signs(g: Graph, convention: str, dest: int) -> dict[tuple[int, int], int]:
    """sigma_ij for every edge key (i, j) with i < j. sigma_ji = -sigma_ij."""
    edges = g.edges
    if convention == "id":
        return {e: 1 for e in edges}
    if convention == "reverse":
        return {e: -1 for e in edges}
    if convention.startswith("random:"):
        rng = np.random.default_rng(int(convention.split(":")[1]))
        return {e: int(s) for e, s in zip(edges, rng.choice([-1, 1], size=len(edges)))}
    if convention == "oracle_hops":
        # Diagnostic only: uses global knowledge of the destination.
        hops = nx.single_source_shortest_path_length(g.to_networkx(), dest)
        return {(i, j): -1 if hops[i] < hops[j] else 1 for i, j in edges}
    raise ValueError(f"unknown orientation convention {convention!r}")


@dataclass(frozen=True)
class Model:
    kind: str
    kappa: float = 0.0  # M0, M1: attenuation scale
    kappa2: float = 0.0  # M1: second attenuation scale
    beta: float = 0.0  # M3H, M3S: phase per unit cost
    theta0: float = -math.pi / 2  # M3S: per-hop phase offset
    orientation: str = "id"  # M3H: edge sign convention
    a: float = 1.0

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise ValueError(f"unknown model kind {self.kind!r}")

    @property
    def k(self) -> int:
        return 1 if self.kind == "M0" else 2

    @property
    def is_complex(self) -> bool:
        return self.kind in ("M2", "M3H", "M3S")

    @property
    def variant(self) -> str:
        if self.kind == "M3H":
            return f"M3H[{self.orientation}]"
        if self.kind == "M3S":
            return f"M3S[theta0={'0' if self.theta0 == 0 else '-pi/2'}]"
        return self.kind

    @property
    def param(self) -> float:
        if self.kind in ("M0", "M1"):
            return self.kappa
        return self.beta if self.kind in ("M3H", "M3S") else 0.0

    def block(self, cost: float, sigma: int) -> np.ndarray:
        """K_ij for one direction of an edge. sigma is the M3-H sign for that direction."""
        if self.kind == "M0":
            return np.array([[self.a * math.exp(-self.kappa * cost)]])
        if self.kind == "M1":
            return self.a * np.diag([math.exp(-self.kappa * cost), math.exp(-self.kappa2 * cost)])
        if self.kind == "M2":
            return self.a * rot(-math.pi / 2)
        if self.kind == "M3H":
            return self.a * rot(-math.pi / 2 + sigma * self.beta * cost)
        return self.a * rot(self.theta0 + self.beta * cost)

    def coupling(self, g: Graph, dest: int) -> np.ndarray:
        K = np.zeros((g.n, g.n, self.k, self.k))
        signs = orientation_signs(g, self.orientation, dest) if self.kind == "M3H" else {}
        for i, j in g.edges:
            sigma = signs.get((i, j), 1)
            K[i, j] = self.block(g.cost(i, j), sigma)
            K[j, i] = self.block(g.cost(i, j), -sigma)
        return K

    def source(self, n: int, dest: int) -> np.ndarray:
        b = np.zeros((n, self.k))
        b[dest] = (1.0, 0.0) if self.is_complex else 1.0
        return b


def to_complex(x: np.ndarray) -> np.ndarray:
    """(n, 2) real state -> complex vector."""
    return x[:, 0] + 1j * x[:, 1]


def complex_matrix(K: np.ndarray) -> np.ndarray:
    """(n, n, 2, 2) rotation-scaling blocks -> complex (n, n) matrix."""
    return K[:, :, 0, 0] + 1j * K[:, :, 1, 0]
```

`src/phaseflow/oracle.py`:

```python
"""Exact solutions used only to validate the local dynamics.

Nothing here may be passed to the engine as an input.
"""
from __future__ import annotations

import numpy as np


def flat(K: np.ndarray) -> np.ndarray:
    n, _, k, _ = K.shape
    return K.transpose(0, 2, 1, 3).reshape(n * k, n * k)


def fixed_point(K: np.ndarray, b: np.ndarray, gamma: float) -> np.ndarray:
    n, _, k, _ = K.shape
    return np.linalg.solve(gamma * np.eye(n * k) - flat(K), b.reshape(n * k)).reshape(n, k)


def abscissa(K: np.ndarray, gamma: float) -> float:
    """Largest real part of the spectrum of -gamma*I + K. Negative means stable."""
    return float(np.linalg.eigvals(flat(K)).real.max() - gamma)


def iteration_radius(K: np.ndarray, gamma: float, alpha: np.ndarray) -> float:
    """Spectral radius of the local iteration matrix. Below 1 means convergent."""
    n, _, k, _ = K.shape
    M = flat(K) - gamma * np.eye(n * k)
    G = np.eye(n * k) + np.repeat(alpha, k)[:, None] * M
    return float(np.abs(np.linalg.eigvals(G)).max())
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `uv run pytest tests/test_models.py -q`
Expected: `10 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/phaseflow/models.py src/phaseflow/oracle.py tests/test_models.py
git commit -m "Add the five coupling models and the exact oracle"
```

---

### Task 3: Local engine

**Files:**
- Create: `src/phaseflow/engine.py`
- Test: `tests/test_engine.py`

**Interfaces:**
- Consumes: `Model.coupling`, `Model.source`, `oracle.fixed_point`, `oracle.abscissa`, `oracle.iteration_radius`.
- Produces: `engine.step_sizes(K, gamma) -> ndarray (n,)`, `engine.step(x, K, b, gamma, alpha) -> ndarray (n, k)`, `engine.run(x, K, b, gamma, rounds, diverge_at=1e12) -> tuple[ndarray, int, bool]` returning `(state, rounds_done, diverged)`.

- [ ] **Step 1: Write the failing tests**

`tests/test_engine.py`:

```python
import numpy as np

from phaseflow import engine, oracle, topology
from phaseflow.models import Model

G = topology.random_costs(topology.lopsided_diamond(), 3)
DEST = 5


def test_converges_to_oracle_fixed_point():
    cases = (
        (Model("M3H", beta=1.0), 0.5),
        (Model("M3S", beta=1.0), 3.5),
        (Model("M0", kappa=2.0), 1.0),
        (Model("M1", kappa=2.0, kappa2=4.0), 1.0),
    )
    for m, gamma in cases:
        K, b = m.coupling(G, DEST), m.source(G.n, DEST)
        x, rounds, diverged = engine.run(np.zeros((G.n, m.k)), K, b, gamma, 4000)
        assert not diverged and rounds == 4000
        assert np.allclose(x, oracle.fixed_point(K, b, gamma), atol=1e-9)


def test_step_size_uses_only_local_coupling_strength():
    K = Model("M2").coupling(G, DEST)
    degree = np.array([len(G.neighbors(i)) for i in range(G.n)])
    assert np.allclose(engine.step_sizes(K, 0.5), 0.5 / (0.25 + degree**2))


def test_update_is_local():
    m = Model("M3S", beta=1.0)
    K, b = m.coupling(G, DEST), m.source(G.n, DEST)
    alpha = engine.step_sizes(K, 3.5)
    x = np.random.default_rng(0).normal(size=(G.n, 2))
    y = x.copy()
    y[4] += 10.0  # node 4 (E) has neighbours 3 (C) and 5 (T) only
    a, c = engine.step(x, K, b, 3.5, alpha), engine.step(y, K, b, 3.5, alpha)
    changed = {i for i in range(G.n) if not np.allclose(a[i], c[i])}
    assert changed == {3, 4, 5}


def test_hermitian_iteration_converges_at_small_damping_on_irregular_graph():
    for beta in (0.0, 1.0):
        K = Model("M3H", beta=beta, orientation="random:1").coupling(G, DEST)
        for gamma in (0.05, 0.3, 2.0):
            assert oracle.iteration_radius(K, gamma, engine.step_sizes(K, gamma)) < 1.0


def test_divergence_below_threshold_is_predicted_and_detected():
    m = Model("M3S", beta=1.0)
    K, b = m.coupling(G, DEST), m.source(G.n, DEST)
    assert oracle.abscissa(K, 0.1) > 0
    assert oracle.iteration_radius(K, 0.1, engine.step_sizes(K, 0.1)) > 1.0
    x, rounds, diverged = engine.run(np.zeros((G.n, 2)), K, b, 0.1, 20000)
    assert diverged and rounds < 20000
    assert np.isfinite(x).all()
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run pytest tests/test_engine.py -q`
Expected: `ImportError: cannot import name 'engine' from 'phaseflow'`.

- [ ] **Step 3: Implement `engine.py`**

```python
"""Local synchronous iteration.

One round: every node sends its k reals to each neighbour, then updates from
what it received. Node i uses only its own state, its own blocks K[i, :] and
its neighbours' states. Nothing from the oracle enters here.
"""
from __future__ import annotations

import numpy as np


def step_sizes(K: np.ndarray, gamma: float) -> np.ndarray:
    """alpha_i = gamma / (gamma^2 + s_i^2), s_i = sum_j ||K_ij||_2 (node-local)."""
    s = np.linalg.norm(K, ord=2, axis=(2, 3)).sum(axis=1)
    return gamma / (gamma**2 + s**2)


def step(x: np.ndarray, K: np.ndarray, b: np.ndarray, gamma: float, alpha: np.ndarray) -> np.ndarray:
    received = np.einsum("ijab,jb->ia", K, x)
    return x + alpha[:, None] * (-gamma * x + received + b)


def run(
    x: np.ndarray, K: np.ndarray, b: np.ndarray, gamma: float, rounds: int, diverge_at: float = 1e12
) -> tuple[np.ndarray, int, bool]:
    """Iterate; stop early and report divergence if the state leaves [-diverge_at, diverge_at]."""
    alpha = step_sizes(K, gamma)
    for t in range(rounds):
        x = step(x, K, b, gamma, alpha)
        if not np.isfinite(x).all() or np.abs(x).max() > diverge_at:
            return x, t + 1, True
    return x, rounds, False
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `uv run pytest tests/test_engine.py -q`
Expected: `5 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/phaseflow/engine.py tests/test_engine.py
git commit -m "Add local synchronous engine with node-local step rule"
```

---

### Task 4: Readout and the convention tests

**Files:**
- Create: `src/phaseflow/readout.py`
- Test: `tests/test_readout.py`

**Interfaces:**
- Consumes: `Model`, `oracle.fixed_point`, topologies.
- Produces: `readout.power_inflow(K, x) -> ndarray (n, n)`, `readout.weights(P) -> ndarray`, `readout.cancellation(K, x, is_complex) -> ndarray (n,)`, `readout.next_hops(W, dest, rtol=1e-9) -> tuple[ndarray[int], ndarray[bool]]` (next hop `-1` for a node with no positive weight, `dest` for the destination), `readout.policy(W, dest) -> ndarray (n, n)` row-stochastic with zero rows for the destination and dead nodes.

These tests are the spec's A1–A5 and tests 1–7 and 10 of its §11. They gate the sweep.

- [ ] **Step 1: Write the failing tests**

`tests/test_readout.py`:

```python
import numpy as np
import pytest

from phaseflow import oracle, readout, topology
from phaseflow.models import Model

ALL = (
    Model("M0", kappa=1.0),
    Model("M1", kappa=1.0, kappa2=2.0),
    Model("M2"),
    Model("M3H", beta=0.7),
    Model("M3S", beta=0.7),
    Model("M3S", beta=0.7, theta0=0.0),
)
LOP = topology.random_costs(topology.lopsided_diamond(), 3)


def solve(m, g, dest, gamma):
    K, b = m.coupling(g, dest), m.source(g.n, dest)
    x = oracle.fixed_point(K, b, gamma)
    return K, b, x, readout.power_inflow(K, x)


def test_a1_two_node_sign_and_forwarding_direction():
    g = topology.path(2)
    for m in ALL:
        for gamma in (1.5, 3.0):  # above every model's threshold on two nodes
            _, _, x, P = solve(m, g, 0, gamma)
            assert P[1, 0] > 0
            assert P[1, 0] == pytest.approx(2 * gamma * (x[1] @ x[1]))
            nxt, _ = readout.next_hops(readout.weights(P), 0)
            assert nxt[1] == 0
    for m in (Model("M2"), Model("M3H", beta=0.7)):
        for gamma in (0.2, 1.0, 3.0):
            _, _, _, P = solve(m, g, 0, gamma)
            assert P[1, 0] > 0 and P[0, 1] == pytest.approx(-P[1, 0])


def test_a2_balance_identity():
    gamma = 3.5
    for m in ALL:
        _, b, x, P = solve(m, LOP, 5, gamma)
        lhs = P.sum(axis=1) + 2 * np.einsum("ia,ia->i", x, b)
        assert np.allclose(lhs, 2 * gamma * np.einsum("ia,ia->i", x, x), atol=1e-12)


def test_inflow_is_a_current_only_for_hermitian_models():
    for m in (Model("M2"), Model("M3H", beta=0.7, orientation="random:1")):
        _, _, _, P = solve(m, LOP, 5, 0.3)
        assert np.allclose(P, -P.T, atol=1e-12)
    _, _, _, P = solve(Model("M3S", beta=0.7), LOP, 5, 3.5)
    assert not np.allclose(P, -P.T, atol=1e-6)
    _, _, _, P = solve(Model("M0", kappa=1.0), LOP, 5, 3.5)
    assert (P >= 0).all()


def test_a3_tree_gauge_invariance_on_observables_not_raw_state():
    g = topology.random_costs(topology.binary_tree(2), 2)
    for dest in (0, 6):
        _, _, x0, P0 = solve(Model("M2"), g, dest, 0.3)
        for orientation in ("id", "reverse", "random:0", "random:1", "oracle_hops"):
            for beta in (0.5, 1.7):
                _, _, x, P = solve(Model("M3H", beta=beta, orientation=orientation), g, dest, 0.3)
                assert np.allclose(np.linalg.norm(x, axis=1), np.linalg.norm(x0, axis=1), atol=1e-12)
                assert np.allclose(P, P0, atol=1e-12)
    _, _, x, _ = solve(Model("M3H", beta=1.7, orientation="id"), g, 0, 0.3)
    _, _, x0, _ = solve(Model("M2"), g, 0, 0.3)
    assert not np.allclose(x, x0, atol=1e-3)  # raw states differ; only observables agree


def test_a4_bipartite_flux_sign_blindness_and_odd_cycle_counterexample():
    def flip_difference(g):
        dest = g.n - 1
        _, _, _, Pa = solve(Model("M3H", beta=0.8, orientation="id"), g, dest, 0.3)
        _, _, _, Pb = solve(Model("M3H", beta=0.8, orientation="reverse"), g, dest, 0.3)
        return np.abs(Pa - Pb).max()

    for g in (topology.ring(8), topology.diamond()):
        assert flip_difference(topology.random_costs(g, 4)) < 1e-12
    for g in (topology.ring(7), topology.lopsided_diamond()):
        assert flip_difference(topology.random_costs(g, 4)) > 1e-3


def test_a5_diamond_split_is_exactly_half_for_m3h():
    for seed in range(3):
        g = topology.random_costs(topology.diamond(), seed)
        for orientation in ("id", "random:0", "oracle_hops"):
            for beta in (0.5, 1.0, 2.0):
                for gamma in (0.1, 1.0):
                    _, _, _, P = solve(Model("M3H", beta=beta, orientation=orientation), g, 6, gamma)
                    W = readout.weights(P)
                    assert W[1, 2] > 0
                    assert W[1, 2] == pytest.approx(W[1, 3], rel=1e-9)


def test_cancellation_is_zero_for_real_models_and_positive_for_coherent_ones():
    for m in (Model("M0", kappa=1.0), Model("M1", kappa=1.0, kappa2=2.0)):
        K, _, x, _ = solve(m, LOP, 5, 3.5)
        assert np.allclose(readout.cancellation(K, x, m.is_complex), 0.0, atol=1e-12)
    m = Model("M3S", beta=1.0)
    K, _, x, _ = solve(m, LOP, 5, 3.5)
    assert readout.cancellation(K, x, m.is_complex)[1] > 1e-3  # node A has three neighbours


def test_next_hops_handles_ties_dead_nodes_and_destination():
    W = np.array(
        [
            [0.0, 0.0, 0.0, 0.0],
            [2.0, 0.0, 2.0, 0.0],  # tie between 0 and 2 -> lowest id
            [0.0, 0.0, 0.0, 0.0],  # dead
            [0.0, 1.0, 3.0, 0.0],
        ]
    )
    nxt, ties = readout.next_hops(W, 0)
    assert list(nxt) == [0, 0, -1, 2]
    assert list(ties) == [False, True, False, False]


def test_policy_rows_sum_to_one_or_zero():
    W = np.array([[0.0, 5.0, 0.0, 0.0], [2.0, 0.0, 2.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 1.0, 3.0, 0.0]])
    Pm = readout.policy(W, 0)
    assert np.allclose(Pm.sum(axis=1), [0.0, 1.0, 0.0, 1.0])
    assert Pm[3, 2] == pytest.approx(0.75)


def test_zero_field_does_not_crash_the_readout():
    m = Model("M3S", beta=1.0)
    K = m.coupling(topology.path(3), 0)
    x = np.zeros((3, 2))
    assert np.isnan(readout.cancellation(K, x, True)).all()
    nxt, ties = readout.next_hops(readout.weights(readout.power_inflow(K, x)), 0)
    assert list(nxt) == [0, -1, -1] and not ties.any()


def test_m1_with_equal_scales_forwards_like_m0():
    _, _, _, P0 = solve(Model("M0", kappa=1.0), LOP, 5, 3.5)
    _, _, _, P1 = solve(Model("M1", kappa=1.0, kappa2=1.0), LOP, 5, 3.5)
    assert np.allclose(P1, 2 * P0)
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run pytest tests/test_readout.py -q`
Expected: `ImportError: cannot import name 'readout' from 'phaseflow'`.

- [ ] **Step 3: Implement `readout.py`**

```python
"""Forwarding readout shared by every model.

P[i, j] = 2 x_i^T K_ij x_j is the rate at which neighbour j feeds |x_i|^2.
For Hermitian models it is the probability current into i from j.
"""
from __future__ import annotations

import numpy as np


def power_inflow(K: np.ndarray, x: np.ndarray) -> np.ndarray:
    return 2.0 * np.einsum("ia,ijab,jb->ij", x, K, x)


def weights(P: np.ndarray) -> np.ndarray:
    return np.maximum(P, 0.0)


def cancellation(K: np.ndarray, x: np.ndarray, is_complex: bool) -> np.ndarray:
    """1 - |sum of neighbour contributions| / sum of |contributions|, per node.

    A complex state is one channel (modulus over both components); a real
    state is one channel per component. NaN where nothing is received.
    """
    contrib = np.einsum("ijab,jb->ija", K, x)
    if is_complex:
        total = np.linalg.norm(contrib.sum(axis=1), axis=1)
        parts = np.linalg.norm(contrib, axis=2).sum(axis=1)
    else:
        total = np.abs(contrib.sum(axis=1)).sum(axis=1)
        parts = np.abs(contrib).sum(axis=(1, 2))
    out = np.full(len(x), np.nan)
    np.divide(total, parts, out=out, where=parts > 0)
    return 1.0 - out


def next_hops(W: np.ndarray, dest: int, rtol: float = 1e-9) -> tuple[np.ndarray, np.ndarray]:
    """Deterministic next hop per node (-1 if no positive weight) and a tie flag."""
    n = W.shape[0]
    nxt = np.full(n, -1, dtype=int)
    ties = np.zeros(n, dtype=bool)
    nxt[dest] = dest
    for i in range(n):
        if i == dest:
            continue
        best = W[i].max()
        if best <= 0:
            continue
        top = np.flatnonzero(W[i] >= best * (1 - rtol))
        nxt[i] = top[0]
        ties[i] = len(top) > 1
    return nxt, ties


def policy(W: np.ndarray, dest: int) -> np.ndarray:
    """Stochastic forwarding p(i -> j) proportional to W[i, j]."""
    Pm = W.copy()
    Pm[dest] = 0.0
    total = Pm.sum(axis=1, keepdims=True)
    return np.divide(Pm, total, out=np.zeros_like(Pm), where=total > 0)
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `uv run pytest tests/test_readout.py -q`
Expected: `11 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/phaseflow/readout.py tests/test_readout.py
git commit -m "Add power-inflow readout and convention tests A1-A5"
```

---

### Task 5: Route metrics

**Files:**
- Create: `src/phaseflow/metrics.py`
- Test: `tests/test_metrics.py`

**Interfaces:**
- Consumes: `Graph`, `readout.next_hops`, `readout.policy`.
- Produces: `metrics.shortest(g, dest) -> tuple[ndarray, list[set[int]]]` (distance per node, set of optimal next hops per node); `metrics.follow(nxt, src, dest) -> tuple[str, list[int]]` with status `"delivered" | "loop" | "dead"`; `metrics.route_metrics(g, dest, W, canc) -> dict` with keys `nexthop_acc, delivered_frac, n_loop_nodes, n_dead_nodes, n_ties, n_dark_nodes, mean_stretch, max_stretch, mean_cancellation, max_cancellation, min_margin`; `metrics.stochastic_metrics(g, dest, W) -> dict` with keys `p_delivery_min, exp_stretch_mean, entropy_mean`; `metrics.circulation(P, cycle) -> float`; `metrics.split_share(W, node, upper, lower) -> float`.

- [ ] **Step 1: Write the failing tests**

`tests/test_metrics.py`:

```python
import math

import numpy as np
import pytest

from phaseflow import metrics, topology


def test_shortest_reference_on_lopsided_diamond():
    g = topology.lopsided_diamond()  # upper A-B-T costs 2, lower A-C-E-T costs 3
    dist, opt = metrics.shortest(g, 5)
    assert list(dist) == [3, 2, 1, 2, 1, 0]
    assert opt[1] == {2} and opt[3] == {4} and opt[5] == set()
    cheap = g.with_cost(1, 3, 0.5).with_cost(3, 4, 0.5).with_cost(4, 5, 0.5)
    assert metrics.shortest(cheap, 5)[1][1] == {3}
    assert metrics.shortest(topology.diamond(), 6)[1][1] == {2, 3}


def test_follow_statuses():
    assert metrics.follow(np.array([0, 0, 1]), 2, 0) == ("delivered", [2, 1, 0])
    assert metrics.follow(np.array([0, 2, 1]), 1, 0)[0] == "loop"
    assert metrics.follow(np.array([0, -1, 1]), 2, 0)[0] == "dead"


def test_route_metrics_on_hand_built_weights():
    g = topology.lopsided_diamond()
    W = np.zeros((6, 6))
    W[0, 1] = 1
    W[1, 3], W[1, 2] = 2, 1  # A prefers the three-hop lower arm: wrong
    W[2, 5] = W[3, 4] = W[4, 5] = 1
    m = metrics.route_metrics(g, 5, W, np.zeros(6))
    assert m["nexthop_acc"] == pytest.approx(4 / 5)
    assert m["delivered_frac"] == 1.0
    assert m["max_stretch"] == pytest.approx(3 / 2)
    assert m["mean_stretch"] == pytest.approx((4 / 3 + 3 / 2 + 1 + 1 + 1) / 5)
    assert m["n_loop_nodes"] == 0 and m["n_dead_nodes"] == 0 and m["n_ties"] == 0
    assert m["min_margin"] == pytest.approx(1 / 3)


def test_loop_dead_and_dark_nodes_are_counted():
    g = topology.path(4)
    W = np.zeros((4, 4))
    W[1, 2] = W[2, 1] = 1  # 1 <-> 2 loop; node 3 has no positive weight
    canc = np.array([0.0, 0.2, 1.0, np.nan])
    m = metrics.route_metrics(g, 0, W, canc)
    assert m["n_loop_nodes"] == 2 and m["n_dead_nodes"] == 1
    assert m["delivered_frac"] == 0.0 and math.isnan(m["mean_stretch"])
    assert m["n_dark_nodes"] == 1
    assert m["max_cancellation"] == 1.0 and m["mean_cancellation"] == pytest.approx(0.6)


def test_tie_counts_as_correct_when_the_chosen_hop_is_optimal():
    g = topology.diamond()
    W = np.zeros((7, 7))
    W[0, 1] = 1
    W[1, 2] = W[1, 3] = 1  # exact tie at A; both arms are optimal
    W[2, 4] = W[4, 6] = W[3, 5] = W[5, 6] = 1
    m = metrics.route_metrics(g, 6, W, np.zeros(7))
    assert m["n_ties"] == 1 and m["nexthop_acc"] == 1.0 and m["min_margin"] == 0.0


def test_stochastic_metrics_simple_chain():
    g = topology.path(3)
    W = np.zeros((3, 3))
    W[1, 0] = W[1, 2] = 1  # node 1: to the destination or back to 2, each 1/2
    W[2, 1] = 1
    s = metrics.stochastic_metrics(g, 0, W)
    assert s["p_delivery_min"] == pytest.approx(1.0)
    # E1 = 1 + E2 / 2 and E2 = 1 + E1 give E1 = 3, E2 = 4; stretch = (3/1 + 4/2) / 2
    assert s["exp_stretch_mean"] == pytest.approx(2.5)
    assert s["entropy_mean"] == pytest.approx(math.log(2) / 2)


def test_stochastic_metrics_closed_loop_without_destination():
    g = topology.path(3)
    W = np.zeros((3, 3))
    W[1, 2] = W[2, 1] = 1
    s = metrics.stochastic_metrics(g, 0, W)
    assert s["p_delivery_min"] == 0.0
    assert s["exp_stretch_mean"] == math.inf


def test_circulation_and_split_share():
    P = np.zeros((3, 3))
    for a, b in ((0, 1), (1, 2), (2, 0)):  # unit current 0 -> 1 -> 2 -> 0
        P[b, a], P[a, b] = 1, -1
    assert metrics.circulation(P, (0, 1, 2)) == pytest.approx(1.0)
    assert metrics.circulation(np.abs(P), (0, 1, 2)) == pytest.approx(0.0)
    W = np.zeros((4, 4))
    W[1, 2], W[1, 3] = 3, 1
    assert metrics.split_share(W, 1, 2, 3) == pytest.approx(0.75)
    assert math.isnan(metrics.split_share(W, 0, 2, 3))
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run pytest tests/test_metrics.py -q`
Expected: `ImportError: cannot import name 'metrics' from 'phaseflow'`.

- [ ] **Step 3: Implement `metrics.py`**

```python
"""Route quality metrics. Shortest paths come from networkx as an independent reference."""
from __future__ import annotations

import math

import networkx as nx
import numpy as np

from .graph import Graph
from .readout import next_hops, policy

DARK = 1 - 1e-6


def shortest(g: Graph, dest: int) -> tuple[np.ndarray, list[set[int]]]:
    """Cost-to-destination per node, and the set of optimal next hops per node."""
    lengths = nx.single_source_dijkstra_path_length(g.to_networkx(), dest)
    dist = np.array([lengths[i] for i in range(g.n)])
    opt = [
        {j for j in g.neighbors(i) if math.isclose(g.cost(i, j) + dist[j], dist[i], rel_tol=1e-9)}
        for i in range(g.n)
    ]
    opt[dest] = set()
    return dist, opt


def follow(nxt: np.ndarray, src: int, dest: int) -> tuple[str, list[int]]:
    """Follow deterministic next hops from src. Status: delivered, loop or dead."""
    route, seen = [src], {src}
    while route[-1] != dest:
        j = int(nxt[route[-1]])
        if j < 0:
            return "dead", route
        if j in seen:
            return "loop", route + [j]
        route.append(j)
        seen.add(j)
    return "delivered", route


def route_metrics(g: Graph, dest: int, W: np.ndarray, canc: np.ndarray) -> dict:
    dist, opt = shortest(g, dest)
    nxt, ties = next_hops(W, dest)
    nodes = [i for i in range(g.n) if i != dest]
    statuses, stretches, margins = [], [], []
    for i in nodes:
        status, route = follow(nxt, i, dest)
        statuses.append(status)
        if status == "delivered":
            stretches.append(sum(g.cost(u, v) for u, v in zip(route, route[1:])) / dist[i])
        if len(g.neighbors(i)) >= 2:
            total = W[i].sum()
            top = np.sort(W[i])[::-1]
            margins.append(float((top[0] - top[1]) / total) if total > 0 else 0.0)
    c = canc[nodes]
    seen = c[~np.isnan(c)]
    return {
        "nexthop_acc": float(np.mean([int(nxt[i]) in opt[i] for i in nodes])),
        "delivered_frac": statuses.count("delivered") / len(nodes),
        "n_loop_nodes": statuses.count("loop"),
        "n_dead_nodes": statuses.count("dead"),
        "n_ties": int(ties[nodes].sum()),
        "n_dark_nodes": int((seen > DARK).sum()),
        "mean_stretch": float(np.mean(stretches)) if stretches else math.nan,
        "max_stretch": float(np.max(stretches)) if stretches else math.nan,
        "mean_cancellation": float(seen.mean()) if len(seen) else math.nan,
        "max_cancellation": float(seen.max()) if len(seen) else math.nan,
        "min_margin": min(margins) if margins else math.nan,
    }


def stochastic_metrics(g: Graph, dest: int, W: np.ndarray) -> dict:
    """Delivery probability, expected-cost stretch and entropy of p(i -> j) ∝ W[i, j]."""
    dist, _ = shortest(g, dest)
    Pm = policy(W, dest)
    nodes = [i for i in range(g.n) if i != dest]
    Q, r = Pm[np.ix_(nodes, nodes)], Pm[nodes, dest]
    eye = np.eye(len(nodes))
    if np.abs(np.linalg.eigvals(Q)).max() < 1 - 1e-12:
        h = np.linalg.solve(eye - Q, r)
    else:  # a closed class that excludes the destination: take the minimal solution
        h = np.zeros(len(nodes))
        for _ in range(20000):
            h = r + Q @ h
    if h.min() > 1 - 1e-9:
        step_cost = np.array([sum(Pm[i, j] * g.cost(i, j) for j in g.neighbors(i)) for i in nodes])
        exp_stretch = float(np.mean(np.linalg.solve(eye - Q, step_cost) / dist[nodes]))
    else:
        exp_stretch = math.inf
    logs = np.log(Pm, out=np.zeros_like(Pm), where=Pm > 0)
    entropy = -(Pm * logs).sum(axis=1)
    return {
        "p_delivery_min": float(h.min()),
        "exp_stretch_mean": exp_stretch,
        "entropy_mean": float(entropy[nodes].mean()),
    }


def circulation(P: np.ndarray, cycle: tuple[int, ...]) -> float:
    """Mean antisymmetric flow along the cycle in the order given."""
    m = len(cycle)
    flows = [(P[cycle[(t + 1) % m], cycle[t]] - P[cycle[t], cycle[(t + 1) % m]]) / 2 for t in range(m)]
    return float(np.mean(flows))


def split_share(W: np.ndarray, node: int, upper: int, lower: int) -> float:
    """Share of node's two-arm forwarding weight placed on the upper arm."""
    total = W[node, upper] + W[node, lower]
    return float(W[node, upper] / total) if total > 0 else math.nan
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `uv run pytest tests/test_metrics.py -q`
Expected: `8 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/phaseflow/metrics.py tests/test_metrics.py
git commit -m "Add route, stretch, loop and stochastic-policy metrics"
```

---

### Task 6: Perturbations, scenarios and the predeclared grid

**Files:**
- Create: `src/phaseflow/perturbations.py`, `src/phaseflow/scenarios.py`
- Test: `tests/test_scenarios.py`

**Interfaces:**
- Consumes: `Graph`, topologies, `Model`, `metrics.shortest` (tests only).
- Produces: `CostChange(round: int, edge: tuple[int, int], cost: float)`; `phases(g, changes, total_rounds) -> list[tuple[int, int, Graph]]` as `(start, end, graph in force)`.
- Produces: `Scenario(name, family, graph, dest, cycle=(), split=())`; `Curve(name, crossover, points)` where `points` is a tuple of `(swept_cost, Scenario)`; `DynamicCase(scenario, changes, total_rounds)`; `Config(scenarios, models, gammas, curves, curve_models, curve_gammas, dynamic_cases, dynamic_models, dynamic_gammas)`; `static_scenarios()`, `model_grid()`, `full_config()`, `quick_config()`.
- `family` is one of `"tree"`, `"diamond"`, `"lopsided"`, `"ring"`. `split` is `(node, upper_neighbour, lower_neighbour)`.

- [ ] **Step 1: Write the failing tests**

`tests/test_scenarios.py`:

```python
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
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run pytest tests/test_scenarios.py -q`
Expected: `ImportError: cannot import name 'scenarios' from 'phaseflow'`.

- [ ] **Step 3: Implement `perturbations.py` and `scenarios.py`**

`src/phaseflow/perturbations.py`:

```python
"""Deterministic link-cost change schedule."""
from __future__ import annotations

from dataclasses import dataclass

from .graph import Graph


@dataclass(frozen=True)
class CostChange:
    round: int
    edge: tuple[int, int]
    cost: float


def phases(g: Graph, changes, total_rounds: int) -> list[tuple[int, int, Graph]]:
    """Split [0, total_rounds) at each change. Returns (start, end, graph in force)."""
    out, start, current = [], 0, g
    for change in sorted(changes, key=lambda c: c.round):
        if not start < change.round < total_rounds:
            raise ValueError(f"change round {change.round} must lie strictly inside ({start}, {total_rounds})")
        out.append((start, change.round, current))
        current = current.with_cost(*change.edge, change.cost)
        start = change.round
    out.append((start, total_rounds, current))
    return out
```

`src/phaseflow/scenarios.py`:

```python
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
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `uv run pytest tests/test_scenarios.py -q`
Expected: `5 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/phaseflow/perturbations.py src/phaseflow/scenarios.py tests/test_scenarios.py
git commit -m "Add cost-change schedule, scenarios and the predeclared grid"
```

---

### Task 7: Static sweep, control checks and transfer table

**Files:**
- Create: `src/phaseflow/experiment0.py`
- Test: `tests/test_experiment0.py`

**Interfaces:**
- Consumes: everything above.
- Produces: `experiment0.evaluate(scn, model, gamma) -> tuple[dict, ndarray | None]`; `experiment0.StaticResult(metrics: DataFrame, readouts: DataFrame, inflows: dict)`; `experiment0.run_static(scns, models, gammas) -> StaticResult`; `experiment0.control_checks(res) -> dict` with keys `tree_gauge_max_rel_diff`, `diamond_m3h_split_max_dev`, `tree_cells_all_correct`; `experiment0.transfer_table(metrics_df) -> DataFrame` with columns `kind, variant, param, gamma, n_scenarios, stable_all, min_acc, mean_acc, min_delivered, max_stretch, mean_cancellation`; `experiment0.summarize_groups(transfer) -> dict`; `experiment0.run_discrimination(cfg) -> DataFrame`.
- Row id columns everywhere: `scenario, family, kind, variant, param, gamma`. `inflows` is keyed by `(scenario, variant, param, gamma)`.

- [ ] **Step 1: Write the failing tests**

`tests/test_experiment0.py`:

```python
import math

import numpy as np
import pandas as pd

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
    assert row["split_upper"] == 0.5


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
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run pytest tests/test_experiment0.py -q`
Expected: `ImportError: cannot import name 'experiment0' from 'phaseflow'`.

- [ ] **Step 3: Implement the static half of `experiment0.py`**

```python
"""Experiment 0 smoke test: static sweep, discrimination curves, dynamic run."""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import engine, metrics, oracle, readout
from .models import Model
from .scenarios import Config, Scenario

STABLE_TOL = 1e-9
ID_COLS = ["scenario", "family", "kind", "variant", "param", "gamma"]
METRIC_COLS = [
    "nexthop_acc", "delivered_frac", "n_loop_nodes", "n_dead_nodes", "n_ties", "n_dark_nodes",
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
        **dict.fromkeys(METRIC_COLS, math.nan),
    }
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
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `uv run pytest tests/test_experiment0.py -q`
Expected: `8 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/phaseflow/experiment0.py tests/test_experiment0.py
git commit -m "Add static sweep, control checks and transfer table"
```

---

### Task 8: Dynamic perturbation run

**Files:**
- Modify: `src/phaseflow/experiment0.py` (append)
- Test: `tests/test_dynamic.py`

**Interfaces:**
- Consumes: `engine.step`, `engine.step_sizes`, `perturbations.phases`, `CostChange`, `Scenario`, `_ids`, `STABLE_TOL`.
- Produces: `experiment0.settle_round(ok: ndarray[bool]) -> float`; `experiment0.DynamicResult(phases, series, states, currents)` (four DataFrames); `experiment0.run_dynamic(scn, changes, model, gamma, total_rounds) -> DynamicResult`; `experiment0.run_all_dynamic(cfg) -> DynamicResult`.
- `phases` columns: id columns plus `phase, start, rounds, stable, diverged, settle_rounds, route_settle_rounds, loop_rounds, final_nexthop_acc, alpha_mean`. `series` columns: id columns plus `phase, round, residual, split_upper, n_correct`. `states`: id columns plus `round, node, x0, x1, magnitude`. `currents`: id columns plus `round, i, j, inflow`.
- A round is logged when it is one of the first 60 of its phase, a multiple of 50 within its phase, or the last round of its phase.

- [ ] **Step 1: Write the failing tests**

`tests/test_dynamic.py`:

```python
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
    assert split[4499] == split[1499] or abs(split[4499] - split[1499]) < 1e-6  # no hysteresis


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
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run pytest tests/test_dynamic.py -q`
Expected: `AttributeError: module 'phaseflow.experiment0' has no attribute 'settle_round'`.

- [ ] **Step 3: Append the dynamic runner to `experiment0.py`**

Add `from .perturbations import phases` to the imports, then append:

```python
SETTLE_TOL = 1e-6
DIVERGE_AT = 1e12
LOG_DENSE = 60
LOG_EVERY = 50


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
        done = hops and not diverged
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
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `uv run pytest tests/test_dynamic.py -q`
Expected: `5 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/phaseflow/experiment0.py tests/test_dynamic.py
git commit -m "Add dynamic cost-perturbation run with settle and loop metrics"
```

---

### Task 9: Writers, plots and the command-line entry point

**Files:**
- Modify: `src/phaseflow/experiment0.py` (append)
- Create: `src/phaseflow/plots.py`, `README.md`
- Test: `tests/test_main.py`

**Interfaces:**
- Consumes: `run_static`, `run_discrimination`, `run_all_dynamic`, `transfer_table`, `summarize_groups`, `control_checks`, `full_config`, `quick_config`.
- Produces: `experiment0.main(out_dir="results/experiment0-smoke", quick=False, make_plots=True) -> dict` and the CLI `python -m phaseflow.experiment0 [--out DIR] [--quick] [--no-plots]`.
- Produces: `plots.make_all(out: Path) -> list[Path]`, reading only the CSVs in `out`.
- Files written: `config.json, summary.json, route_metrics.csv, static_readouts.csv, transfer.csv, discrimination.csv, dynamic_metrics.csv, dynamic_series.csv, states.csv, currents.csv, plots/*.png`.

Plot rules (from the dataviz guidance): one y-axis per panel; model identity uses one fixed colour per variant across every figure plus a line style; ordered parameters use a single-hue ramp; grid lines are solid hairlines; text is never drawn in a series colour; a legend is present whenever there are two or more series.

- [ ] **Step 1: Write the failing test**

`tests/test_main.py`:

```python
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
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `uv run pytest tests/test_main.py -q`
Expected: `AttributeError: module 'phaseflow.experiment0' has no attribute 'main'`.

- [ ] **Step 3: Append the writers and CLI to `experiment0.py`**

Add these imports at the top of the file: `argparse`, `dataclasses`, `json`, `platform`, `subprocess`, `from importlib.metadata import version`, `from pathlib import Path`, and `from .scenarios import full_config, quick_config`. Then append:

```python
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
```

- [ ] **Step 4: Implement `plots.py`**

```python
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
SEQUENTIAL = LinearSegmentedColormap.from_list("phaseflow_blue", ["#cde2fb", "#3987e5", "#0d366b"])
SEQUENTIAL.set_bad(MASK)
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
        fig, axes = plt.subplots(rows, cols, figsize=(3.3 * cols, 2.7 * rows), sharex=True, sharey=True, squeeze=False)
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
        fig.legend(handles=handles, loc="lower center", ncol=len(handles), bbox_to_anchor=(0.5, -0.02))
        fig.supxlabel("swept cost on the lower arm (arms cost the same at the vertical line)", y=0.06, color=MUTED)
        fig.supylabel("share of A's forwarding weight on the upper arm", color=MUTED)
        fig.suptitle(f"{curve}: a correct model is below 0.5 left of the line and above 0.5 right of it")
        fig.tight_layout(rect=(0, 0.08, 1, 0.96))
        made.append(_save(fig, plots / f"discrimination_{curve}.png"))
    return made


def plot_accuracy_heatmaps(transfer: pd.DataFrame, plots: Path) -> list[Path]:
    """Mean next-hop accuracy over the cyclic scenarios, per parameter cell."""
    panels = {}
    arbitrary = transfer[(transfer["kind"] == "M3H") & (transfer["variant"] != ORACLE)]
    for variant in ("M0", "M1", "M2", "M3H[id]", ORACLE, "M3S[theta0=-pi/2]", "M3S[theta0=0]"):
        sub = transfer[transfer["variant"] == variant]
        if len(sub):
            panels[variant] = sub.pivot(index="gamma", columns="param", values="mean_acc")
    if len(arbitrary):
        worst = arbitrary.groupby(["gamma", "param"])["mean_acc"].min().unstack()
        panels["M3H worst arbitrary orientation"] = worst
    if not panels:
        return []
    cols = min(4, len(panels))
    rows = math.ceil(len(panels) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(3.2 * cols, 2.8 * rows), squeeze=False)
    image = None
    for ax, (title, grid) in zip(axes.flat, panels.items()):
        image = ax.imshow(np.ma.masked_invalid(grid.to_numpy(dtype=float)), cmap=SEQUENTIAL, vmin=0, vmax=1, origin="lower", aspect="auto")
        ax.set_xticks(range(len(grid.columns)), [f"{c:g}" for c in grid.columns])
        ax.set_yticks(range(len(grid.index)), [f"{g:g}" for g in grid.index])
        ax.grid(False)
        ax.set_title(title)
        ax.set_xlabel("β or κ")
        ax.set_ylabel("γ")
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    fig.suptitle("Mean next-hop accuracy over cyclic scenarios (grey: unstable on at least one)")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.colorbar(image, ax=axes, shrink=0.8, label="share of nodes forwarding along a shortest path")
    return [_save(fig, plots / "accuracy_heatmaps.png")]


def plot_orientation_spread(metrics: pd.DataFrame, plots: Path) -> list[Path]:
    """M3-H next-hop accuracy per cyclic scenario, across orientations and parameters."""
    df = metrics[(metrics["kind"] == "M3H") & (metrics["family"] != "tree") & metrics["stable"]]
    if df.empty:
        return []
    names = list(dict.fromkeys(df["scenario"]))
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(1.1 * len(names) + 2, 3.6))
    groups = (("arbitrary orientations", df[df["variant"] != ORACLE], "#2a78d6"), ("oracle_hops (diagnostic)", df[df["variant"] == ORACLE], "#eb6834"))
    for offset, (label, sub, colour) in zip((-0.18, 0.18), groups):
        if sub.empty:
            continue
        x = sub["scenario"].map(names.index) + offset + rng.uniform(-0.12, 0.12, len(sub))
        ax.scatter(x, sub["nexthop_acc"], s=28, color=colour, alpha=0.55, edgecolors=SURFACE, linewidths=0.8, label=label)
    ax.set_xticks(range(len(names)), names, rotation=30, ha="right")
    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel("next-hop accuracy")
    ax.set_title("M3-H: one dot per orientation × β × γ")
    ax.legend(loc="lower left")
    return [_save(fig, plots / "m3h_orientation_spread.png")]


def plot_cancellation(metrics: pd.DataFrame, plots: Path) -> list[Path]:
    """Does more destructive interference go with better or worse forwarding?"""
    df = metrics[(metrics["family"] != "tree") & metrics["stable"]]
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    drawn = 0
    for variant in ("M3H[id]", "M3S[theta0=-pi/2]", "M3S[theta0=0]"):
        sub = df[df["variant"] == variant]
        if sub.empty:
            continue
        ax.scatter(sub["mean_cancellation"], sub["nexthop_acc"], s=28, color=SERIES[variant][0], alpha=0.6, edgecolors=SURFACE, linewidths=0.8, label=variant)
        drawn += 1
    if not drawn:
        plt.close(fig)
        return []
    ax.set_xlabel("mean cancellation at a node (0: contributions add, 1: they cancel)")
    ax.set_ylabel("next-hop accuracy")
    ax.set_ylim(-0.03, 1.03)
    ax.set_title("Interference against forwarding quality, one dot per cyclic cell")
    ax.legend(loc="lower left")
    return [_save(fig, plots / "cancellation_vs_accuracy.png")]


def plot_dynamic(series: pd.DataFrame, phases: pd.DataFrame, plots: Path) -> list[Path]:
    """Distance to the fixed point and A's split over time, through the cost changes."""
    made = []
    for (scenario, gamma), df in series.groupby(["scenario", "gamma"]):
        fig, (top, bottom) = plt.subplots(2, 1, figsize=(7.5, 5.4), sharex=True)
        for variant, (colour, style) in SERIES.items():
            sub = df[df["variant"] == variant].sort_values("round")
            if sub.empty:
                continue
            top.semilogy(sub["round"], sub["residual"].clip(lower=1e-16), color=colour, linestyle=style, label=variant)
            bottom.plot(sub["round"], sub["split_upper"], color=colour, linestyle=style)
        starts = phases[(phases["scenario"] == scenario) & (phases["start"] > 0)]["start"].unique()
        for ax in (top, bottom):
            for start in starts:
                ax.axvline(start, color=MUTED, linewidth=0.8)
        top.set_ylabel("relative distance to the oracle fixed point")
        bottom.set_ylabel("share of A's weight on the upper arm")
        bottom.set_xlabel("round (vertical lines: cost raised, then restored)")
        bottom.set_ylim(-0.03, 1.03)
        fig.legend(loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.955))
        fig.suptitle(f"{scenario}, γ = {gamma:g}: unstable models have no distance curve")
        fig.tight_layout(rect=(0, 0, 1, 0.9))
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
```

- [ ] **Step 5: Write `README.md`**

````markdown
# PhaseFlow

Exploratory study of whether interference in a locally propagated complex
field is a useful routing signal on classical networks. Start with
`docs/phaseflow_kickoff.md` and `docs/phaseflow_research_charter.md`.

## Reproduce Experiment 0 (smoke test)

```bash
uv sync
uv run pytest -q
uv run python -m phaseflow.experiment0
```

Outputs go to `results/experiment0-smoke/`. The run is deterministic.
Design: `docs/superpowers/specs/2026-10-02-experiment0-smoke-design.md`.
Findings and the open decision: `docs/experiment-0.md`.
````

- [ ] **Step 6: Run the test and confirm it passes**

Run: `uv run pytest tests/test_main.py -q`
Expected: `1 passed`.

- [ ] **Step 7: Run the whole suite**

Run: `uv run pytest -q`
Expected: every test passes, no warnings about invalid values or overflow.

- [ ] **Step 8: Commit**

```bash
git add src/phaseflow/experiment0.py src/phaseflow/plots.py tests/test_main.py README.md
git commit -m "Add output writers, plots and the experiment entry point"
```

---

### Task 10: Run the experiment and record the evidence

**Files:**
- Create: `results/experiment0-smoke/*`, `docs/experiment-0.md`, `docs/literature.md`

No new code. This task produces the evidence and the decision document required by the charter §5 and spec §14.

- [ ] **Step 1: Run the full experiment**

Run: `uv run python -m phaseflow.experiment0`
Expected: prints the control checks; `tree_gauge_max_rel_diff` below `1e-10`, `diamond_m3h_split_max_dev` below `1e-9`.

- [ ] **Step 2: Look at every plot**

Open each PNG in `results/experiment0-smoke/plots/`. Check for overlapping labels, clipped legends and empty panels. Fix `plots.py` if any figure is unreadable, re-run, and re-check.

- [ ] **Step 3: Check the predictions against the data**

From `route_metrics.csv`, `transfer.csv`, `discrimination.csv` and `dynamic_metrics.csv`, mark each of Pr1–Pr4 in the spec confirmed or refuted, with the numbers that decide it. Record counterexamples by scenario name.

- [ ] **Step 4: Write `docs/literature.md`**

Save the prior-art note. Every listed reference carries the URL that was retrieved to confirm it; anything recalled but not retrieved goes under "Unverified leads".

- [ ] **Step 5: Write `docs/experiment-0.md`**

Sections, in this order:

1. What was run (one paragraph, pointing at the spec and `config.json`).
2. Control checks (A1–A6 results).
3. M3-H findings.
4. M3-S findings.
5. Real-valued controls (M0, M1) and the phase-disabled model (M2).
6. Predictions Pr1–Pr4: confirmed or refuted.
7. The charter's five questions, answered separately for M3-H and M3-S.
8. Limitations of this smoke test.
9. Recommendation. The decision is left to the project owner.

Every number in the document must come from a file in `results/experiment0-smoke/`.

- [ ] **Step 6: Commit**

```bash
git add results/experiment0-smoke docs/experiment-0.md docs/literature.md
git commit -m "Record Experiment 0 smoke results, literature note and decision draft"
```
