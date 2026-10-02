# Experiment 0 smoke test — design

**Date:** 2026-10-02
**Status:** Exploratory. Nothing here is confirmatory; see `docs/phaseflow_research_charter.md` §2.
**Inputs:** `docs/phaseflow_kickoff.md` (§19 first milestone), `docs/phaseflow_research_charter.md`.

## 1. Purpose

Find out, on a tree, a diamond and a ring, whether either of two phase-field
formulations is mathematically coherent as a routing primitive and whether
either shows an effect that depends on interference rather than on extra
state, damping, or a weak baseline.

The two formulations are separate exploratory hypotheses and are reported
separately. Results from one are never combined with results from the other
into a single claim.

- **M3-H** — Hermitian / chiral field. Strong mathematical controls.
- **M3-S** — complex-symmetric / accumulated path-phase field. Closest to the
  original "routes interfere by cost" intuition.

The output of this milestone is evidence plus a draft `docs/experiment-0.md`
answering the charter's five questions. The decision to select a branch,
revise the question once, or stop belongs to the project owner.

## 2. Scope

In scope: tree, diamond, ring topologies; models M0, M1, M2, M3-H, M3-S; one
unified local engine; exact fixed-point oracle; forwarding readout and route
extraction; one deterministic cost perturbation (raise, then restore); CSV,
JSON and plot outputs; unit tests for conventions.

Out of scope: random meshes, packet simulation, UI, Bellman-Ford / ECMP /
Q-routing / AntNet baselines, multi-frequency fields, asymmetric link costs,
confirmatory statistics.

## 3. Unified formulation

Every model is the same linear local system. Only the per-edge coupling block
differs. This is deliberate: it removes per-model readout or engine choices
as a way to favour one model.

Each node `i` holds a real state `x_i ∈ R^k` for one destination `d`
(`k = 1` or `2`). For neighbours `j ∈ N(i)` there is a real `k×k` block
`K_ij`, meaning "what `i` receives from `j`":

```text
dx_i/dt = -γ x_i + Σ_{j∈N(i)} K_ij x_j + b_i          b_i = 0 for i ≠ d
```

A complex number `p + iq` is the block `[[p, -q], [q, p]]`, so complex models
are literally two reals per node. The kickoff form `dψ/dt = -(γI + iH)ψ + b`
is this system with `K = -iH`.

**Local iteration** (synchronous rounds; one round = every node broadcasts
its `k` reals to each neighbour, then updates):

```text
x_i ← x_i + α_i ( -γ x_i + Σ_j K_ij x_j + b_i )
α_i = γ / (γ² + s_i²),      s_i = Σ_{j∈N(i)} ‖K_ij‖₂
```

`α_i` uses only quantities node `i` already has, and is recomputed by the
node when one of its link costs changes. For Hermitian models the iteration
converges for every `γ > 0` on every graph (each eigenvalue `μ` of the scaled
operator satisfies `|1 - μ|² ≤ 1 - Re μ`); for non-negative real couplings it
is a regular splitting, so it converges whenever the continuous system is
stable. M3-S has no such guarantee; the oracle reports the iteration's
spectral radius per cell. The kickoff's §5 update omitted the factor `-i`;
this form restores it.

**Oracle** (validation only, never visible to nodes): exact fixed point
`x* = (γI - K)⁻¹ b`, spectral abscissa `max Re λ(K) - γ` (continuous
stability), and spectral radius of the iteration matrix (discrete
convergence).

## 4. Models

`a = 1` throughout. `c_ij = c_ji > 0` is the link cost. `R(θ)` is the 2×2
rotation by `θ`.

| ID | State | Block `K_ij` | Cost enters via | Source `b_d` |
|---|---|---|---|---|
| M0 | 1 real | `a·e^{-κ c_ij}` | attenuation | `1` |
| M1 | 2 reals | `a·diag(e^{-κ₁ c_ij}, e^{-κ₂ c_ij})` | attenuation, two scales | `(1, 1)` |
| M2 | 2 reals (complex) | `a·R(-π/2)` | not at all | `(1, 0)` |
| M3-H | 2 reals (complex) | `a·R(-π/2 + σ_ij β c_ij)`, `σ_ji = -σ_ij` | loop flux | `(1, 0)` |
| M3-S | 2 reals (complex) | `a·R(θ₀ + β c_ij)`, same both directions | walk phase | `(1, 0)` |

- **M0** is cost-weighted diffusion. Its fixed point is a sum over walks
  weighted by `e^{-κ·cost}`, a soft-min of path cost. It is the real-valued
  image of M3-S under `iβ → -κ`.
- **M1** is the matched two-real control: same state size, same payload, same
  per-edge information as the complex models, same engine and readout, but
  the block is non-negative and diagonal, so contributions can only add.
  `κ₁ = κ₂` reduces it to M0.
- **M2** is M3 with `β = 0`. It keeps the coherent `-i` dynamics (walks of
  different hop count still interfere) but has no cost information.
- **M3-H**: `H_ij = a·e^{iσ_ij β c_ij}`, `H_ji = conj(H_ij)`, `K = -iH`.
  Stable for every `γ > 0`. Needs an edge orientation `σ`; see §4.1.
- **M3-S**: `S_ij = S_ji = a·e^{iβ c_ij}`, `K = e^{iθ₀}S`. `θ₀ = -π/2`
  is the kickoff form `K = -iS`; `θ₀ = 0` removes the per-hop phase so walk
  phase is exactly `β × cost`. `θ₀` is part of the cost-to-phase rule
  `φ = βc + θ₀` and both values are run. Stable only when
  `γ > max Re λ(K)`; the Gershgorin bound `γ > a·deg_max` is sufficient and
  locally computable.

### 4.1 M3-H orientation

A symmetric cost needs a sign convention for which direction gets `+βc`.
Conventions run, none tuned per graph:

- `id`: `σ_ij = +1` when `i < j`.
- `reverse`: the negation of `id`.
- `random:s` for seeds `s = 0..4`: an independent fair coin per edge.
- `oracle_hops` — **diagnostic only**: `+βc` in the direction away from the
  destination by hop count. This uses global destination knowledge, so it is
  excluded from every M3-H verdict. It answers one question: would M3-H work
  if the orientation problem were solved?

### 4.2 Information and budget

| | M0 | M1 | M2 | M3-H | M3-S |
|---|---|---|---|---|---|
| Reals stored per node per destination | 1 | 2 | 2 | 2 | 2 |
| Reals sent per neighbour per round | 1 | 2 | 2 | 2 | 2 |
| Uses incident link costs | yes | yes | no | yes | yes |
| Global protocol constants | `γ, κ` | `γ, κ₁, κ₂` | `γ` | `γ, β` | `γ, β, θ₀` |
| Extra shared knowledge | none | none | none | per-edge sign `σ` | none |

All models get the same number of rounds and the same step rule. M0 sends
half the payload of the others; this is reported, not corrected.

## 5. Readout and route

One readout for every model, the **local power inflow**:

```text
P_i←j = 2 x_iᵀ K_ij x_j          (= 2 Re(ψ_i* K_ij ψ_j) for complex models)
```

It is the rate at which neighbour `j` feeds `|x_i|²`, from the exact identity

```text
d|x_i|²/dt = -2γ|x_i|² + Σ_j P_i←j + 2 x_iᵀ b_i
```

Node `i` computes it from its own state, its own blocks and the states it
received. At a fixed point `Σ_j P_i←j = 2γ|x_i|²` for `i ≠ d`.

What it means per model:

- **M2, M3-H**: `P_i←j = 2 Im(H_ij ψ_i* ψ_j)`, the kickoff's `J_ij`. It is
  antisymmetric, so it is a genuine edge current. Positive `J_ij` is flow
  **from `j` into `i`**. The field is injected at the destination and
  dissipated everywhere, so it is not globally conserved and it flows
  outward from `d`.
- **M3-S**: `P` is not antisymmetric and is not a current. It measures how
  far `j`'s contribution reinforces (positive) or cancels (negative) the
  field at `i`. Both `P_i←j` and `P_j←i` may be positive.
- **M0, M1**: every `P_i←j ≥ 0`; for fixed `i` it is proportional to `j`'s
  contribution `K_ij x_j`.

**Forwarding.** A packet at `i` bound for `d` moves upstream, toward the
neighbour that feeds `i` most:

```text
weights        w_i→j = max(P_i←j, 0)
deterministic  next(i) = argmax_j w_i→j     (ties: lowest id, tie recorded)
stochastic     p(i→j) = w_i→j / Σ_j w_i→j
```

**Route.** From `s`, follow `next(·)` until reaching `d` (delivered),
revisiting a node (loop), or reaching a node with no positive weight (dead).
This is the definition of "route" for all five models.

**Interference measure.** `cancellation_i = 1 - |Σ_j K_ij x_j| / Σ_j |K_ij x_j|`.
It is 0 when contributions add in phase and 1 at perfect cancellation. It is
computed per channel: a complex state is one channel (modulus), a real state
is one channel per component. It is therefore identically 0 for M0 and M1. A
node with cancellation above `1 - 10⁻⁶` is *dark*.

## 6. Pre-experiment analysis

Derived before the run and checked in a throwaway numeric probe. Each becomes
a unit test.

- **A1. Sign.** Two nodes, destination 0: `P_1←0 = 2γ|x_1|² > 0` and, for
  Hermitian models, `P_0←1 = -P_1←0`. Node 1 forwards to node 0.
- **A2. Balance.** The identity in §5 holds at the fixed point for every
  model, to machine precision.
- **A3. Tree gauge invariance (M3-H).** On a tree, `|ψ_i|` and `P` equal
  those of M2 for every `β` and every orientation. Raw `ψ` differs. Compare
  gauge-invariant quantities only.
- **A4. Bipartite flux-sign blindness (M3-H).** If the graph is bipartite,
  reversing every phase (`φ → -φ`) leaves `|ψ_i|` and all currents
  unchanged. Proof: with `D = diag(±1)` by sublattice, `D H* D = -H*`, so
  `Dψ*` solves the reversed system up to a global sign, and `D_i D_j = -1`
  on every edge restores the sign of `J`. Odd cycles are not covered and do
  differ.
- **A5. Corollary: the kickoff diamond cannot test M3-H.** Its only cycle
  has six nodes and a mirror symmetry exchanging the arms. By gauge
  invariance only the loop flux `Φ` matters, not how cost is distributed
  along the arms. Mirror symmetry then gives `split(Φ) = 1 - split(-Φ)`, and
  A4 gives `split(Φ) = split(-Φ)`. So node
  A splits exactly 50/50 for every cost assignment, `β`, `γ` and
  orientation. The suite therefore adds odd-cycle graphs (§7).
- **A6. Stability.** M3-H is stable for all `γ > 0`. M3-S and the real
  models are stable only above a graph-dependent threshold.

Predictions recorded before the run, to be confirmed or refuted:

- **Pr1.** On trees every stable model forwards every node correctly.
- **Pr2.** M0 forwards correctly at large `κ` and is biased toward fewer
  hops at small `κ`.
- **Pr3.** M3-S on equal-hop arms discriminates cost only weakly, through
  multi-pass walks, and the direction of its preference depends on `θ₀` and
  `β`. (The probe saw this at one parameter cell.)
- **Pr4.** M3-H on odd cycles is orientation-dependent, and no arbitrary
  orientation is correct across cost assignments.

## 7. Topologies and scenarios

One destination per scenario; every other node is evaluated.

| Family | Scenario | Costs | Tests |
|---|---|---|---|
| Tree | `path6` | unit | sign, hop decay |
| Tree | `path6_rand` | U[0.5, 2], seed 1 | gauge control |
| Tree | `btree7_root`, `btree7_leaf` | U[0.5, 2], seed 2 | gauge control, branching |
| Diamond | `diamond_equal` | unit | tie |
| Diamond | `diamond_lower+0.2`, `diamond_lower+1.0` | C–E raised | cost discrimination |
| Diamond | `diamond_upper+1.0` | B–D raised | mirror check |
| Lopsided diamond (5-cycle) | `lop_unit` | unit: 2 hops vs 3 | hops and cost agree |
| Lopsided diamond | `lop_lower_cheap` | lower edges 0.5 | hops and cost disagree |
| Lopsided diamond | `lop_near_tie` | lower edges 0.7 | near tie |
| Ring | `ring8_unit`, `ring7_unit` | unit | even vs odd cycle |
| Ring | `ring8_heavy`, `ring7_heavy` | edge 1–2 cost 4 | detour |

The diamond is the kickoff's (S–A, arms A–B–D–T and A–C–E–T). The lopsided
diamond has arms A–B–T and A–C–E–T. Both are additions required by A5, as is
`ring7`.

**Discrimination curves.** Split fraction at A (share of A's forwarding
weight on the upper arm) against a swept cost: diamond C–E cost from 0.5 to
2.5; lopsided lower-arm edge cost from 0.3 to 1.5 (true crossover at 2/3).

## 8. Parameter grid

Declared here, before the run. The whole grid is reported; no cell is
selected and presented alone.

- `γ ∈ {0.1, 0.3, 1.0, 2.0, 3.5}`
- M0: `κ ∈ {0.25, 0.5, 1, 2, 4}`
- M1: `(κ₁, κ₂) = (κ, 2κ)` for the same five `κ`
- M2: no further parameter
- M3-H: `β ∈ {0.25, 0.5, 1, 2, 4}` × the eight orientations of §4.1
- M3-S: `β ∈ {0.25, 0.5, 1, 2, 4}` × `θ₀ ∈ {-π/2, 0}`

Unstable cells are recorded as unstable and their readouts are not
evaluated. Discrimination curves use `γ ∈ {0.3, 3.5}` and
`β, κ ∈ {0.5, 1, 2}`.

**Transfer summary** (descriptive, same procedure for every model): for each
parameter cell, whether it is stable on all cyclic scenarios, and its
minimum and mean next-hop accuracy and delivery over them. Report the full
distribution over cells and the best single cell.

## 9. Dynamic perturbation

On `diamond` and `lopsided`, from `x = 0`:

```text
phase 1   base costs                                   2500 rounds
phase 2   diamond: B–D 1 → 3;  lopsided: B–T 1 → 3     2500 rounds
phase 3   restore                                      2500 rounds
```

Run at `γ ∈ {0.3, 1.0, 3.5}`, `β = κ = 1`, M1 at `(1, 2)`, M3-H orientation
`id`, M3-S at both `θ₀`. A run whose state norm exceeds `10¹²` is stopped and
recorded as diverged.

Per phase: rounds until the relative distance to the oracle fixed point
stays below `10⁻⁶`; rounds until the next-hop map stops changing; rounds
spent with a forwarding loop; final next-hop accuracy.

Settle times are comparable across models only as far as the shared step
rule is fair to each. They are reported in rounds and in intrinsic time
`Σ α`, and are treated as descriptive.

## 10. Metrics

Per scenario × model × parameter cell, at the oracle fixed point:

- stability: spectral abscissa, discrete spectral radius, both flags;
- `nexthop_acc`: share of nodes whose next hop lies on a shortest path;
- `delivered_frac`, `n_loop_nodes`, `n_dead_nodes`, `n_dark_nodes`, `n_ties`;
- `mean_stretch`, `max_stretch` over delivered nodes;
- stochastic policy: minimum delivery probability, mean expected-cost
  stretch, mean forwarding entropy;
- `mean_cancellation`, `max_cancellation`;
- `min_margin`: smallest gap between the top two normalised weights;
- `split_upper` at node A (diamond families);
- `circulation`: mean antisymmetric part of `P` around the cycle.

Shortest paths come from `networkx` as an independent reference.

## 11. Tests that gate the sweep

1. Two-node sign convention and forwarding direction, all models (A1).
2. Balance identity, all models (A2).
3. `P` antisymmetric for M2 and M3-H; not for M3-S.
4. Tree gauge invariance on `|ψ|` and `P`, with raw `ψ` shown to differ (A3).
5. Bipartite flux-sign blindness, and its failure on odd cycles (A4).
6. Diamond 50/50 split for M3-H (A5).
7. Orientation: `σ` antisymmetric, `H` Hermitian; `S` symmetric and
   non-Hermitian.
8. Engine: converges to the oracle fixed point; a node's update is
   unaffected by non-neighbour state; real-block form equals a direct
   complex solve.
9. Stability: M3-H stable at small `γ`; M3-S diverges below its threshold
   and converges above it, as the oracle predicts.
10. M1 with `κ₁ = κ₂` gives the same forwarding as M0.
11. Metrics: shortest-path next hops, loop and dead-node detection, stretch,
    absorbing-chain expected cost, on hand-checked graphs.
12. Perturbation: after a cost change the engine re-converges to the new
    oracle fixed point.
13. Determinism: two runs with the same configuration give identical output.

## 12. Outputs

```text
results/experiment0-smoke/
  config.json            grid, seeds, package versions, git commit
  summary.json           transfer tables, control-check residuals, counts
  route_metrics.csv      static sweep, one row per scenario × model × cell
  static_readouts.csv    per directed edge P and weights for the static sweep
  transfer.csv           per parameter cell, aggregated over cyclic scenarios
  discrimination.csv     split-vs-cost curves
  dynamic_metrics.csv    per-phase settle and route metrics
  dynamic_series.csv     dynamic run, residual and split over time (subsampled)
  states.csv             dynamic run, per node state over time (subsampled)
  currents.csv           dynamic run, per directed edge P over time (subsampled)
  plots/
```

Run with `uv run python -m phaseflow.experiment0`. The run is deterministic.
Dynamic time series keep the first 60 rounds of each phase, every 50th round
after that, and the last round of each phase.

Plots: discrimination curves per model; next-hop accuracy heatmaps over the
grid with unstable cells masked; M3-H orientation spread; cancellation
against accuracy; dynamic residual, split and state magnitude over time.

## 13. Code layout

```text
pyproject.toml
src/phaseflow/
  graph.py          undirected weighted graph, cost updates
  topology.py       path, binary tree, diamond, lopsided diamond, ring
  models.py         the five coupling builders, orientation conventions
  engine.py         local synchronous iteration and step rule
  oracle.py         exact fixed point and spectra
  readout.py        power inflow, weights, cancellation, forwarding policy
  metrics.py        routes, accuracy, stretch, loops, stochastic policy
  perturbations.py  cost-change schedule
  scenarios.py      scenarios, parameter grid, curves, dynamic cases
  experiment0.py    static sweep, dynamic run, summaries, CSV/JSON writers
  plots.py
tests/
```

The kickoff proposed one file per field under `fields/`. One `models.py`
replaces them because the models differ only in a coupling block.
Dependencies: `numpy`, `networkx`, `pandas`, `matplotlib`, `pytest`.

## 14. Decision document

After the run, draft `docs/experiment-0.md` answering the charter's five
questions from the recorded outputs, with M3-H and M3-S in separate sections,
counterexamples included, and predictions Pr1–Pr4 marked confirmed or
refuted. It ends with a recommendation; the decision is the owner's.

## 15. Choices open to challenge

- **M1 uses the shared readout.** A bespoke readout for the real control
  (for example the ratio of a cost field to a hop field) might be stronger.
  If either M3 beats M1, the first follow-up is to try to break that result
  with a stronger real readout and with a sign-stripped twin of M3 (the same
  operator with every block entry replaced by its absolute value).
- **The step rule is shared, not optimal per model.** Settle-time
  comparisons are descriptive only.
- **`oracle_hops` is a diagnostic**, not a candidate convention.
- **Asymmetric costs** would give M3-H an orientation for free
  (`φ_ij ∝ c_ij - c_ji` is antisymmetric by construction). Out of scope
  here; noted as a possible revision of the question.
- **Added topologies** (lopsided diamond, ring7) go beyond the kickoff's
  §19 list. They are required by A5, not chosen after seeing results.
