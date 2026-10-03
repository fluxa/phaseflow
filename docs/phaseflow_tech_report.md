# PhaseFlow: can interference in a distributed complex field route packets? A negative result

**Technical report, 2026-10-03. Status: archived.**

**One sentence.** On trees, diamonds and rings, a locally propagated
complex-valued field whose phase encodes link cost does not produce a
forwarding rule that a one-real-per-node attenuation field does not already
produce better; the Hermitian variant is provably blind to cost ordering on
bipartite graphs and steers by a chirality unrelated to cost on odd cycles,
and the symmetric variant's cost response changes sign with its own
parameters.

---

## 1. The question

Quantum mechanics offers structures that classical distributed systems might
borrow: superposition, phase, interference. PhaseFlow asked a narrow version
of that idea. Each node of an ordinary network holds a small complex state
per destination, exchanges it only with neighbours, and the destination
injects a constant source. Route alternatives then reinforce or cancel. Does
the resulting interference pattern carry routing information that a matched
real-valued local field does not?

Two formulations were tested as separate hypotheses.

- **M3-H, Hermitian (chiral) field.** `dψ/dt = -(γI + iH)ψ + b`, with
  `H_ij = a·e^{iφ_ij}`, `H_ji = conj(H_ij)`, and `φ` derived from link cost.
  This is a driven, damped chiral continuous-time quantum walk run on
  classical hardware. Its steady-state probability current is the forwarding
  signal.
- **M3-S, complex-symmetric field.** The same dynamics with
  `H_ij = H_ji = a·e^{iβc_ij}`, so every walk accumulates phase
  `β × walk cost` and two routes interfere according to their cost
  difference. `H` is not Hermitian.

The project made no claim of quantum speed-up, entanglement or novelty in
advance. A complex number is two reals; any benefit had to come from the
interference structure, not from the extra state. The working rule, carried
over from an earlier project, was: narrow hypothesis, nearest prior work,
strongest fair baseline, matched controls, seeded experiments, predeclared
grid, and an honest archive if the mechanism fails.

## 2. Models, controls and the shared readout

Every model was built as the same linear local system so that nothing but
the per-edge coupling differs:

```text
dx_i/dt = -γ x_i + Σ_{j∈N(i)} K_ij x_j + b_i          x_i ∈ R^k, k ∈ {1, 2}
```

A complex number `p + iq` is the block `[[p, -q], [q, p]]`. `K_ij` is what
node `i` receives from node `j`. The local iteration is forward Euler with a
node-local step `α_i = γ/(γ² + s_i²)`, `s_i = Σ_j ‖K_ij‖`; it uses nothing a
node does not already have. An exact oracle (`x* = (γI − K)⁻¹b`, spectra)
was used only to validate the local dynamics, never as an input to them.

| Model | State | Block `K_ij` | Cost enters via |
|---|---|---|---|
| M0 | 1 real | `a·e^{−κc}` | attenuation |
| M1 | 2 reals | `a·diag(e^{−κc}, e^{−2κc})` | attenuation, two scales |
| M2 | complex | `a·R(−π/2)` | not at all (phase disabled) |
| M3-H | complex | `a·R(−π/2 + σ_ij βc)`, `σ_ji = −σ_ij` | loop flux |
| M3-S | complex | `a·R(θ₀ + βc)`, same both ways | walk phase |

M0 is cost-weighted diffusion: its fixed point sums all walks weighted by
`e^{−κ·cost}`, a soft-min of path cost and the real-valued image of M3-S
under `iβ → −κ`. M1 is the matched two-real control: same state size,
payload and per-edge information as the complex models, but non-negative and
diagonal, so contributions can only add. M2 is M3 with `β = 0`: coherent
dynamics, no cost information. M3-H needs a sign convention `σ` for each
undirected edge; eight conventions were run, one of them an explicitly
labelled oracle that uses global destination knowledge.

**One readout for every model.** From the exact identity

```text
d|x_i|²/dt = −2γ|x_i|² + Σ_j P_i←j + 2 x_iᵀb_i,      P_i←j = 2 x_iᵀ K_ij x_j,
```

`P_i←j` is the rate at which neighbour `j` feeds node `i`. For Hermitian
models it is the probability current `2 Im(H_ij ψ_i* ψ_j)`, antisymmetric,
flowing outward from the destination. For M3-S and the real models it is a
reinforcement measure, not a current. A packet at `i` moves upstream, to
`argmax_j max(P_i←j, 0)`. Two-node, balance, antisymmetry and tree-gauge
properties of this readout are unit-tested. One correction to the kickoff
brief came out of this: positive `J_ij` is flow from `j` into `i`, not the
reverse; the forwarding rule survives because packets travel against the
field's flow.

## 3. Two results obtained before the experiment ran

**Bipartite flux-sign blindness (M3-H).** Let `D = diag(±1)` by sublattice
of a bipartite graph. Then `D H* D = −H*`, so if `ψ` solves
`(γ + iH)ψ = b` with real `b`, then `Dψ*` solves the system with every phase
reversed, up to a global sign. Since `D_iD_j = −1` on every edge, every
current `J_ij` is unchanged. Hence on a bipartite graph the steady-state
currents of the driven, damped chiral walk are even in the loop flux: the
sign of the flux, which is the only thing that could distinguish "upper arm
cheaper" from "lower arm cheaper", is invisible. The result reproduces, for
this driven steady state, what Zimborás et al. (2013) and Lu et al. (2016)
proved for transition probabilities. Numerically it holds to `1e-15`.

**Corollary: the standard diamond cannot test M3-H.** The kickoff's diamond
(arms `A–B–D–T` and `A–C–E–T`) has one six-cycle and a mirror symmetry
exchanging the arms. By gauge invariance only the flux `Φ` matters; mirror
symmetry gives `split(Φ) = 1 − split(−Φ)`, evenness gives
`split(Φ) = split(−Φ)`, so node A splits exactly 50/50 for every cost
assignment, `β`, `γ` and orientation. The experiment therefore added odd-cycle
graphs: a lopsided diamond (arms `A–B–T` and `A–C–E–T`, a 5-cycle) and a
7-ring.

A smaller technical result: with the node-local step rule above, the
Hermitian iteration converges for every `γ > 0` on every graph (each scaled
eigenvalue `μ` satisfies `|1 − μ|² ≤ 1 − Re μ`); for non-negative real
couplings it is a regular splitting and converges whenever the continuous
system is stable. M3-S has no such guarantee.

## 4. Experiment 0

**Scenarios.** Fifteen, one destination each: path and binary tree with unit
and seeded random costs (4); the diamond with equal, slightly unequal and
strongly unequal arms (4); the lopsided diamond with the three-hop arm dearer,
cheaper and nearly tied (3); rings of 8 and 7 nodes with unit costs and with
one heavy edge (4).

**Grid, declared before the run.** `γ ∈ {0.1, 0.3, 1, 2, 3.5}`;
`κ` or `β ∈ {0.25, 0.5, 1, 2, 4}`; eight M3-H orientations; two M3-S values
of `θ₀` (`−π/2`, the kickoff form, and `0`). 61 model cells × 5 `γ` × 15
scenarios = 4575 cells at the exact fixed point; 692 are unstable and carry
no readout. Plus discrimination curves (node A's split as one arm's cost is
swept) and a dynamic run through "raise one cost, then restore it" with the
local iteration.

**Metrics.** Next-hop accuracy against Dijkstra, delivery, loops, dark nodes,
stretch, a stochastic-policy expected stretch, a cancellation measure, and
stability flags from the oracle. Four predictions were written down before
the run.

**What is post-hoc.** After the run, two things changed and are labelled as
such throughout. Ties are scored fairly (a node whose top weights tie scores
the share of tied hops that are optimal), because the planned lowest-id
tie-break let node numbering decide some tallies. And a set of comparisons
(against the cost-blind M2, on "cost-sensitive" scenarios, by loop flux) was
added; `src/phaseflow/analysis.py` regenerates all of them. The orientation
set also turned out to contain duplicates: `random:4` equals `id` on graphs
with at most seven edges, and `id` equals the oracle on the diamond graphs.
Because on a single-cycle graph only the loop flux matters, the 1925 M3-H
scenario-cells reduce to 995 distinct ones, and the oracle sits inside the
M3-H counts on seven of eleven cyclic scenarios. This favours M3-H; the
counts below are as run.

## 5. Results

**Controls behaved as derived.** M3-H equals M2 on every tree to `1e-15`
(gauge invariance); M3-H splits the diamond 50/50 to `7e-15`; M3-H and M2 are
stable in every cell, M3-S in 340 of 750, M0 and M1 in 234 of 375.

**Headline tallies** (tie-fair accuracy; "cost-sensitive" means the six
cyclic scenarios on which the cost-blind M2 gets a node wrong):

| | Cells correct on every cyclic scenario | Fully correct on cost-sensitive scenarios | Against M2 at the same `γ`: better / equal / worse |
|---|---|---|---|
| M0 | 5 of 13 usable | 70 of 93 | 81 / 81 / 9 |
| M1 | 5 of 13 usable | 71 of 93 | 82 / 80 / 9 |
| M3-H | 0 of 175 | 15 of 1050 | 15 / 1743 / 167 |
| M3-H, oracle orientation | 0 of 25 | 4 of 150 | 4 / 251 / 20 |
| M3-S, `θ₀ = −π/2` | 0 of 11 usable | 14 of 69 | 19 / 78 / 30 |
| M3-S, `θ₀ = 0` | 0 of 10 usable | 9 of 65 | 11 / 73 / 33 |

M0's worst usable cell (mean accuracy 0.954) equals M3-S's best. At `γ ≥ 2`,
where only 7 of 1342 cyclic scenario-cells are unstable, mean accuracy is
M1 0.981, M0 0.979, M2 0.924, M3-H 0.924, M3-S 0.921 and 0.909, so the
ordering is not an artefact of dropping unstable cells.

**M3-H.** At `γ ≥ 1` its accuracy equals M2's in all 1320 cyclic
scenario-cells: the damping that makes it fast removes any effect of phase
on decisions. At `γ ≤ 0.3` phase does change decisions, and the 15 cells
that are fully correct on a cost-sensitive scenario are all on the heavy
7-ring at `γ = 0.3`. They are a chirality effect:

| Sign of loop flux | Cells | Accuracy at `γ = 0.3` | at `γ = 0.1` |
|---|---|---|---|
| positive | 19 | 1.000 | 0.833 |
| zero | 5 | 0.833 (= M2) | 0.833 |
| negative | 16 | 0.500 – 0.667 | 0.333 |

Moving the heavy edge to its mirror image with the orientations unchanged
reverses the table: 23 of 23 negative-flux cells become perfect and no
positive-flux cell is. The sign is set by the orientation and by how
`β × cost` wraps around `2π`, not by which side is dear. On the lopsided
diamond a scan over every loop flux, `γ` from 0.02 to 3.5 and lower-arm cost
from 0.3 to 1.5 gives a smallest upper-arm share of 0.646 at node A: M3-H
cannot prefer a three-hop arm however cheap it is. Its worst route is a node
adjacent to the destination on the heavy 7-ring forwarding the long way
round, stretch 9. It produced no forwarding loop and no dark node in 3000
cells.

**M3-S.** Usable only at strong damping (19 of its 21 usable cells have
`γ ≥ 2`), and 15 stable cells cannot be reached by the local iteration at
the shared step. It does respond to cost on equal-hop arms, but the response
flips sign with its own parameters:

| `θ₀` | `β` | Curve points on the correct side (of 20) | Range of A's upper-arm share |
|---|---|---|---|
| `−π/2` | 0.5 | 20 | 0.495 – 0.506 |
| `−π/2` | 1 | 20 | 0.491 – 0.649 |
| `−π/2` | 2 | 0 | 0.000 – 0.567 |
| `0` | 0.5 | 0 | 0.488 – 0.507 |
| `0` | 1 | 0 | 0.386 – 0.506 |
| `0` | 2 | 20 | 0.451 – 1.000 |

Where hop counts differ it never prefers the cheaper arm (upper-arm share
never below 0.72 on the lopsided diamond). All 23 of its fully correct
cost-sensitive cells are on diamond scenarios; none is on the lopsided
diamond or either heavy ring. Cancellation is present (mean 0.21 per node)
and does not track accuracy.

**Dynamics.** Every stable run re-converges after the cost change and
returns to its original routing after the restore; no stable run forms a
forwarding loop in any round. M3-H settles in 803–905 rounds at `γ = 0.3`,
96–111 at `γ = 1`, 16–20 at `γ = 3.5`; the regime where phase matters is the
slow one. M3-S diverges at `γ = 0.3` and `γ = 1` and settles in 19–36 rounds
at `γ = 3.5`. M0 at `κ = 1` diverges at `γ = 0.3` and fails the lopsided cost
change at `γ = 3.5` (hop bias) while passing it at `γ = 1`.

**Predictions.** Pr1 (every stable model routes trees correctly) was
refuted for M0, M1 and M3-S: M0 forms two-node loops within 0.06 of its
stability edge. Pr2 (M0 correct at large `κ`, hop-biased at small `κ`) held,
with large `γ` also causing hop bias. Pr3 (M3-S discriminates weakly with a
rule-dependent sign) held except that at `β = 2` the response is not weak.
Pr4 (no M3-H orientation is correct across costs) held and was sharpened to
the flux-sign statement above.

## 6. What favoured the complex models

Two findings go the other way and are part of the record.

- **Proportional forwarding.** Under `p(i→j) ∝ max(P, 0)`, mean expected-cost
  stretch over stable cyclic cells is M2 1.035, M3-H (`γ ≥ 1`) 1.040,
  M3-S 1.14–1.24, M1 1.253, M0 1.266. M0 is worse than M2 in 128 of 171
  scenario-cells. The reason is structural: a current is antisymmetric, so a
  downstream neighbour gets zero weight, whereas M0's inflows are all
  positive and some traffic flows backwards. The cost-blind M2 is best of
  all, so the advantage belongs to the Hermitian structure, not to
  cost-derived phase; it may equally be read as the shared readout
  handicapping the real controls under this policy.
- **One M3-S rule is consistent on equal-hop arms.** `θ₀ = −π/2`, `β = 1` is
  cost-correct on all four diamond scenarios at `γ = 2` and `3.5`, on all 20
  curve points, and in the dynamic diamond. The same rule scores 0.800,
  0.833 and 0.857 on the three scenarios where hop count and cost disagree,
  where M0's good cells score 1.

Neither rescues a branch: the first is not a phase effect, and the second
is confined to one parameter pair and to graphs where the arms have equal
hop counts.

## 7. Why it fails: mechanism

- **M3-H.** Symmetric costs enter a chiral field only through signed loop
  fluxes, which need an edge orientation. On bipartite graphs the sign is
  provably invisible. On odd cycles the flux does steer forwarding, but by a
  chirality that depends on the orientation and on `2π` wrap-around, not on
  which route is cheaper; the mirrored ring shows the same cells going from
  right to wrong when only the dear side moves. Hop count sets the magnitude
  structure and phase can only modulate it, never override it, which is why
  the three-hop arm is never chosen at any flux.
- **M3-S.** Two arms with equal hop count deliver equal-magnitude
  contributions whose relative phase is `β × Δcost`. The readout projects
  each contribution onto the sum, and to first order that projection is even
  in the phase difference: it reports `|Δcost|` modulo `2π/β`, not its sign.
  The sign information the curves do show comes from higher-order walks and
  from the per-hop phase `θ₀`, which is why it flips with `β` and `θ₀`.
  Reading the phase itself instead would recover `β × path cost` for a
  dominant path, but that is a distance-vector protocol carried in phase, in
  which interference is noise rather than signal; it was not tested and is a
  different hypothesis.
- **The real control.** M0's fixed point is a sum over walks weighted by
  `e^{−κ·cost}`, a soft-min that approaches the shortest-path next hop as
  `κ` grows, at the price of needing `γ` above its spectral radius and of a
  hop bias of `ln(γ/a)/κ` per hop when `γ ≠ a`. That is the known
  randomized-shortest-path / logit-assignment construction; M3-S is its
  continuation to imaginary inverse temperature, and the experiment says the
  continuation loses the ordering information the real version keeps.

## 8. Prior work and what remains

The nearest prior work (details and URLs in `docs/literature.md`):

- Zimborás et al. (2013) and Lu et al. (2016): the chiral-walk Hamiltonian,
  with the tree and bipartite gauge results that §3 reproduces for the
  driven steady state.
- Matsuoka, Ohno, Segawa (2024): the stationary state of a sourced,
  absorbing quantum walk used as a shortest-path finder, proposed as a
  quantum-inspired classical algorithm, with weaker selectivity once cycles
  offer alternative routes; the same weakening seen here.
- Higuchi, Segawa (2019); Higuchi, Sabri, Segawa (2020): a constantly driven
  walk whose stationary state is the classical electrical current.
- Akamatsu (1996); Saerens et al. (2009); Kivimäki et al. (2014);
  Todorov (2006): the real control's literature.
- Ranzani and Aumentado (2015); Metelmann and Türeci (2018): driven-damped
  coupled-mode networks with complex couplings, the M3-H equation in
  classical hardware with phases designed by hand.

No positive novelty claim remains. The negative one, that steady-state
current in a driven, damped chiral walk cannot order routes by symmetric
cost, was not found stated anywhere, but it follows closely from known gauge
results and is demonstrated here only at smoke-test scale.

## 9. Limitations

- One readout. The phase-reading alternative described in §7 was not run.
- One cycle per graph, at most eight nodes, five values each of `γ` and
  `β`/`κ`. The lopsided M3-H result is not grid-limited (full flux scan);
  the others are.
- The orientation set contained duplicates and an oracle coincidence (§4),
  both favouring M3-H.
- M1 adds nothing to M0 under the shared readout, so it is not a stronger
  baseline, and the shared readout may handicap the real controls under
  proportional forwarding.
- Stability was judged on the continuous system; 15 stable M3-S cells are
  unreachable by the local iteration as configured.
- Exploratory throughout: only the grid and four predictions were fixed in
  advance, and the headline comparisons are post-hoc and labelled as such.
- The spec and plan were implemented without a separate owner review round;
  an independent reviewer later checked the branch and its findings, and
  corrections from that review are recorded in `docs/experiment-0.md` §9.

## 10. Decision

Both branches are stopped as formulated and the project is archived. The
kickoff's kill criteria met at this scale: no transfer (A); a simpler real
method matches or beats it on deterministic forwarding (B, E); destructive
interference at low damping is mostly harmful for M3-H (C, 167 worse against
15 better); and for M3-H the damping needed for speed removes the feature
entirely (D).

The one revision considered and not pursued is asymmetric link costs in
M3-H, where a phase proportional to `c_ij − c_ji` needs no orientation and
chirality would become a function of cost. It would help only on
non-bipartite graphs and only at low damping, and a real field with
asymmetric weights is the obvious control. The odds were judged low.

## Appendix: reproducing this report

```bash
uv sync
uv run pytest -q                           # 72 tests
uv run python -m phaseflow.experiment0     # predeclared sweep, about 8 s
uv run python -m phaseflow.analysis        # post-hoc tables -> posthoc.json
```

| Where | What |
|---|---|
| `docs/phaseflow_kickoff.md`, `docs/phaseflow_research_charter.md` | the original question and the working rules |
| `docs/superpowers/specs/2026-10-02-experiment0-smoke-design.md` | the design, including the pre-run derivations and predictions |
| `docs/experiment-0.md` | the full findings, with every number tied to a results file |
| `docs/literature.md` | prior-art note with retrieved URLs |
| `src/phaseflow/` | models, engine, oracle, readout, metrics, scenarios, runner, post-hoc analysis |
| `results/experiment0-smoke/` | `config.json`, `summary.json`, `posthoc.json`, CSVs, plots |

The run is deterministic; `config.json` records the grid, seeds, package
versions and the code commit.
