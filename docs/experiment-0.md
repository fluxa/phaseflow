# Experiment 0 (smoke test) — findings and open decision

**Date:** 2026-10-02
**Status:** Exploratory. Draft for the project owner; the decision in §10 is
not taken here.
**Design:** `docs/superpowers/specs/2026-10-02-experiment0-smoke-design.md`
**Data:** `results/experiment0-smoke/`, produced by code at commit `e6dd4e8`.
Every number below comes from `route_metrics.csv`, `transfer.csv`,
`discrimination.csv`, `dynamic_metrics.csv`, `summary.json` or `posthoc.json`
in that directory.

This document was revised after an independent review of the branch. The
review found two wrong statements, duplicated orientations, and tallies that
depended on node numbering. §9 lists what changed.

## Summary

Neither formulation produced a general rule that turns link cost into phase
and yields correct forwarding.

- No M3-H or M3-S parameter cell forwards correctly on every cyclic scenario.
  Five cells of the real-valued scalar control M0 do.
- On the six scenarios where a cost-blind field gets a node wrong, M3-H is
  fully correct in 15 of 1050 stable cells, M3-S in 23 of 134, and M0 in 70
  of 93.
- Where cost-derived phase changes a forwarding decision, it makes it worse
  more often than better: 167 against 15 for M3-H, 63 against 30 for M3-S.

Two things in favour of the complex models are real and are reported in §3
to §6: under proportional forwarding the Hermitian current wastes less than
M0, and one M3-S rule is cost-correct on equal-hop arms in every test. Neither
rescues a branch, for reasons given there.

M3-H and M3-S are reported separately and are never combined.

## 1. What was run, and how it is scored

- **Static sweep.** 15 scenarios (4 tree, 4 diamond, 3 lopsided diamond, 4
  ring) × 61 model cells × 5 values of `γ` = 4575 cells, each evaluated at the
  exact fixed point. 692 cells are unstable and carry no readout.
- **Discrimination curves.** Node A's forwarding split as one lower-arm cost
  is swept, on the diamond (21 points) and the lopsided diamond (25 points).
- **Dynamic run.** Local iteration from zero through "raise one cost, then
  restore it", 2500 rounds per phase, on the diamond and lopsided diamond.
- **Tests.** 72 unit tests pass, including the sign, balance, gauge and
  symmetry conventions that gate the sweep.

**Accuracy is tie-fair.** A node whose top weights tie scores the share of
tied hops that are optimal. The planned rule (lowest node id wins) made M2
and M3-H right on two diamond scenarios and wrong on a third purely by
numbering. Figures under the planned rule are in `posthoc.json` under keys
ending `_lowest_id_ties`.

**Which analyses are post-hoc.** The grid, the metrics and predictions
Pr1–Pr4 were fixed before the run. Everything in `posthoc.json` was defined
afterwards: tie-fair scoring, the "cost-sensitive" subset, the comparisons
against M2, the flux-sign tables, the mirrored ring and the flux scan.
`uv run python -m phaseflow.analysis` regenerates it.

**"Cost-sensitive" scenarios** are the cyclic scenarios on which the
cost-blind model M2 gets at least one node wrong: `diamond_lower+0.2`,
`diamond_lower+1.0`, `diamond_upper+1.0`, `lop_lower_cheap`, `ring8_heavy`,
`ring7_heavy`.

## 2. Control checks

| Check | Result |
|---|---|
| A1 sign: on two nodes, the non-destination node forwards to the destination | holds for all five models (unit test) |
| A2 balance identity at the fixed point | holds for all models (unit test) |
| A3 tree gauge invariance, M3-H vs M2 | largest relative difference in `P` over all 800 M3-H tree cells: `1.3e-15` |
| A4 bipartite flux-sign blindness, M3-H | exact on ring-8 and the diamond; fails on ring-7 and the lopsided diamond, as it should (unit test) |
| A5 diamond split, M3-H | largest deviation from 0.5 over all 800 M3-H diamond cells: `7.3e-15` |
| A6 stability | M3-H and M2 stable in every cell (3000 of 3000, 75 of 75); M3-S in 340 of 750; M0 and M1 in 234 of 375 each |

The local iteration converges in every stable M0, M1, M2 and M3-H cell. It
diverges in 15 of the 340 stable M3-S cells, all at `γ = 1`: the shared step
rule is too large for M3-S near its stability edge. Those cells are still
evaluated at their fixed point.

## 3. M3-H — Hermitian / chiral field

**Forwarding interpretation: coherent.** `P` is a true edge current. M3-H is
stable for every `γ`, produced no forwarding loop and no dark node in 3000
cells, and forwards every node correctly on every tree (800 of 800 cells),
exactly as M2 does.

**The orientation set was weaker than planned.** On a graph with one cycle
only the wrapped loop flux `Φ` matters (unit test), and the eight
conventions turned out not to be eight different things:

- `random:4` is edge-for-edge identical to `id` on every graph with seven or
  fewer edges, and `random:2` to `random:3` on the lopsided diamond.
- `id` is identical to the diagnostic `oracle_hops` on all seven diamond and
  lopsided scenarios, because those graphs are numbered toward the
  destination.

So the 1925 "arbitrary orientation" scenario-cells contain 995 distinct
`(scenario, Φ, γ)` cells, the oracle orientation sits inside the M3-H counts
on seven of eleven cyclic scenarios, and the M3-H dynamic run and the `id`
discrimination curve used the oracle orientation. All of this favours M3-H.
The counts below are as run, duplicates included.

**Cost-to-phase rule: none found.**

- *Diamond.* Node A splits exactly 50/50 for every cost, `β`, `γ` and
  orientation. This is the theorem in spec A5, confirmed numerically.
- *All cyclic scenarios.* None of the 175 cells is correct on every cyclic
  scenario. Mean accuracy per cell ranges from 0.724 to 0.924. The maximum is
  M2's value, so no M3-H cell beats the cost-blind model on average.
- *Against M2 at the same `γ`.* Over 1925 scenario-cells: better in 15, equal
  in 1743, worse in 167. At `γ ≥ 1` accuracy equals M2's in all 1320 cyclic
  scenario-cells (oracle included): the damping that makes M3-H fast removes
  any effect of phase on decisions. This is the kickoff's Kill D.
- *Cost-sensitive scenarios.* Fully correct in 15 of 1050 cells, all on
  `ring7_heavy` at `γ = 0.3`.
- *Those 15 cells are a chirality effect, not a cost effect.* On
  `ring7_heavy` at `γ = 0.3`, accuracy is decided by the sign of `Φ`:

  | Sign of `Φ` | Cells (with oracle) | Accuracy at `γ = 0.3` | at `γ = 0.1` |
  |---|---|---|---|
  | positive | 19 | 1.000 | 0.833 |
  | zero | 5 | 0.833 (M2's value) | 0.833 |
  | negative | 16 | 0.500 – 0.667 | 0.333 |

  The sign is set by the orientation and by how `β × cost` wraps around
  `2π`, not by which side of the ring is dear. Moving the heavy edge to its
  mirror image, with the orientations unchanged, reverses the table: 23 of 23
  negative-flux cells are perfect and no positive-flux cell is.
- *Lopsided diamond.* A scan over every loop flux, `γ` from 0.02 to 3.5 and
  lower-arm cost from 0.3 to 1.5 gives a smallest upper-arm share of 0.646.
  By gauge invariance this covers every orientation and every `β`: M3-H
  cannot prefer the three-hop arm at node A, however cheap it is. It is fully
  correct on `lop_lower_cheap` in 0 of 200 cells.

**Oracle orientation (diagnostic).** 4 of 150 cells correct on the
cost-sensitive scenarios; against M2, better in 4, equal in 251, worse in 20.
The orientation problem is not the only obstacle.

**Counterexample.** On `ring7_heavy` at `γ = 0.1`, a node adjacent to the
destination forwards the long way round, across the heavy edge: a path-cost
stretch of 9. Stretch exceeds 5 in 88 M3-H cells and in no cell of any other
model.

**Speed.** In the dynamic run M3-H settles in 803–905 rounds at `γ = 0.3`,
96–111 at `γ = 1` and 16–20 at `γ = 3.5`. Phase affects decisions only in the
slow regime.

## 4. M3-S — complex symmetric / accumulated path-phase field

**Forwarding interpretation: valid but weaker.** `P` is not antisymmetric
(unit test), so it is a reinforcement measure, not a current. It forwards
correctly on trees in 95 of 96 stable cells; the exception loops
(`btree7_leaf`, `β = 2`, `γ = 1`).

**Stability.** Stable in 340 of 750 cells. Of the 50 parameter cells, 21 are
stable on every cyclic scenario, 19 of them at `γ ≥ 2`. In the dynamic run it
diverges at `γ = 0.3` and `γ = 1` and settles in 19–36 rounds at `γ = 3.5`.

**Cost-to-phase rule: none found.**

- *All cyclic scenarios.* None of the 21 usable cells is correct on every
  cyclic scenario. Mean accuracy ranges from 0.880 to 0.954; M2's is 0.924.
- *Against M2 at the same `γ`.* `θ₀ = -π/2`: better in 19, equal in 78, worse
  in 30. `θ₀ = 0`: better in 11, equal in 73, worse in 33.
- *Cost-sensitive scenarios.* Fully correct in 23 of 134 stable cells (14 of
  69 and 9 of 65). All 23 are on the three diamond scenarios. None is on
  `lop_lower_cheap`, `ring8_heavy` or `ring7_heavy`.
- *Diamond curve.* M3-S responds to cost on equal-hop arms, but the sign of
  the response depends on the rule:

  | `θ₀` | `β` | swept points on the correct side | range of A's upper-arm share |
  |---|---|---|---|
  | `-π/2` | 0.5 | 20 of 20 | 0.495 – 0.506 |
  | `-π/2` | 1 | 20 of 20 | 0.491 – 0.649 |
  | `-π/2` | 2 | 0 of 20 | 0.000 – 0.567 |
  | `0` | 0.5 | 0 of 20 | 0.488 – 0.507 |
  | `0` | 1 | 0 of 20 | 0.386 – 0.506 |
  | `0` | 2 | 20 of 20 | 0.451 – 1.000 |

  Doubling `β` or changing `θ₀` turns a rule that is right at every point
  into one that is wrong at every point.
- *Lopsided diamond curve.* A never puts less than 0.72 on the two-hop arm.
  Hop count dominates cost.

**The most favourable M3-S result.** The rule `θ₀ = -π/2`, `β = 1` is
cost-correct on equal-hop arms in every test: all four diamond scenarios at
`γ = 2` and `γ = 3.5`, all 20 curve points, and the dynamic diamond, where it
re-routes correctly after the cost rise and returns after the restore. The
same rule scores 0.800 on `lop_lower_cheap`, 0.833 on `ring7_heavy` and 0.857
on `ring8_heavy`, where M0's good cells score 1.

**Interference is present but not informative.** Mean cancellation per node
is 0.21 across stable cyclic cells (0.49 for M2 and M3-H, 0 for M0 and M1).
It does not track accuracy (`plots/cancellation_vs_accuracy.png`).

## 5. Controls

**M0, scalar attenuation.** Stable on every cyclic scenario in 13 of 25
cells. Five are correct on every cyclic scenario: `(κ, γ)` = (2, 1), (2, 2),
(4, 1), (4, 2), (4, 3.5). Its worst usable cell has mean accuracy 0.954,
equal to M3-S's best. Fully correct on the cost-sensitive scenarios in 70 of
93 stable cells. Against M2: better in 81, equal in 81, worse in 9.

M0's weaknesses, for the record:

- It is unstable at low `γ` unless `κ` is large.
- At strong damping it is biased toward fewer hops. No M0 cell on the
  lopsided curve is both stable and correct on both sides of the crossover.
  It fails the lopsided cost change in the dynamic run at `γ = 3.5` and
  passes it at `γ = 1`.
- Within 0.06 of its stability edge it forms forwarding loops: 5 stable cells,
  4 of them on trees.
- Under proportional forwarding it wastes more than the Hermitian models
  (§6).

**M1, two real states.** Its accuracy differs from M0 in one row of 375.
Under the shared readout the second real state adds nothing. M1 shows that
extra state alone does not help, but it is not a stronger baseline than M0.

**M2, phase disabled.** Cost-blind by construction. Mean accuracy 0.924 at
every `γ`; correct exactly where hop count and cost agree.

**At matched damping.** At `γ ≥ 2`, where only 7 of 1342 cyclic scenario-cells
are unstable, mean accuracy is M1 0.981, M0 0.979, M2 0.924, M3-H 0.924,
M3-S 0.921 (`θ₀ = -π/2`) and 0.909 (`θ₀ = 0`). The ordering is not an
artefact of dropping unstable cells.

## 6. Proportional forwarding and the dynamic run

**Proportional forwarding** (`p(i → j) ∝ max(P, 0)`). Every stable cell
delivers with probability 1. Mean expected-cost stretch over stable cyclic
cells:

| M2 | M3-H, `γ ≥ 1` | M3-S, `θ₀ = 0` | M3-H, `γ < 1` | M3-S, `θ₀ = -π/2` | M1 | M0 |
|---|---|---|---|---|---|---|
| 1.035 | 1.040 | 1.141 | 1.224 | 1.235 | 1.253 | 1.266 |

M0 has higher expected stretch than M2 in 128 of 171 scenario-cells. So the
statement "M0 is better on route quality" holds for the deterministic next
hop only. The reason is structural: a current is antisymmetric, so a
downstream neighbour gets zero weight, whereas M0's inflows are all positive
and some traffic flows backwards. The cost-blind M2 has this property too and
is the best of all, and M3-H is higher than M2 in 697 scenario-cells and
lower in 153. The advantage belongs to the Hermitian structure, not to
cost-derived phase. It may also mean the shared readout handicaps the real
controls under this policy.

**Dynamic run: accuracy after the cost rise** (lowest-id ties, as logged).
"div." means the run diverged.

| | diamond `γ=0.3` | `γ=1` | `γ=3.5` | lopsided `γ=0.3` | `γ=1` | `γ=3.5` |
|---|---|---|---|---|---|---|
| M0, M1 (`κ = 1`) | div. | 1.000 | 1.000 | div. | 1.000 | 0.800 |
| M2 | 0.833 | 0.833 | 0.833 | 0.800 | 0.800 | 0.800 |
| M3-H (oracle orientation) | 0.833 | 0.833 | 0.833 | 0.600 | 0.800 | 0.800 |
| M3-S, `θ₀ = -π/2` | div. | div. | 1.000 | div. | div. | 0.800 |
| M3-S, `θ₀ = 0` | div. | div. | 0.833 | div. | div. | 0.800 |

Every stable run returns to its original routing after the restore, and no
stable run has a forwarding loop in any round. M2 and M3-H's 0.833 on the
diamond is the exact 50/50 tie broken toward the wrong arm. Their
`route_settle_rounds` there (up to 1396) measure how long a floating-point
tie takes to fall inside the tie tolerance, not a routing transient.

## 7. Predictions

| | Prediction | Outcome |
|---|---|---|
| Pr1 | Every stable model forwards every tree node correctly | **Refuted** for M0 and M1 (59 of 63 stable tree cells) and M3-S (95 of 96). Confirmed for M2 and M3-H. |
| Pr2 | M0 is correct at large `κ`, hop-biased at small `κ` | Confirmed, with the addition that large `γ` also causes hop bias. |
| Pr3 | M3-S on equal-hop arms discriminates weakly, with a sign that depends on `θ₀` and `β` | Sign dependence confirmed. "Weakly" is wrong at `β = 2`, where the split reaches 0 or 1. |
| Pr4 | M3-H on odd cycles is orientation-dependent and no arbitrary orientation is correct across costs | Confirmed, and sharpened: the dependence is on the sign of the loop flux, which is unrelated to which route is cheaper. |

## 8. The charter's five questions

**1. Which branch has a coherent general cost-to-phase rule and forwarding
interpretation?**

- *M3-H:* the forwarding interpretation is coherent. There is no cost-to-phase
  rule. On bipartite graphs none can exist for symmetric costs (A4). On odd
  cycles the flux biases forwarding by chirality, and the chirality is not a
  function of which route is cheaper.
- *M3-S:* the readout is well defined but is not a current, and the model is
  usable only at strong damping. No rule transfers: the sign of the cost
  response flips with `β` and `θ₀`, and hop count overrides cost on unequal
  arms.

**2. Is any advantage tied to phase or interference rather than state,
damping or a weak baseline?**

None was found that is tied to cost-derived phase. M3-H's stability, freedom
from loops and low proportional-forwarding stretch are all shared by M2, so
they come from the Hermitian structure. M3-S's one consistent success is
confined to equal-hop arms and to one `(θ₀, β)`.

**3. What did the tree, diamond and ring reveal?**

- *Tree:* M3-H is exactly M2. M0 loops near its stability edge.
- *Diamond:* M3-H is exactly cost-blind. M3-S responds to cost with a
  rule-dependent sign.
- *Lopsided diamond:* neither branch ever prefers a cheaper arm that has more
  hops, and for M3-H this holds for every flux.
- *Ring:* on the even ring M3-H's accuracy is identical to M2's in every
  cell. On the odd ring the flux sign decides the outcome, and stretch
  reaches 9.

**4. Closest prior work and remaining novelty claim?**

See `docs/literature.md`. The M3-H tree and bipartite results follow from
Zimborás et al. (2013) and Lu et al. (2016). The real control is the
randomized-shortest-path and logit-assignment literature. The closest
positive claim in print, Matsuoka et al. (2024), reports the same weakening
once cycles offer alternative routes. No positive novelty claim remains. A
possible negative one is that steady-state current in a driven, damped chiral
walk cannot order routes by symmetric cost. The search found no paper stating
it, but it is close to a corollary of known results.

**5. Stop, revise once, or preregister a larger test?**

See §10.

## 9. Limitations, deviations and review changes

**Limitations.**

- **One readout.** Only the power-inflow readout was tested. A readout on the
  phase of each neighbour's contribution was not. For a single dominant path
  that phase is `β × path cost`, so it would be a distance-vector protocol
  carried in phase, with interference as noise rather than signal. That is a
  different hypothesis.
- **One cycle per graph, at most eight nodes.**
- **Coarse grid** for M3-S and for M3-H on rings. The lopsided M3-H result is
  not grid-limited (flux scan, §3).
- **Orientations were not independent** (§3).
- **M1 is not a strong control** (§5), and the shared readout may handicap M0
  and M1 under proportional forwarding (§6).
- **Stability is judged on the continuous system.** Fifteen stable M3-S cells
  cannot be reached by the local iteration as configured.
- **Settle times are descriptive.** The step rule is shared, not tuned.
- **Exploratory.** Only the grid and Pr1–Pr4 were fixed in advance.

**Deviations from the spec.**

- No "state magnitude over time" plot was produced (spec §12).
- Settle time is reported in rounds only; the per-run mean step `alpha_mean`
  is recorded but intrinsic time `Σα` is not (spec §9).
- Tie-fair accuracy, the loop-flux column and `analysis.py` were added after
  review.
- Three large CSVs (`states.csv`, `currents.csv`, `static_readouts.csv`) are
  not committed; they regenerate in eight seconds.
- The spec and plan were not separately reviewed by the owner before the run.

**Changed after the independent review.**

- "Scattered with no pattern" for the 15 correct M3-H cells was wrong; the
  pattern is the flux sign (§3).
- "Flux reversed gives a split of 1.0" described the wrong curve and is
  replaced by the flux scan (§3).
- Duplicate orientations and the oracle coincidence are now stated with their
  consequences (§3).
- Tallies that depended on node numbering are replaced by tie-fair ones. The
  cost-sensitive set grew from four scenarios to six.
- Proportional-forwarding and dynamic routing outcomes, both in the spec,
  were missing and are now in §6.
- "M0 matches or beats both branches on route quality" is restricted to the
  deterministic next hop.

## 10. Recommendation

Stop both branches as formulated, and archive Experiment 0 as a negative
result. Kill criteria from the kickoff met at smoke-test scale:

- **Kill A, no transfer:** no cell of either branch is correct across the
  cyclic scenarios.
- **Kill B and E, a simpler real method matches it:** on deterministic
  forwarding M0 is better with one real per node.
- **Kill C, destructive interference is mostly pathological:** M3-H at low
  damping is worse than M2 in 167 scenario-cells and better in 15.
- **Kill D, damping removes the feature:** for M3-H, at `γ ≥ 1` phase changes
  no decision.

If one revision is wanted before archiving, the candidate with a specific
rationale is **asymmetric link costs in M3-H**. A phase proportional to
`c_ij − c_ji` is antisymmetric by construction, so it needs no orientation
convention, and the chirality it produces would then be a function of cost.
The flux-sign result in §3 shows that chirality does steer forwarding on odd
cycles when it has the right sign. Three cautions: it helps only on
non-bipartite graphs, only at low damping where M3-H is slow, and a real
field with asymmetric weights is the obvious control and may match it. I
would put the odds of a positive result as low.

I do not recommend a larger held-out test of the current formulations.

The decision is the project owner's.
