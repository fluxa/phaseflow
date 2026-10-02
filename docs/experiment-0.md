# Experiment 0 (smoke test) — findings and open decision

**Date:** 2026-10-02
**Status:** Exploratory. Draft for the project owner; the decision in §9 is
not taken here.
**Design:** `docs/superpowers/specs/2026-10-02-experiment0-smoke-design.md`
**Data:** `results/experiment0-smoke/` (code at commit `da0025a`; every
number below comes from the files there)

## Summary

Neither formulation produced a general rule that turns link cost into phase
and yields correct forwarding. On the four scenarios where a cost-blind field
gets a node wrong, M3-H is fully correct in 15 of 700 stable cells and M3-S
in 7 of 91; the real-valued scalar control M0 is fully correct in 41 of 61.
No M3-H or M3-S parameter cell is correct on every cyclic scenario; five M0
cells are. Where cost-derived phase changes a forwarding decision at all, it
makes it worse more often than better.

M3-H and M3-S are reported separately below and are never combined.

## 1. What was run

- **Static sweep.** 15 scenarios (4 tree, 4 diamond, 3 lopsided diamond, 4
  ring) × 61 model cells × 5 values of `γ` = 4575 cells, each evaluated at the
  exact fixed point. 692 cells are unstable and carry no readout.
- **Discrimination curves.** Node A's forwarding split as one lower-arm cost
  is swept, on the diamond (21 points) and the lopsided diamond (25 points).
- **Dynamic run.** Local iteration from zero through "raise one cost, then
  restore it", 2500 rounds per phase, on the diamond and lopsided diamond.
- **Tests.** 58 unit tests pass, including the sign, balance, gauge and
  symmetry conventions that gate the sweep.

Three files are not committed because they are large and regenerate in
eight seconds: `states.csv`, `currents.csv`, `static_readouts.csv`.

## 2. Control checks

| Check | Result |
|---|---|
| A1 sign: on two nodes, the non-destination node forwards to the destination | holds for all five models (unit test) |
| A2 balance identity at the fixed point | holds to `1e-12` for all models (unit test) |
| A3 tree gauge invariance, M3-H vs M2 | largest relative difference in `P` over all 800 M3-H tree cells: `1.3e-15` |
| A4 bipartite flux-sign blindness, M3-H | exact on ring-8 and the diamond; fails on ring-7 and the lopsided diamond, as it should (unit test) |
| A5 diamond split, M3-H | largest deviation from 0.5 over all 800 M3-H diamond cells: `7.3e-15` |
| A6 stability | M3-H and M2 stable in every cell (3000 of 3000, 75 of 75); M3-S in 340 of 750; M0 and M1 in 234 of 375 each |

The local iteration converges in every stable M0, M1, M2 and M3-H cell. It
diverges in 15 of the 340 stable M3-S cells, all at `γ = 1`: the shared step
rule is too large for M3-S near its stability edge.

## 3. M3-H — Hermitian / chiral field

**Forwarding interpretation: coherent.** `P` is a true edge current. M3-H is
stable for every `γ`, produced no forwarding loop and no dark node in 3000
cells, and forwards every node correctly on every tree (800 of 800 cells),
exactly as M2 does.

**Cost-to-phase rule: none found.**

- *Diamond.* Node A splits exactly 50/50 for every cost, `β`, `γ` and
  orientation. This is the theorem in spec A5, confirmed numerically. M3-H's
  accuracy on the diamond is therefore set by the tie-break: 1.0 when the
  lower arm is dearer, 0.833 when the upper arm is.
- *All cyclic scenarios.* None of the 175 cells with an arbitrary orientation
  is correct on every cyclic scenario. Mean next-hop accuracy per cell ranges
  from 0.726 to 0.926. The maximum is the value the cost-blind M2 reaches at
  every `γ`, and 111 of the 175 cells sit exactly on it.
- *Against M2 at the same `γ`.* Over 1925 scenario-cells: better in 15, equal
  in 1743, worse in 167. At `γ ≥ 1` accuracy equals M2's in every
  scenario-cell, so phase improves nothing there. At `γ ≤ 0.3` phase does change decisions, mostly for the
  worse.
- *The four cost-sensitive scenarios* (`diamond_upper+1.0`, `lop_lower_cheap`,
  `ring8_heavy`, `ring7_heavy`, where M2 gets a node wrong): fully correct in
  15 of 700 cells. All 15 are `ring7_heavy` at `γ = 0.3`, scattered across
  six orientations and all five values of `β` with no pattern, and none of those
  cells is also correct on `lop_lower_cheap`.
- *Lopsided diamond curve.* A never puts less than 0.646 of its weight on the
  two-hop arm, however cheap the three-hop arm becomes. With the flux
  reversed, a circulating current removes the three-hop arm entirely (split
  exactly 1.0 at every swept cost).

**Oracle orientation (diagnostic, not a candidate).** Giving M3-H the
"away from the destination" orientation does not rescue it: 4 of 100 cells
correct on the cost-sensitive scenarios, and against M2 it is better in 4,
equal in 251, worse in 20. The orientation problem is not the only obstacle.
On the diamond and lopsided scenarios the `id` convention happens to coincide
with `oracle_hops`, because the destination has the highest node id.

**Counterexamples.** On `ring7_heavy` at `γ = 0.1`, a node adjacent to the
destination forwards the long way round, across the heavy edge: a path-cost
stretch of 9. Stretch
exceeds 5 in 88 M3-H cells and in no cell of any other model.

**Damping.** The regime where phase affects decisions is the slow one. In
the dynamic run M3-H settles in 803–905 rounds at `γ = 0.3`, 96–111 at
`γ = 1` and 16–20 at `γ = 3.5`.

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
  cyclic scenario. Mean accuracy ranges from 0.880 to 0.954; four cells are
  above M2's 0.926.
- *Against M2 at the same `γ`.* `θ₀ = -π/2`: better in 8, equal in 103, worse
  in 16. `θ₀ = 0`: better in 6, equal in 94, worse in 17.
- *The four cost-sensitive scenarios:* fully correct in 7 of 91 stable cells,
  all on `diamond_upper+1.0`. Zero on `lop_lower_cheap`, `ring8_heavy` and
  `ring7_heavy`.
- *Diamond curve.* M3-S does respond to cost on equal-hop arms, but the sign
  of the response depends on the rule:

  | `θ₀` | `β` | swept points on the correct side | range of A's upper-arm share |
  |---|---|---|---|
  | `-π/2` | 0.5 | 20 of 20 | 0.495 – 0.506 |
  | `-π/2` | 1 | 20 of 20 | 0.491 – 0.649 |
  | `-π/2` | 2 | 0 of 20 | 0.000 – 0.567 |
  | `0` | 0.5 | 0 of 20 | 0.488 – 0.507 |
  | `0` | 1 | 0 of 20 | 0.386 – 0.506 |
  | `0` | 2 | 20 of 20 | 0.451 – 1.000 |

  Doubling `β` or changing `θ₀` turns a rule that is right at every point
  into one that is wrong at every point. Where it is right at small `β`, the
  margin is under one percentage point.
- *Lopsided diamond curve.* A never puts less than 0.72 on the two-hop arm.
  Hop count dominates cost.

**Interference is present but not informative.** Mean cancellation per node
is 0.21 across stable cyclic cells (0.49 for M2 and M3-H, 0 for M0 and M1).
It does not track accuracy (`plots/cancellation_vs_accuracy.png`).

## 5. Controls

**M0, scalar attenuation.** Stable on every cyclic scenario in 13 of 25
cells. Five are correct on every cyclic scenario: `(κ, γ)` = (2, 1), (2, 2),
(4, 1), (4, 2), (4, 3.5). Its worst usable cell has mean accuracy 0.954,
which equals M3-S's best. Fully correct on the cost-sensitive scenarios in 41
of 61 stable cells. Against M2: better in 52, equal in 112, worse in 7.

M0's weaknesses, for the record:

- It is unstable at low `γ` unless `κ` is large.
- At strong damping it is biased toward fewer hops. On the lopsided curve
  only one of its four stable cells is right at every point.
- Within 0.06 of its stability edge it forms forwarding loops: 5 stable cells,
  4 of them on trees.

**M1, two real states.** Its accuracy differs from M0 in one row of 375.
Under the shared readout the second real state adds nothing. M1 as designed
shows that extra state alone does not help, but it is not a stronger baseline
than M0.

**M2, phase disabled.** Cost-blind by construction. Mean accuracy 0.926 at
every `γ`; correct exactly where hop count and cost agree.

## 6. Predictions

| | Prediction | Outcome |
|---|---|---|
| Pr1 | Every stable model forwards every tree node correctly | **Refuted** for M0 and M1 (59 of 63 stable tree cells) and M3-S (95 of 96). Confirmed for M2 and M3-H. |
| Pr2 | M0 is correct at large `κ`, hop-biased at small `κ` | Confirmed, with the addition that large `γ` also causes hop bias. |
| Pr3 | M3-S on equal-hop arms discriminates weakly, with a sign that depends on `θ₀` and `β` | Sign dependence confirmed. "Weakly" is wrong at `β = 2`, where the split reaches 0 or 1. |
| Pr4 | M3-H on odd cycles is orientation-dependent and no arbitrary orientation is correct across costs | Confirmed. The oracle orientation fails too. |

## 7. The charter's five questions

**1. Which branch has a coherent general cost-to-phase rule and forwarding
interpretation?**

- *M3-H:* the forwarding interpretation is coherent. There is no cost-to-phase
  rule. On bipartite graphs none can exist for symmetric costs, by A4; on the
  odd cycles tested, none of 25 `(β, γ)` cells × 8 orientations orders routes
  by cost.
- *M3-S:* the readout is well defined but is not a current, and the model is
  usable only at strong damping. No rule transfers: the sign of the cost
  response flips with `β` and `θ₀`, and hop count overrides cost on unequal
  arms.

**2. Is any advantage tied to phase or interference rather than state,
damping or a weak baseline?**

No advantage was found for either branch. M3-H's unconditional stability and
freedom from loops are shared by M2, so they come from the Hermitian
structure, not from cost-derived phase. M0 matches or beats both branches on
route quality with half the state.

**3. What did the tree, diamond and ring reveal?**

- *Tree:* M3-H is exactly M2. M0 loops near its stability edge.
- *Diamond:* M3-H is exactly cost-blind. M3-S responds to cost with a
  rule-dependent sign.
- *Lopsided diamond:* neither branch ever prefers a cheaper arm that has more
  hops.
- *Ring:* on the even ring M3-H's accuracy is identical to M2's in every
  cell. On the odd ring phase changes decisions and produces stretch up to 9.

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

See §9.

## 8. Limitations

- **One readout.** Only the power-inflow readout was tested. A readout on the
  phase of each neighbour's contribution was not. For a single dominant path
  that phase is `β × path cost`, so such a readout would be a distance-vector
  protocol carried in phase, with interference as noise rather than signal.
  That is a different hypothesis.
- **One cycle per graph, at most eight nodes.** Graphs with several
  independent loops were not run.
- **Coarse grid.** Five values each of `γ` and `β`. A narrow window between
  grid points would be missed.
- **M1 is not a strong control** (§5). This matters less than planned because
  neither branch beat M0.
- **Settle times are descriptive.** The shared step rule is not tuned per
  model and is too large for M3-S near its stability edge.
- **Exploratory.** Nothing here was preregistered beyond the grid and the
  four predictions in the spec.

## 9. Recommendation

Stop both branches as formulated, and archive Experiment 0 as a negative
result. Three kill criteria from the kickoff are met at smoke-test scale:

- **Kill A, no transfer:** no cell of either branch is correct across the
  cyclic scenarios.
- **Kill B and E, a simpler real method matches it:** M0 is better with one
  real per node.
- **Kill C, destructive interference is mostly pathological** (M3-H at low
  damping: 167 scenario-cells worse than M2 against 15 better).

If one revision is wanted before archiving, the candidate with a specific
rationale is **asymmetric link costs in M3-H**. A phase proportional to
`c_ij − c_ji` is antisymmetric by construction, so it needs no orientation
convention, and chirality is the one thing the Hermitian model has that a
symmetric real field lacks. Two cautions: it only helps on non-bipartite
graphs, and a real field with asymmetric weights is the obvious control and
may match it. I would put the odds of a positive result as low.

I do not recommend a larger held-out test of the current formulations.

The decision is the project owner's.
