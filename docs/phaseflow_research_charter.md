# PhaseFlow — Research Charter and Experiment 0 Clarifications

**Status:** Exploratory. Complements `phaseflow-kickoff.md`; it does not replace it.

## 1. Why this project exists

PhaseFlow began as a broader question: can mathematical structures inspired by quantum mechanics (especially superposition, phase, and interference) provide useful behavior in ordinary classical networking, with no quantum computer and no changes to standard routers? After considering several directions, we narrowed the idea to a distributed complex-valued routing field in which path alternatives can reinforce or cancel.

We are not trying to build a networking product immediately. We are looking for a defensible, potentially novel mechanism worth researching further. If the mechanism survives mathematical analysis, prior-art scrutiny, matched controls, and repeatable tests, a future outcome could be publishable research. Novelty and publication are ambitions to investigate, never assumptions or success metrics on their own.

The motivating question is not simply whether complex-valued routing works. It is whether phase-dependent interference yields a causally identifiable, practically meaningful advantage over strong classical alternatives at comparable resource budgets.

## 2. How we work: the Persiston precedent

Our earlier Persiston project followed the chain:

```text
curiosity
-> narrow hypothesis
-> literature review
-> simulator
-> controls
-> repeated seeded experiments
-> pre-registered decision criterion
-> accepted negative result and archived
```

Its central local-regulation hypothesis failed when a simpler periodic re-seed matched performance. Retaining the negative result rather than shifting the goalposts is the working standard for PhaseFlow too.

Rules:

1. **Search nearest prior work before claiming novelty.** Include mathematical analogues and classical routing methods, not just papers labeled "quantum-inspired."
2. **Separate exploratory tests from confirmatory tests.** The smoke test is exploratory and may inform a later hypothesis. Freeze the selected model, controls, metrics, parameters, effect thresholds, and seeds/holdout protocol before confirmatory runs.
3. **Use strong and fair baselines.** A matched two-real-state method is necessary but not sufficient; later compare against appropriate graph diffusion, distributed shortest-path, and adaptive multipath methods at comparable state and message budgets.
4. **Investigate why something works.** Establish mechanism and ablations, not merely attractive plots or a favorable performance number.
5. **Publishability requires a genuinely distinct contribution.** A clear negative result or mathematical limitation is valuable if rigorously demonstrated, but does not by itself prove novelty.
6. **Archive honestly when the criterion fails.** Preserve scripts, seeds, plots, negative findings, limitations, and design decisions.

## 3. Clarification: two exploratory hypotheses, not one interchangeable model

Explore both in the smoke test (choice A), then make a selection or stop.

### M3-H — Hermitian / chiral phase field

- `H_ji = conj(H_ij)`. Driven, damped dynamics: `dψ/dt = -(γI + iH)ψ + b_d` with `γ > 0`.
- Phase effects reside in gauge-invariant loop flux. On a tree the phaseful model is gauge-equivalent to zero-phase for appropriate source transformation; compare gauge-invariant outputs, not raw node phases.
- Need an orientation convention that does not secretly inject a solution or global destination knowledge. Random or arbitrary edge orientation must be treated as part of robustness testing, not optimized per test graph.
- The edge-current convention in the original brief must be corrected: with `J_ij = 2 Im(H_ij ψ_i* ψ_j)`, positive `J_ij` contributes incoming coherent flow from `j` to `i` in `d|ψ_i|²/dt`. The driven, damped model includes injection and dissipation, so current is not globally conserved. Specify packet-forwarding direction separately and test it on a two-node example.

### M3-S — Complex symmetric / accumulated path-phase field

- `H_ji = H_ij = a exp(iβc_ij)` for symmetric link cost `c_ij`. In general this `H` is not Hermitian.
- Walk amplitudes accumulate edge-cost phase, potentially allowing cost-difference interference. Caution: the `-iH` evolution also introduces hop-dependent factors, and the steady state combines all walks, including repeated/cyclic walks; this is not automatically a clean shortest-path estimator.
- The Hermitian continuity/current interpretation does not apply unchanged. Define a mathematically valid, explicitly operational routing readout, and validate it on simple graphs before comparing performance.
- Establish a stable parameter region from the operator spectrum and check transients; arbitrary `γ > 0` is insufficient in general. Do not use Hermitian tree gauge equivalence as this model's negative control.

Neither branch is automatically "PhaseFlow." The smoke test compares coherence and operational usefulness; the promising branch, if any, becomes the main formal hypothesis. Report both outcomes, including unfavorable ones. Avoid combining the strongest result from each variant into a single post-hoc claim.

## 4. Experiment 0 and non-negotiable controls

- Start only with a tree, diamond, and ring; a small deterministic dynamic link-cost perturbation; M0 scalar, M1 carefully matched two-real-state, M2 phase-disabled, and M3-H/M3-S.
- Match state sizes, update opportunities, and communication budgets. Document exactly what information each baseline is allowed to know.
- Calculate exact fixed points as a validation oracle, never as information available to simulated local nodes. Distinguish convergence from routing usefulness.
- Unit-test orientation/sign convention, numerical stability, Hermitian gauge invariance on trees, and the forwarding/readout implementation before any broad sweeps.
- State explicitly what "route" means for each model. A field that is visually interesting but cannot define robust forwarding is not a routing result.
- Keep source code, seeds, configuration, data outputs, baseline details, and experiment decisions reproducible from the first milestone.

## 5. Decision after the smoke test

Create a short `docs/experiment-0.md` recording:

1. Which branch, if either, has a coherent general cost-to-phase rule and forwarding interpretation?
2. Is any advantage genuinely tied to phase/interference rather than more state, damping, or baseline weakness?
3. What did the tree/diamond/ring reveal, including counterexamples?
4. What is the closest prior work and the precise remaining novelty claim, if any?
5. Should we stop, revise the question once with explicit rationale, or preregister a larger held-out graph test?

Only after this decision should we create random-mesh sweeps, a packet simulator, or UI. The eventual experiment should pursue an identifiable new mechanism or useful limitation, not a new branding of an established algorithm.

## 6. Success means learning, not keeping the project alive

PhaseFlow may produce a useful new routing primitive, a rigorous characterization of why interference fails to encode routing costs robustly, or a well-supported negative result that tells us to archive it. All are legitimate scientific outcomes. A convincing result that survives strong controls and prior-art comparison is the only reason to invest in the next phase.
