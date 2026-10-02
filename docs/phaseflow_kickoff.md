# PhaseFlow — Research Kickoff

> **PhaseFlow explores whether distributed complex-valued interference fields can produce useful adaptive routing behavior on ordinary classical networks.**

**Status:** Phase 0 complete enough to justify Experiment 0  
**Scope:** classical networking only; no quantum hardware required  
**Primary objective:** determine quickly whether phase interference provides a routing signal that matched real-valued local dynamics do not.

## 1. Original intuition

Quantum mechanics suggests structures that may inspire classical distributed systems: superposition, phase, interference, measurement and decoherence. PhaseFlow focuses narrowly on **phase + interference**.

Each node maintains a small complex-valued state:

```text
ψ = r e^(iφ)
```

Nodes exchange state only with neighbors. Multiple route signals combine locally. Because amplitudes carry phase, route families may constructively reinforce or destructively cancel.

The question is not whether this can be simulated. It can.

The question is whether this gives a useful distributed routing primitive.

## 2. Research question

> **Can a fully classical distributed network use locally propagated complex amplitudes and interference to form and adapt forwarding paths more effectively than matched real-valued or probabilistic local routing methods?**

Stronger initial hypothesis:

> **In cyclic / multipath topologies, can interference suppress undesirable route structures and reorganize after topology or cost changes using comparable or lower control-state exchange than classical local baselines?**

## 3. What we are NOT claiming

PhaseFlow is not:
- a quantum network;
- a quantum-computer algorithm;
- a claim of quantum speedup;
- real entanglement or nonlocality;
- faster-than-light communication;
- a replacement for shortest-path routing;
- a claim that complex numbers are inherently better;
- a novelty claim before a deeper literature review.

A complex number is just two real numbers. Any benefit must come from the **interference structure**, not extra state.

## 4. Phase 0 findings

### Quantum walks already exist

Continuous-time quantum walks evolve complex amplitudes across graph nodes. Chiral quantum walks attach complex phases to edges and can strongly alter transport through constructive and destructive interference.

### Phase matters on cycles, not trees

A key result from chiral quantum-walk theory is that on a tree, edge phases can be gauged away. Meaningful phase variables are associated with independent loops.

For a connected graph:

```text
independent loop phases ≈ E - N + 1
```

Therefore:

> **Trees are the negative control.**

If phase appears to provide a systematic advantage on trees, the implementation or comparison is suspect.

### Pure coherent dynamics do not converge naturally

A pure evolution such as:

```text
dψ/dt = -iHψ
```

generally oscillates. Quantum stochastic walks solve this in physics by mixing coherent and dissipative dynamics. PhaseFlow should borrow the principle while remaining a simple classical distributed algorithm.

### Existing routing neighbors

Important classical neighbors/baselines:
- distributed Bellman-Ford;
- scalar diffusion / potential fields;
- probabilistic multipath;
- Q-routing;
- AntNet.

Existing “quantum-inspired routing” also exists, but much of it uses quantum-inspired optimization rather than a continuously distributed interference field as the routing primitive.

Current status:

> **Adjacent work exists; novelty is unproven, but an interference-based classical routing protocol remains worth testing.**

## 5. Candidate model

Use a driven, damped complex field:

```text
dψ/dt = -(γI + iH)ψ + b_d
```

where:
- `ψ` is the complex node-state vector for one destination;
- `H` describes graph connectivity and edge coupling;
- `γ > 0` provides damping/decoherence;
- `b_d` is a persistent source/reference emitted by destination `d`.

At steady state:

```text
ψ* = (γI + iH)^(-1) b_d
```

No node computes that inverse globally. Each node updates only from its neighbors.

A first local discrete implementation may look like:

```text
ψ_i(t+1) =
    (1 - α γ) ψ_i(t)
    + α Σ_j H_ij(t) ψ_j(t)
    + α b_i
```

with a stable discretization chosen during implementation.

## 6. Edge representation

For edge `i <-> j`:

```text
H_ij = a_ij exp(i φ_ij)
```

where:
- `a_ij` controls coupling magnitude;
- `φ_ij` is derived from link cost or condition.

Candidate inputs:
- latency;
- queue delay;
- loss;
- congestion;
- reliability;
- composite normalized link cost.

The central unknown is:

```text
φ_ij = f(c_ij)
```

**Experiment 0 must test whether one general rule `f` transfers across unseen graphs.**

Graph-specific phase tuning does not count as success.

## 7. Forwarding signal

A structurally meaningful edge-flow candidate is:

```text
J_ij = 2 Im(H_ij ψ_i* ψ_j)
```

Interpretation:
- positive current: flow favored from `i` to `j`;
- negative current: opposite direction;
- near zero: little preferred transport.

Candidate forwarding policy:

```text
P(i -> j) ∝ max(J_ij, 0)
```

normalized over neighbors.

A deterministic max-current policy can also be tested.

## 8. Why interference might matter

Classical probabilities only add:

```text
p1 + p2 >= p1
```

Complex amplitudes can cancel:

```text
ψ1 + ψ2 ≈ 0
```

Potentially useful property:

> **Bad cyclic or competing route structures may destructively interfere without explicitly enumerating and scoring every route.**

Potentially fatal property:

> **Good alternate routes may cancel accidentally.**

Experiment 0 should determine which dominates.

## 9. Experiment 0

### Goal

> **Does phase interference contain useful routing information that a matched real-valued local system does not?**

Do not build a packet simulator or polished UI yet.

## 10. Topology suite

### Trees — negative control
- path;
- binary tree;
- random tree.

Expected:

```text
phase-aware ≈ phase-disabled
```

### Diamond

Minimal two-route topology:

```text
      B ---- D
     /        \
S -- A          T
     \        /
      C ---- E
```

Test equal, slightly unequal and changing branch costs.

### Rings / ladders
- simple cycle;
- ladder;
- joined cycles.

### Random cyclic mesh

Generate connected graphs with controlled:

```text
N
mean degree
cycle rank = E - N + 1
```

This is the main generalization test.

## 11. Initial algorithms

### M0 — scalar local field
One real scalar per node/destination.

### M1 — matched two-real-state control
Two real variables per node, same message/state budget as the complex model, but no complex multiplication/interference.

**Critical control.**

### M2 — complex magnitude, phase disabled

```text
φ_ij = 0
```

### M3 — PhaseFlow
Complex states + cost-derived phase + interference.

## 12. Later routing baselines

Only if Experiment 0 survives:
- Dijkstra as a static optimal reference;
- distributed Bellman-Ford;
- ECMP / probabilistic multipath;
- Q-routing;
- AntNet.

The intended regime is:

```text
dynamic cost
dynamic topology
multipath connectivity
local-only information
```

## 13. Perturbations

Each graph should experience the same deterministic schedule:

```text
t0  initialize
t1  settle
t2  increase cost on one branch
t3  restore branch
t4  remove one edge
t5  add alternate edge
t6  change congestion elsewhere
```

Use identical perturbations across models and seeded randomness where needed.

## 14. Metrics

### Stability
- settling time;
- residual oscillation;
- sensitivity to damping.

### Route quality

```text
path cost
path stretch = selected path cost / shortest path cost
```

### Adaptation

```text
T_recover
```

after a link/cost change.

### Loop behavior
- forwarding loops;
- circulation around cycles;
- current trapped in cyclic structures.

### Multipath behavior
- useful path diversity;
- load spread;
- forwarding entropy.

### Control-state cost

```text
state updates sent
bytes/state exchanged
```

Normalize comparisons by equivalent state/message budgets.

### Useful-forwarding correlation

For each node:

```text
does the strongest positive current point toward a low-cost route?
```

Measure across all nodes and graphs.

## 15. Primary generalization test

Generate many unseen cyclic graphs.

For each:

1. assign random positive link costs;
2. derive phases with the same frozen `φ = f(c)`;
3. evolve the field;
4. derive forwarding decisions;
5. compare with shortest-path structure;
6. perturb link costs;
7. measure adaptation.

Do **not** optimize `f` independently per graph.

Tune on one graph family if necessary, then freeze parameters and evaluate on held-out topology families.

## 16. Kill criterion

Stop or fundamentally reframe PhaseFlow if any of the following holds robustly:

### Kill A — no transfer
No single reasonable phase rule `φ = f(c)` produces useful forwarding current across held-out cyclic graphs.

### Kill B — matched real control is equivalent
M1 reaches the same route quality, stability and adaptation at the same state/communication budget.

### Kill C — destructive interference is mostly pathological
Good alternate routes cancel often enough that useful behavior requires graph-specific tuning.

### Kill D — damping destroys the feature
The damping needed for stability removes the measurable interference benefit.

### Kill E — simple multipath matches it
A simple scalar/probabilistic method matches the behavior with lower complexity.

A visually interesting interference pattern is **not success**.

## 17. Continue criterion

Proceed to a real routing simulator only if PhaseFlow demonstrates all of:

1. **Transfer:** one phase rule works on held-out cyclic topology families.
2. **Specificity:** the effect is absent or strongly reduced on trees.
3. **Interference dependence:** M3 materially outperforms M1/M2 on at least one meaningful metric.
4. **Stability:** no persistent pathological oscillation under a usable damping regime.
5. **Mechanism:** results support an interference-related explanation.
6. **Practical signal:** the effect is large enough to justify the extra state/math.

Only then preregister quantitative Experiment 1 criteria.

## 18. Suggested implementation

```text
phaseflow/
  pyproject.toml
  src/phaseflow/
    graph.py
    fields/
      scalar.py
      two_real.py
      complex_no_phase.py
      phaseflow.py
    current.py
    topology.py
    perturbations.py
    metrics.py
    experiment0.py
  tests/
  docs/
    kickoff.md
    literature.md
    experiment-0.md
  results/
```

Python is sufficient.

Suggested libraries:

```text
numpy
scipy
networkx
pandas
matplotlib
pytest
```

No ML frameworks initially.

## 19. First coding milestone

Implement only:

```text
diamond graph
tree control
single ring
scalar field
matched two-real control
phase-disabled complex field
PhaseFlow field
edge current
path extraction
one dynamic edge-cost change
plots of state/current over time
```

Expected output:

```text
results/experiment0-smoke/
  summary.json
  currents.csv
  states.csv
  route_metrics.csv
  plots/
```

Then decide whether the mechanism is coherent before adding random meshes.

## 20. Deferred idea — multi-frequency PhaseFlow

A single phase aliases:

```text
exp(iφ) = exp(i(φ + 2π))
```

A later extension could run:

```text
ψ^(1), ψ^(2), ..., ψ^(m)
```

with:

```text
φ_k = β_k c
```

The combined signature might encode path cost more robustly, analogous to multi-frequency interferometry.

**Do not implement this unless single-field PhaseFlow first demonstrates useful behavior.**

## 21. Starting references

These are nearest-neighbor references, not an exhaustive novelty review.

1. **Enhanced quantum transport in chiral quantum walks**  
   https://link.springer.com/article/10.1007/s11128-024-04331-y

2. **Controlled transport in chiral quantum walks on graphs**  
   https://arxiv.org/abs/2308.12516

3. **Quantum-classical distance as a tool to design optimal chiral quantum walks**  
   https://doi.org/10.1103/PhysRevA.105.032425

4. **Quantum stochastic walks for portfolio optimization** — useful as an example of coherent + dissipative graph dynamics producing a stationary state.  
   https://doi.org/10.1038/s44335-025-00050-4

5. **QSWalk: A Mathematica package for quantum stochastic walks on arbitrary graphs**  
   https://doi.org/10.1016/j.cpc.2017.03.014

6. **Packet Routing in Dynamically Changing Networks: A Reinforcement Learning Approach** — Q-routing baseline.  
   https://papers.nips.cc/paper_files/paper/1993/hash/4ea06fbc83cdd0a06020c35d50e1e89a-Abstract.html

7. **Ant colonies for adaptive routing in packet-switched communications networks** — AntNet baseline.  
   https://doi.org/10.1007/BFB0056909

8. **Hybrid Classical–Quantum Optimization of Wireless Routing Using QAOA and Quantum Walks** — recent neighboring work focused on hybrid optimization.  
   https://arxiv.org/abs/2604.01250

Deeper literature search terms:

```text
complex-valued distributed routing
phase-based routing
wave-based path finding
interference graph algorithms
complex consensus routing
chiral quantum walk classical simulation routing
quantum stochastic walk routing
oscillatory distributed routing
wavefront routing graph
phasor routing
interference-based path selection
```

## 22. Research discipline

Follow the Persiston-style process:

```text
intuition
-> nearest prior work
-> narrow falsifiable claim
-> strongest baseline
-> matched controls
-> seeded experiments
-> preregister success/failure
-> pursue or archive
```

Do not preserve the project because the idea is attractive.

Preserve it only if the mechanism survives controls.

## 23. One-sentence definition

> **PhaseFlow investigates whether locally exchanged complex phases can make interference itself a useful adaptive routing primitive in ordinary classical multipath networks.**

## 24. Tagline

> **Let routes interfere.**
