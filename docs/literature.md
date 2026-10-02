# PhaseFlow — prior-art note

**Search date:** 2026-10-02. **Status:** a focused search, not an exhaustive review.

How to read this note:

- Every entry under "Verified references" was retrieved at the URL given, and
  its title, authors and year were confirmed there. Caveats are stated per
  entry where the page was a preprint, an author copy, or had no abstract.
- Six entries were re-fetched independently as a spot check and matched:
  Zimborás et al. (arXiv:1208.4049), Matsuoka et al. (arXiv:2411.12191),
  Higuchi et al. (arXiv:2002.05261), Hatano et al. (arXiv:2409.01873),
  Rosales (arXiv:2605.27425), Sanfui and Dutta (arXiv:2609.39890). The rest
  rest on the single retrieval made during the search.
- "Unverified leads" were recalled or seen only in snippets. Do not cite them
  without checking.
- Absence here means only that these searches did not find something.

The two formulations are those of
`docs/superpowers/specs/2026-10-02-experiment0-smoke-design.md`: **M3-H**
(Hermitian, chiral, driven and damped) and **M3-S** (complex symmetric, walk
phase `β × cost`). The **real control** is the damped attenuation field whose
steady state sums walks weighted by `exp(-κ·cost)`.

## 1. Verified references

### 1.0 Chiral walks on trees and bipartite graphs

Both papers confirm the structure behind the M3-H controls. On trees the
phases have no effect at all. On bipartite graphs the phases do change
transition probabilities, but symmetrically under phase reversal, so no
directional bias is possible.

**Zimborás, Faccin, Kádár, Whitfield, Lanyon, Biamonte (2013), "Quantum
Transport Enhancement by Time-Reversal Symmetry Breaking", Sci. Rep. 3, 2361.**
https://arxiv.org/abs/1208.4049

- Introduction: "time-reversal asymmetry cannot affect the site-to-site
  transport in some simple cases, such as linear chains and trees—this is
  proven in the Methods Section."
- Discussion: "we prove that the site-to-site transfer probability is
  insensitive to phases in tree graphs. For bipartite graphs the phases can
  have an effect, however, the dynamics still remains time-reversal symmetric."
- Methods: "(i) phases on tree graphs can be transformed out completely and
  (ii) the sum of phases along loops is invariant under gauge
  transformations." The proof covers the Lindblad open-system case.
- Supplementary: "Directing transport in bipartite graphs is forbidden", by
  the sublattice sign flip; also "Even chiral loops suppress transport as
  θ = π results in p = 0 for all times".
- Relation to PhaseFlow: these statements concern site-to-site transfer
  probability and sink occupancy, not steady-state current under a constant
  source. The M3-H bipartite result (spec A4) uses the same sublattice
  transformation, so it is a corollary of known structure, not a new symmetry.

**Lu, Biamonte, Li, Li, Johnson, Bergholm, Faccin, Zimborás, Laflamme, Baugh,
Lloyd (2016), "Chiral Quantum Walks", Phys. Rev. A 93, 042302.**
https://arxiv.org/abs/1405.6209

- Table I classifies linear chains, trees and bipartite graphs ("only even
  cycles") as probability time-symmetric for all phases; non-bipartite graphs
  are not.
- Section 2 gives the gauge form `H → Λ†HΛ = −H`, "which immediately implies
  the probability time-symmetry condition".
- Caveat that matters here: "Note that this disallows self-edges". The
  bipartite result needs uniform diagonal terms. Node-dependent damping on a
  bipartite graph with cycles would fall outside it.

### 1.1 Classical routing or path-finding using complex amplitudes or interference

No classical network routing protocol was found whose forwarding rule comes
from a distributed complex-valued interference field. The nearest items:

**Matsuoka, Ohno, Segawa (2024), "A mathematical framework for maze solving
using quantum walks".** https://arxiv.org/abs/2411.12191 (search results list
it as Quantum Inf. Process. 24, 93 (2025); not retrieved directly).

- Finds shortest paths from the stationary state of a non-unitary Grover walk
  with absorbing holes and self-loops, on unweighted graphs.
- Framed as classical: "this approach enables efficient computations using
  GPUs ... suggesting its potential as a classical algorithm inspired by
  quantum mechanics."
- On graphs with alternative routes it reports weaker selectivity: on ladder
  graphs "the amplitude of the navigation vector may not be zero even at a
  point outside of the shortest route", and for a sink-at-goal variant with a
  cycle "the position of the quantum walker does not tells us the shortest
  route. Probably, it misleads us."
- Differences: discrete-time coined walk, unweighted edges, no cost phases,
  no chiral element, centralised computation, readout by amplitude
  localisation rather than current. This is the closest item overall.

**Matsuoka, Yuki, Lavička, Segawa (2021), "Maze Solving by a Quantum Walk
with Sinks and Self-Loops: Numerical Analysis", Symmetry 13, 2263.** Title,
authors and year confirmed at https://github.com/kyuki-rp/qw-maze-solving;
the publisher page returned 403, so the description rests on the 2024 paper,
which notes the method "is not suitable for mazes with cycles that have an
odd number of nodes".

**Howard, Dasore, Dhungana, Kuttala, Murphy, Soo, Haque (2026), "Hybrid
Classical–Quantum Optimization of Wireless Routing Using QAOA and Quantum
Walks".** https://arxiv.org/abs/2604.01250

- States the framing conceptually: "quantum walks treat the network itself as
  the substrate of evolution", and "desirable links may contribute
  differently to phase evolution than undesirable ones".
- A discussion paper aimed at quantum hardware subroutines. No distributed
  classical protocol, no current-based forwarding rule, no numerical results.

**Quantum graphs (physics counterpart of M3-S).** Lawrie, Gnutzmann, Tanner
(2023), "Closed form expressions for the Green's function of a quantum graph
– a scattering approach", https://arxiv.org/abs/2309.11251, building on Barra,
Gaspard (2001), "Transport and dynamics on open quantum graphs", Phys. Rev. E
65, 016205, https://arxiv.org/abs/cond-mat/0107285. A wave accumulates phase
`k × edge length` and the steady-state Green's function is a sum over
trajectories. Only the abstracts were confirmed. Not a routing method.

**Complex-weighted network theory (background for M3-H).**

- Tian, Lambiotte (2023), "Structural Balance and Random Walks on Complex
  Networks with Complex Weights", https://arxiv.org/abs/2307.01813.
- Böttcher, Porter (2024), "Complex networks with complex weights", Phys.
  Rev. E 109, 024314, https://arxiv.org/abs/2212.06257.
- Dong, Qiu (2014), "Complex Laplacians and Applications in Multi-Agent
  Systems", https://arxiv.org/abs/1406.1862. Complex-valued distributed
  consensus for formation control, not path selection.
- Fanuel, Alaíz, Fernández, Suykens (2018), "Magnetic Eigenmaps for the
  Visualization of Directed Networks", ACHA 44, 189–199,
  https://arxiv.org/abs/1606.08266.

**Keyword matches that are not this mechanism.**

- Kirst, Timme, Battaglia (2016), "Dynamic information routing in complex
  networks", Nat. Commun. 7, 11061, https://arxiv.org/abs/1510.05033.
  Oscillator phase states and transfer entropy, not cost-based paths.
- Catanuto, Toumpis, Morabito (2007), "Opti{c,m}al: Optical/Optimal Routing in
  Massively Dense Wireless Networks", https://arxiv.org/abs/cs/0608012.
  Geometrical optics, no interference.
- Adamatzky (2016), "Physical maze solvers. All twelve prototypes implement
  1961 Lee algorithm", https://arxiv.org/abs/1601.04672.

### 1.2 Driven or sourced quantum walks for transport, routing or search

**Higuchi, Segawa (2019), "Dynamical system induced by quantum walk", J.
Phys. A 52.** https://arxiv.org/abs/1812.04730. Grover walk with a constant
external inflow: "for any connected internal graph, a stationary state exists,
moreover a perfect transmission to the opposite tail always occurs in the long
time limit." Prior art for a constantly sourced walk with a steady-state
readout; discrete-time, no cost phases.

**Higuchi, Sabri, Segawa (2020), "Electric circuit induced by quantum walk",
J. Stat. Phys. 181, 603–617.** https://arxiv.org/abs/2002.05261. With constant
input, the stationary state combines a reversible measure with the electric
current of the internal graph. The steady state of a sourced walk is tied
here to the classical electrical-flow solution, which is the real-valued
competitor. Related: Higuchi, Sabri, Segawa (2023), Quantum Stud. Math.
Found. 10, 307–316,
https://link.springer.com/article/10.1007/s40509-023-00296-2; Konno, Segawa,
Štefaňák (2021), https://arxiv.org/abs/2105.03111.

**Chiral-walk routing (unitary quantum state transfer).**

- Bottarelli, Frigerio, Paris (2023), "Quantum routing of information using
  chiral quantum walks", AVS Quantum Sci. 5, 025001,
  https://arxiv.org/abs/2302.08472.
- Ragazzi, Cavazzoni, Benedetti, Bordone, Paris (2025), "Scalable Structure
  For Chiral Quantum Routing", Entropy 27, 498,
  https://arxiv.org/abs/2502.12739.
- Yu, Cai (2023), "Controlled transport in chiral quantum walks on graphs",
  https://arxiv.org/abs/2308.12516.
- Annoni, Frigerio, Paris (2024), "Enhanced quantum transport in chiral
  quantum walks", Quantum Inf. Process. 23, 117,
  https://arxiv.org/abs/2308.14747.
- Sanfui, Dutta (2026), "Quantum State Routing and Perfect State Transfer on
  Signed Graphs under Environmental Noise", https://arxiv.org/abs/2609.39890.

All use the M3-H Hamiltonian, but the phase is a designed control knob on
small engineered graphs and the readout is transfer fidelity at a chosen
time. None derives the phase from link cost, uses a steady state, or is a
distributed classical protocol.

**Classical-linear realisations of the M3-H equation.**

- Ranzani, Aumentado (2015), "Graph-based analysis of nonreciprocity in
  coupled-mode systems", New J. Phys. 17, 023024,
  https://arxiv.org/abs/1406.4922: "nonreciprocity arises as a consequence of
  multi-path interference and dissipation in key ancillary modes."
- Metelmann, Türeci (2018), "Nonreciprocal Signal Routing in an Active Quantum
  Network", Phys. Rev. A 97, 043833, https://arxiv.org/abs/1703.04052.
- Koch, Houck, Le Hur, Girvin (2010), "Time-reversal-symmetry breaking in
  circuit-QED-based photon lattices", Phys. Rev. A 82, 043811,
  https://www.osti.gov/etdeweb/biblio/21454822.

Driven, damped, linear coupled-mode networks with complex couplings are the
same mathematics as M3-H. This literature engineers the phases by hand to
route signals in hardware; it does not let link costs set them.

**Other open-system transport on networks.**

- Caruso, Crespi, Ciriolo, Sciarrino, Osellame (2015), "Fast Escape from
  Quantum Mazes in Integrated Photonics", https://arxiv.org/abs/1501.06438.
  Maze exit is reached best "by partially suppressing the presence of
  interference".
- Hatano, Katsura, Kawabata (2026), "Quantum transport on Bethe lattices with
  non-Hermitian sources and a drain", J. Phys. A 59, 175002,
  https://arxiv.org/abs/2409.01873.
- Sánchez-Burillo, Duch, Gómez-Gardeñes, Zueco (2012), "Quantum Navigation
  and Ranking in Complex Networks", Sci. Rep. 2, 605,
  https://arxiv.org/abs/1202.3471.
- Liu, Balu (2017), "Steady states of continuous-time open quantum walks",
  Quantum Inf. Process. 16, 173, https://arxiv.org/abs/1604.05652.
- Razzoli, Paris, Bordone (2021), "Transport efficiency of continuous-time
  quantum walks on graphs", Entropy 23, 85, https://arxiv.org/abs/2011.13794.
- Reviews: Mülken, Blumen (2011), Phys. Rep. 502, 37–87,
  https://arxiv.org/abs/1101.2572; Manzano, Hurtado (2018), Adv. Phys. 67,
  1–67, https://arxiv.org/abs/1707.07895; Rebentrost et al. (2009),
  "Environment-Assisted Quantum Transport", New J. Phys. 11, 033003,
  https://arxiv.org/abs/0807.0929.

**Quantum path-finding via electrical flow (needs quantum hardware).**
Jeffery, Kimmel, Piedrafita (2023), "Quantum Algorithm for Path-Edge
Sampling", https://arxiv.org/abs/2303.03319; Wesołowski, Piddock (2024),
"Advances in quantum algorithms for the shortest path problem",
https://arxiv.org/abs/2408.10427.

### 1.3 Real-valued analogues: the strongest classical competitors

The real control is covered by existing work in at least four independent
literatures.

- **Akamatsu (1996), "Cyclic flows, Markov process and stochastic traffic
  assignment", Transp. Res. B 30(5), 369–386.**
  https://ideas.repec.org/a/eee/transb/v30y1996i5p369-386.html. Logit
  assignment over all paths, cyclic ones included, computed "by some matrix
  operations". This is the `exp(-θ·cost)` walk sum via a matrix inverse.
- **Randomized shortest paths.** Saerens, Achbany, Fouss, Yen (2009), Neural
  Comput. 21(8), 2363–2404 (record confirmed via PubMed; abstract truncated).
  Kivimäki, Shimbo, Saerens (2014), Physica A,
  https://arxiv.org/abs/1212.1666: a free-energy distance interpolating
  between shortest-path and resistance distance. Leleux, Courtain, Guex,
  Saerens (2021), https://arxiv.org/abs/2007.00419: routing policies that
  converge to the least-cost DAG as temperature goes to zero. M3-S is formally
  this construction at imaginary inverse temperature; no such continuation
  appears in the abstracts retrieved.
- **Todorov (2006), "Linearly-solvable Markov decision problems", NIPS 19.**
  https://papers.neurips.cc/paper/3002-linearly-solvable-markov-decision-problems.
  An exponential transform makes the Bellman equation linear; shortest paths
  are a limiting case. This is soft-min as a linear, neighbour-local
  recursion.
- **Chebotarev, Deza (2013), "A topological interpretation of the walk
  distances".** https://arxiv.org/abs/1111.0284. Distances from logarithms of
  entries of `(I − tA)⁻¹`.
- **Zero-temperature and Varadhan limits.** Litvinov (2012),
  https://arxiv.org/abs/1203.0522: "the Hamilton-Jacobi-Bellman equation is
  treated as a result of the Maslov dequantization applied to the Schrödinger
  equation." Steinerberger (2018), https://arxiv.org/abs/1801.02183. Crane,
  Weischedel, Wardetzky (2013), "Geodesics in Heat",
  https://arxiv.org/abs/1204.6216.
- **Potential, field and electrical routing.** Basu, Lin, Ramanathan (2003),
  "Routing Using Potentials: A Dynamic Traffic-Aware Routing Algorithm",
  SIGCOMM, https://conferences.sigcomm.org/sigcomm/2003/papers.html. Lenders,
  May, Plattner, "Density-based Anycast",
  https://lenders.ch/publications/journals/ton08.pdf (author copy; venue not
  confirmed). Kelner, Maymounkov (2009), "Electric routing and concurrent flow
  cutting", https://arxiv.org/abs/0909.2859. Banirazi, Jonckheere,
  Krishnamachari (2019), "Heat-Diffusion", https://arxiv.org/abs/1902.05649.
  Steinerberger (2020), "A Spectral Approach to the Shortest Path Problem",
  https://arxiv.org/abs/2004.01163.
- **Physarum.** Bonifaci, Mehlhorn, Varma (2012), "Physarum Can Compute
  Shortest Paths", https://arxiv.org/abs/1106.0423. Tero et al. (2010),
  Science 327, 439–442,
  https://ora.ox.ac.uk/objects/uuid:01616eb5-3b21-4848-8b63-f909f34a83cc. Li,
  Torres, Thomas, Rossi, Shen (2011), Swarm Intell. 5, 183–223,
  https://link.springer.com/article/10.1007/s11721-011-0063-y. Sun (2017),
  https://arxiv.org/abs/1712.02910.
- **Walk-sum background.** Estrada, Hatano (2008), "Communicability in complex
  networks", https://arxiv.org/abs/0707.0756. Malioutov, Johnson, Willsky
  (2006), "Walk-Sums and Belief Propagation in Gaussian Graphical Models",
  JMLR 7, 2031–2064, https://www.jmlr.org/papers/v7/malioutov06a.html. Newman
  (2003), https://arxiv.org/abs/cond-mat/0309045.

### 1.4 "Quantum-inspired routing": optimiser or interference field?

Every example retrieved is a quantum or quantum-inspired optimiser over
candidate routes. None holds a distributed interference field at the nodes.

- Uthayakumar, Jayaraman, Raja, Shabbir (2025), "QSEER", Sensors,
  https://pmc.ncbi.nlm.nih.gov/articles/PMC12473431/.
- Chen, Date, Chancellor, Atiquzzaman, Sreenan (2022), IEEE Trans. Quantum
  Eng. 3, https://arxiv.org/abs/2110.06321. Centralised quantum annealing.
- Rosales (2026), https://arxiv.org/abs/2605.27425. Classical annealing and
  tensor networks over an effective routing Hamiltonian.
- Bennett, Matwiejew, Marsh, Wang (2021), "Quantum walk-based vehicle routing
  optimisation", https://arxiv.org/abs/2109.14907. The walk is over the
  solution space, not the network graph.

### 1.5 Can interference order two routes by cost? Dark nodes

**Sign blindness of two-terminal interference.** Bedkihal, Segal (2014),
Phys. Rev. B 90, 235411, https://arxiv.org/abs/1407.5690: "The linear
conductance of a two-terminal Aharonov-Bohm interferometer is an even function
of the applied magnetic flux, as dictated by the Onsager-Casimir symmetry."
This concerns total conductance, not individual bond currents.

**Bond-current direction in a ring with unequal arms is not a reliable
arm-length indicator.**

- Rai, Hod, Nitzan (2010), "Circular Currents in Molecular Wires", J. Phys.
  Chem. C 114, 20583–20594, https://arxiv.org/abs/1006.1729: "the direction of
  the circular current can depend on the imposed bias"; circular currents can
  exceed the net current.
- Yadalam, Harbola (2016), Phys. Rev. B 94, 115424,
  https://arxiv.org/abs/1606.00638: "current through one of the branches can
  be completely suppressed".
- Benjamin, Jayannavar (2002), Int. J. Mod. Phys. B 16, 1787,
  https://arxiv.org/abs/cond-mat/0112186. Current magnification: one arm's
  current can exceed the transport current while the other flows against it.

**Dark nodes and trapping.** Sett, Pan, Falloon, Wang (2019), "Zero Transfer
in Continuous Time Quantum Walks", https://arxiv.org/abs/1902.11115. Vidal,
Mosseri, Douçot (1998), "Aharonov-Bohm cages", Phys. Rev. Lett. 81, 5888,
https://arxiv.org/abs/cond-mat/9806068. Le, Donati, Severini, Caruso (2018),
"How to Suppress Dark States in Quantum Networks and Bio-Engineered
Structures", https://arxiv.org/abs/1707.07482. Caruso, Chin, Datta, Huelga,
Plenio (2009), J. Chem. Phys. 131, 105106, https://arxiv.org/abs/0901.4454.
Markussen, Stadler, Thygesen (2010), Nano Lett. 10, 4260–4265,
https://orbit.dtu.dk/en/publications/the-relation-between-structure-and-quantum-interference-in-single.

No paper was found that poses the exact question "can steady-state
interference of two routes reveal which is cheaper?" The elementary point
that `|e^{iβc₁} + e^{iβc₂}|²` depends only on `cos(β(c₁ − c₂))`, and is
therefore even in the cost difference, is our own remark, not a citation.

## 2. Unverified leads

- Koch, Hillery, "Finding paths in tree graphs with a quantum walk", Phys.
  Rev. A 97, 012308 (2018), arXiv:1710.05084.
- Reitzner, Hillery, Koch, "Finding paths with quantum walks or quantum
  walking through a maze", Phys. Rev. A 96, 032323 (2017).
- Hillery, "Finding more than one path through a simple maze with a quantum
  walk", J. Phys. A 54, 095301 (2021).
- Chebotarev, "The walk distances in graphs", Discrete Appl. Math. (2012).
- Françoisse, Kivimäki, Mantrach, Rossi, Saerens, "A bag-of-paths framework
  for network data analysis", Neural Networks 90 (2017).
- Kalantari, Shayman, "Routing in wireless ad hoc networks by analogy to
  electrostatic theory", IEEE ICC 2004.
- Baumann, Heimlicher, Lenders, May, "HEAT: Scalable Routing in Wireless Mesh
  Networks Using Temperature Fields" (2007).
- Toumpis, "Mother nature knows best", Computer Networks 52 (2008).
- Steinbock, Tóth, Showalter, "Navigating complex labyrinths: optimal paths
  from chemical waves", Science 267, 868 (1995).
- Kottos, Smilansky, quantum-graph scattering papers (1997–1999).
- Büttiker (1986), four-terminal phase-coherent conductance and reciprocity.
- Jayannavar, Deo, current-magnification papers (mid-1990s, Phys. Rev. B).
- Dial (1971) and Bell (1995), logit traffic assignment.
- Tassiulas, Ephremides (1992), backpressure routing; Intanagonwiwat et al.
  (2000), directed diffusion; Katz (1953), centrality.

## 3. Closest prior work

1. **Matsuoka, Ohno, Segawa (2024), with the 2021 predecessor.** Stationary
   state of a sourced, absorbing quantum walk read as a shortest-path finder,
   proposed as a quantum-inspired classical algorithm; weaker once cycles
   offer alternative routes.
2. **Zimborás et al. (2013) and Lu et al. (2016).** The M3-H Hamiltonian, with
   the tree and bipartite gauge results that the smoke test reproduces in the
   driven setting.
3. **Higuchi, Segawa (2019) and Higuchi, Sabri, Segawa (2020).** A constantly
   driven walk whose stationary state is the classical electrical current.
4. **Akamatsu (1996), Saerens et al. (2009), Kivimäki et al. (2014), Todorov
   (2006).** The real control is this literature; M3-S is formally its
   imaginary-parameter continuation.
5. **Ranzani and Aumentado (2015), Metelmann and Türeci (2018).**
   Driven-damped coupled-mode networks with complex couplings: the M3-H
   equation in classical hardware, with phases designed by hand.

## 4. What these searches did not find

- A distributed classical routing protocol in which nodes hold
  per-destination complex amplitudes, exchange them with neighbours, and
  forward along probability current.
- Chiral edge phases derived from link cost. Every chiral-walk routing paper
  retrieved treats the phase as a designed control parameter.
- M3-S as such: the all-walks resolvent at imaginary `θ`, used for routing.
  Quantum graphs supply the physics but not the routing use.
- The bipartite phase-reversal invariance stated for steady-state currents
  under a constant source. It follows from the Zimborás and Lu gauge argument,
  so it is weak as a standalone novelty claim.
- A paper asking whether steady-state interference of two routes can order
  them by cost. Adjacent results (two-terminal flux evenness, bias-dependent
  bond-current direction, the Matsuoka ladder result) all point toward "not
  reliably".

Not uncovered: the real control, potential and electrical routing, Physarum
routing, soft-min as a linear recursion, and "quantum-inspired routing" as
optimisers. All have direct prior work.

## 5. Search terms used

"wave-based routing", "phase-based routing", "phasor routing",
"interference-based routing / path selection", "complex-valued distance
vector", "quantum-inspired distance vector / Bellman-Ford", "probability
current routing protocol", "quantum walk routing wireless sensor network",
"chiral quantum walk routing", "edge weight encoded as phase", "path integral
/ sum over paths routing", "Wick rotation randomized shortest paths", "complex
temperature randomized shortest paths", "stationary phase shortest path graph
complex edge weights", "Helmholtz shortest path graph", "quantum stochastic
walk steady-state current source sink", "non-Hermitian continuous-time quantum
walk source drain", "Aharonov-Bohm ring unequal arm lengths circulating
current direction", "bond current flux reversal bipartite", "dark states
quantum networks".

Pages that could not be retrieved (403, captcha or empty): the PMC and Nature
pages for Zimborás (arXiv full text used instead), ScienceDirect (Akamatsu,
Toumpis, Mülken), ACM DL (Basu), IEEE Xplore (Kalantari, HEAT), MDPI (Matsuoka
2021), MIT Press and Semantic Scholar (Saerens), PubMed (Steinbock).
