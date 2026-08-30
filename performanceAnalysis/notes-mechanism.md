# Mechanism map: merger-tree evolution cost structure

Code read: bleeding-edge master @ d6a011b (2026-08-29). File:line references to that tree.

## The evolve loop (mergerTreeEvolverStandard, source/merger_trees/evolver/standard.F90)

- `standardEvolve` (l.291): outer loop repeats **full tree walks** until a walk evolves
  nothing. Each walk constructs a `mergerTreeWalkerAllNodes` and visits **every node
  in the tree** (post-order DFS: satellites -> children -> parent;
  walkers/all_nodes/_class.F90 l.101).
- Per visited node: `nodeIsEvolvable` O(1) check (l.723). If evolvable:
  `timeEvolveTo` (l.762) computes the max time this node may reach, then
  `mergerTreeNodeEvolver_%evolve` advances it (ODE). After a successful advance the
  evolver *re-evolves the same node* (`nodeNext => node`, l.517) until blocked; with
  `backtrackToSatellites=true` it also revisits the node's satellites (l.514).
- Node promotion continues up the branch within the same walk (l.501), so purely
  branch-internal evolution does NOT multiply walk count.

## Per-call time-limit computation (standardTimeEvolveTo, l.762)

Scans performed *each time a node is considered*:
1. Tree events list, node events list (usually short).
2. `mergerTreeEvolveTimestep_%timeEvolveTo` — for the reference DMO model a `multi`
   of: `simple` (0.1/H cap; computes expansionFactor+expansionRate each call),
   `satellite` (merger-time sync), `satelliteDestruction` (destruction-time cap).
3. Satellite-in-host limit: satellite may reach host time + min(0.1/H, 1 Gyr)
   [timestepHostRelative/Absolute]; `fractionTimestepSatelliteMinimum` (0.75 in the
   reference model, ***default 0***) defers a satellite that can only take a small
   fraction of that.
4. **Sibling scan**: O(#siblings) for primary progenitors (l.946-959).
5. **Hosted-satellite scan**: O(#satellites of this node) (l.963-984). For the
   final host halo this is O(S_total) *every time the host is considered*.
6. **Mergee scan**: O(#mergees registered on this node) (l.990-1004). Note stray
   comment at l.989 ("need to also figure out how to speed up satellite evolve
   checks").

## Merger synchronization (timesteps/satellite.F90)

A satellite about to merge must wait until its merge target has been evolved to
within `timeOffsetMaximumAbsolute` (=0.01 Gyr) [or relative 0.001*t] of the merge
time; conversely the host's step is capped at each mergee's merge time (evolver
l.990). **Every merger/destruction event therefore forces the host+system to land
steps within a ~10 Myr window** around the event, fragmenting the host's ~0.1/H
(~ few hundred Myr at z~0) natural stepping.

Satellite destruction (timesteps/satellite_destruction.F90): satellites record a
destruction time (mass below threshold; reference model ties threshold to
m_res); the timestep is capped to exactly that time, and the destruction is
processed as an end-of-step task.

## Per-ODE-call fixed overhead (node_evolver/standard.F90 standardEvolve, l.385)

Every `evolve` call, regardless of step size: 2x OMP tree-lock cycles,
differentialEvolutionPre/Post over all node operators, Calculations_Reset,
serializationOffsets x3, serializeCount x2, serializeValues x2 + saved copies,
scale initialization + nodeOperator scales pass + serializeScales, **fresh
odeSolver construction** (GSL stepper allocation), initial-step-size estimation
(with `reuseODEStepSize=false`, the reference setting), deserialize, analytics
solve. This is a large constant per call; calls with tiny dt cost nearly as much
as calls with large dt.

## Cost model hypothesis

Let S = peak/typical live subhalo count in the final host (∝ M/m_res at fixed
physics), N = built tree node count (∝ 1/m_res, possibly slightly steeper), W =
number of full tree walks, E = number of merger+destruction events (∝ S).

- Component (a): per-walk overhead = O(live nodes) evolvability checks + O(S)
  host satellite scan + O(S) ODE calls (satellites each advance or get checked).
  Total ≈ W x S.
- Component (b): W itself. With backtracking, branch evolution doesn't drive W;
  the host/satellite ping-pong and **event synchronization** do. If each of the E
  events costs an extra walk (or a fixed number of them), W ≈ W0 + c*E ∝ S at high
  resolution.
- Total ≈ (W0 + cS) x (aS + b) → linear in S at low S, **quadratic at high S**;
  measured slope passes through ~1.5 in the transition.
- Component (c): tree build + branch-segment evolution ∝ N (linear), plus any
  superlinearity in N(m_res) itself (tree construction temporal refinement).

## Instrumentation available

- `metaData/treeTiming` (mergerTreeOperator=treeProcessingTimer): per-tree
  timeConstruct, timeEvolve, countNodes (built-tree node count), treeMass.
- `metaData/evolverProfiler` (mergerTreeEvolveProfiler=simple +
  mergerTreeNodeEvolver profileOdeEvolver=true): histograms vs step size of step
  count, RHS evaluation count, CPU time; per-property limiter hit counts.
  NOTE: adds per-step error-analysis cost — do not mix with timing runs.
- verbosityLevel=warn prints "Evolving tree [N]" walk-counter lines (throttled,
  ~10% accuracy on final count) and per-tree run-time summary lines.

## Empirical results (batch 1, reference DMO subhalo model, M_host=1e13, 4 trees/point)

| m_res | nodes/tree | tCon [s] | tEvo [s] | walks (max) |
|-------|-----------|----------|----------|-------------|
| 1e10  |     1604  |  0.41    |   0.18   |    8  |
| 3.16e9|     4890  |  0.39    |   0.57   |   10  |
| 1e9   |    14484  |  0.63    |   1.79   |   16  |
| 3.16e8|    44122  |  1.69    |   5.54   |   52  |
| 1e8   |   132989  |  4.15    |  17.10   |  195  |

- Node count: slope 0.94-0.97 per decade -> N ∝ (1/m_res)^0.95. The tree builder
  does NOT superlinearly refine branches in time (hypothesis B ruled out at these
  resolutions; nodes/branch stays ~constant).
- tEvolve slope: 0.98-0.99 -> still LINEAR to M/m_res = 1e5. The break must lie
  deeper (batch 2 probes 3.16e7, 1e7).
- Walk count: 8 -> 195, slope rising to ~1.15 in last decade (superlinear).
- Per-walk reports (verbosity=warn) show two regimes: healthy walks evolving
  ~80% of live nodes, and stalled walks evolving <3% (e.g. 1262 live, 21
  evolved; sequences of ~15+ such walks) - the merger-synchronization cascades.
- Total node visits (walks x live nodes): ~2.4k/tree at 1e9 -> ~245k/tree at 1e8:
  scaling ∝ (1/m_res)^2.0. At 1e8 visits are still only ~5% of evolve cost
  (per-visit cost ~ 1-5 us vs per-evolution cost ~140 us) -> the quadratic term
  is present but subdominant; extrapolation puts its takeover 1-2 decades deeper,
  consistent with the user's observed transition to ~(1/m_res)^1.5.

## Strategy sketch: event-driven (priority queue) evolver

The `timeEvolveTo` machinery already identifies, for every blocked node, WHICH
node blocks it (`nodeLock`, currently used only for deadlock reporting). An
event-driven evolver inverts this: maintain per-node waiter lists; a queue of
evolvable nodes ordered by time; popping a node evolves it to its limit,
re-registers it against its new blocker, and wakes its waiters (recompute
their limits, reinsert). Tree mutations (merger/destruction/promotion) wake
affected neighbors. Complexity O((evolutions + wakes) log N) replacing
O(walks x N_live); wakes ~ O(evolutions) since each advance wakes only direct
waiters. This eliminates the stalled-walk regime entirely (walks with 1200
visits per ~20 evolutions).

## Live gdb samples during a 1e7 tree (first-walk evolution phase)

Hot stack family (5/6 samples): satelliteTidalHeating rate evaluation ->
massDistributionSphericalHeated%radiusInitial (Brent root-find per call) ->
tidalSpecificEnergy(Terms) -> NFW velocityDispersion1D (log/fma) and
darkMatterOnlyGet (memmove of mass distribution objects). I.e., the dominant
LINEAR-term payload is heated-profile structure solves inside satellite ODE
RHS evaluations. Candidates for constant-factor reduction: memoize/interpolate
heated-profile radius solutions; avoid per-call mass-distribution object
copies (memmove in darkMatterOnlyGet).

## Variant: fractionTimestepSatelliteMinimum=0 (the CODE DEFAULT) at m_res=1e8

tEvolve = 142.1 s/tree vs 17.1 s/tree baseline -> 8.3x slower, with IDENTICAL
trees, node counts, and walk counts (195). The whole penalty is satellites
taking a micro ODE call on every walk visit (host+eps ratchet) instead of
deferring: each call pays the full fixed overhead (serialization, fresh GSL
solver, initial-step-size search) for a ~10 Myr-or-less step. The deferral
mechanism (reference value 0.75) is the difference between linear scaling to
M/m=3e5+ and an early, steep superlinear regime. Production configs using the
default 0.0 would see the break decades earlier - very likely the user's
observed 1/m^1.5.
