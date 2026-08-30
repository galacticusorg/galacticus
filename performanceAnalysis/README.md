# Run-time scaling of merger-tree evolution with mass resolution

Analysis of why CPU time per merger tree scales as ~1/m_res at low resolution
but steepens toward ~1/m_res^1.5 at high resolution, and what can be done about
it. All measurements use the prebuilt `bleeding-edge` binary (master @ d6a011b),
the reference dark-matter-only subhalo model
(`parameters/reference/darkMatterOnlySubHalos.xml` physics: orbiting satellites,
Chandrasekhar dynamical friction, Zentner2005 tidal stripping, Gnedin1999 tidal
heating, merging at 0.01 r_vir, destruction below m_res), an Eisenstein & Hu
transfer function, fixed random seed, single OpenMP thread, and (unless noted) a
10^13 M_sun host with 4 trees per point. Scripts to reproduce everything are in
`scripts/`; per-run records in `data/results.jsonl`.

## 1. Summary of findings

1. **In the reference configuration, evolution time is linear in 1/m_res far
   deeper than expected** - slope 0.98-1.01 per decade all the way to
   M/m_res = 3x10^5, at a constant ~137 us per built tree node. The break
   begins at M/m_res ~ 10^6 (slope 1.10 and rising).

2. **The superlinear machinery is present and measurable long before it
   dominates.** The evolver advances trees by repeated full-tree walks; the
   number of walks grows superlinearly (8 walks at M/m=10^3 up to 3,659 at
   M/m=10^6; local slope ~1.5), and total node visits (walks x live nodes)
   grow as (1/m_res)^2.06. These visits are individually cheap (~us), which is
   why they only overtake the linear ODE work (~130 us/step) at very high
   resolution - but they inevitably do, and the measured slope passes through
   ~1.5 during the crossover.

3. **The proximate cause of the extra walks is event synchronization, not
   branch time-resolution.** Walk logs show "stalled" walks that visit ~1,300
   live nodes to evolve ~20 (sometimes 0): merger events demand the merge
   target be within 10 Myr (`timeOffsetMaximumAbsolute`) of the merge time, so
   each merger fragments the host's natural ~0.1/H stepping into event-sized
   pieces, and every fragment costs a full O(live nodes) walk. Merger +
   destruction event counts scale as the subhalo count S ~ (M/m_res)^0.96
   (measured: 133 subhalos/tree at z=0 for M/m=10^4; 1,204 at 10^5).

4. **Hypothesis check:** (a) "subhalos build up and are walked every time
   step" - confirmed, with the refinement that the walk itself is cheap and
   the real costs are (i) the O(walks x live-nodes) visit product, and (ii)
   whether each visited satellite *evolves* (pays ~130 us ODE-call overhead)
   or *defers*. (b) "branches become more finely resolved in time" - ruled
   out: built node count scales as (1/m_res)^0.95 (slightly sublinear), and
   nodes/branch stays constant.

5. **The single biggest lever is `fractionTimestepSatelliteMinimum`.** The
   reference file sets 0.75, but the **code default is 0.0**. With the default,
   every stalled-walk visit to a satellite takes a micro ODE call with full
   per-call overhead: at M/m=10^5 this alone is **8.3x slower evolution**
   (142 s vs 17 s per tree), with identical trees and walk counts
   [defaults-series exponent: see section 7]. A production configuration
   using the default would see the ~1/m^1.5 (and eventually ~1/m^2) regime
   one to two decades earlier in m_res. `backtrackToSatellites`, the merger
   sync window (x10 either way), and `reuseODEStepSize` are all neutral at
   these scales (<= 4%).

6. Per-satellite ODE cost is dominated by tidal-heating physics
   (`tidalTensorPathIntegrated` limits most steps; RHS evaluations are
   dominated by heated-profile Brent root finds and NFW velocity
   dispersions) - a constant-factor (not scaling) cost, and the natural
   target for making the *linear* term cheaper.

## 2. How tree evolution spends its time (code map)

File/line references to master @ d6a011b.

### The walk loop

`mergerTreeEvolverStandard::standardEvolve`
(`source/merger_trees/evolver/standard.F90:291`) evolves a tree by repeatedly
walking **every node** (post-order DFS, satellites -> children -> parent;
`source/merger_trees/walkers/all_nodes/_class.F90:101`) until a walk evolves
nothing. Per visited node it runs `nodeIsEvolvable` (O(1), l.723) and, if
evolvable, `timeEvolveTo` (l.762) followed by an ODE evolution
(`mergerTreeNodeEvolver`). A node that advances is re-evolved until blocked
(l.517); `backtrackToSatellites` also revisits its satellites (l.514).
Promotions continue up a branch within one walk (l.501), so branch-internal
evolution does not multiply walk count - events do.

### Per-visit scans (all O(S) for the final host)

`standardTimeEvolveTo` runs, each time a node is considered:

- the timestep classes (`simple` computes expansion factor + rate each call;
  `satellite` merger sync; `satelliteDestruction`),
- a sibling scan for primary progenitors (l.946),
- a **hosted-satellite scan** - O(#satellites) per host visit (l.963),
- a **mergee scan** - O(#mergees) (l.990; the stray comment at l.989 - "need
  to also figure out how to speed up satellite evolve checks" - suggests this
  was already suspected).

Further O(S) linked-list operations: satellite insertion via `lastSatellite()`
(walks the whole list; `source/satellites/promotion.F90` and
`node_evolver/standard.F90:1293`), removal via `removeFromHost` (list walk),
and **every trunk promotion re-parents the entire satellite list** and appends
mergee lists (`standardPromote`, `node_evolver/standard.F90:1284-1320`). With
trunk promotions ~ 1/m_res and S ~ 1/m_res these are intrinsically quadratic,
though with tiny constants (pointer writes).

### Event synchronization

Mergers trigger as mid-ODE interrupts when r < 0.01 r_vir
(`nodes/operators/physics/satellite_merging/radius_trigger.F90` sets
timeOfMerging = now); the satellite then waits until the merge target is
within `timeOffsetMaximumAbsolute` = 0.01 Gyr (`evolve/timesteps/satellite.F90`),
and the target's step is capped at the merge time (evolver l.990). Each merger
therefore forces the host to land a step in a ~10 Myr window - fragmenting its
~0.1/H (~1 Gyr at z=0) natural step into event-sized pieces, each of which
costs one full-tree walk. Destructions (mass < threshold,
`satellite_destruction/mass_threshold.F90`) interrupt only the satellite's own
step. Event counts scale as S.

### Per-ODE-call fixed overhead

Every `nodeEvolver%evolve` call (`node_evolver/standard.F90:385`) pays: tree
locks, pre/post hooks over all node operators, `Calculations_Reset` (event-hook
fan-out to every subscribed class), 3x serialization-offset passes, 2x
serialize + saved copies, scales pass, **a freshly allocated GSL ODE driver**,
and (with `reuseODEStepSize=false`) an initial-step-size search. Measured
~130 us per ODE step, independent of step size - a tiny 10 Myr step costs the
same as a 1 Gyr step. This constant is what makes forced micro-steps so
expensive.

## 3. Headline measurements

10^13 M_sun host, 4 trees per point, per-tree means:

| m_res [M_sun] | M/m_res | nodes | t_construct [s] | t_evolve [s] | walks | slope(nodes) | slope(t_evolve) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1e10   | 1e3    | 1,604     | 0.41  | 0.18   | 8     | -    | -    |
| 3.16e9 | 3.16e3 | 4,890     | 0.39  | 0.57   | 10    | 0.97 | 1.00 |
| 1e9    | 1e4    | 14,484    | 0.63  | 1.79   | 16    | 0.94 | 0.99 |
| 3.16e8 | 3.16e4 | 44,122    | 1.69  | 5.54   | 52    | 0.97 | 0.98 |
| 1e8    | 1e5    | 132,989   | 4.15  | 17.10  | 195   | 0.96 | 0.98 |
| 3.16e7 | 3.16e5 | 398,922   | 12.47 | 54.65  | 645   | 0.95 | 1.01 |
| 1e7    | 1e6    | 1,205,244 | 36.73 | 194.36 | 3,659 | 0.96 | 1.10 |

- t_evolve per built node: 112-137 us, constant to M/m = 3x10^5; 161 us at
  10^6 - the first departure.
- Walk count slope reaches ~1.5/decade in the last step; walks x median live
  nodes (total visits) scales as (1/m_res)^2.06.
- Construction stays ~linear (slope 0.94-0.96) throughout; at every measured
  resolution construction is 19-24% of total - it does not drive the break.
- ODE-step diagnostics (evolve profiler): steps, RHS evaluations, and in-step
  CPU all linear through M/m = 3x10^5 with ~1.7 evaluations/step and ~130
  us/step; step sizes are limited overwhelmingly by
  `satellite:orbiting:tidalTensorPathIntegrated`; ~33% of steps end in
  interrupts. The small-step tail (dt < 1 Myr) grows slowly (11% -> 16% of
  steps from M/m = 1e4 -> 3x10^5).

### 10^12 M_sun host series (host-mass invariance + deeper probe)

| m_res [M_sun] | M/m_res | nodes | t_evolve [s] | walks | slope(t_evolve) |
|---:|---:|---:|---:|---:|---:|
| 1e8    | 1e4    | 13,440    | 1.65   | 10    | -    |
| 1e7    | 1e5    | 121,062   | 15.84  | 159   | 0.98 |
| 3.16e6 | 3.16e5 | 369,846   | 52.28  | 864   | 1.04 |
| 1e6    | 1e6    | 1,112,300 | 178.66 | 3,325 | 1.07 |

At matched M/m_res the 10^12 and 10^13 series agree to better than 10% in
evolve time, node count, and walk count: **cost is a function of M_host/m_res,
essentially independent of host mass**. Both series depart from slope 1.0 at
M/m_res ~ 3x10^5 - 10^6. (A cluster-mass host therefore reaches the break at
proportionally coarser m_res: 10^15 M_sun hosts hit M/m = 10^6 already at
m_res = 10^9.)

### Two-component cost model

All twelve (M/m_res, t_evolve) points from both host-mass series are fit to
0.019 dex rms by

    t_evolve = A (M/m_res) + B (M/m_res)^2,  A = 1.71e-4 s, B = 1.34e-11 s

(single-thread seconds on the test machine; the *ratio* B/A = 7.9e-8 is the
transferable quantity). The local logarithmic slope is then 1.1 at M/m_res =
1.4e6 (matching the measured onset), **1.5 at M/m_res = 1.3e7**, and 1.75 at
3.8e7, tending to 2 asymptotically. A measured slope of ~1.5 is the signature
of being mid-crossover, not a true asymptotic exponent - the asymptote is
quadratic.

Control: per-node *construction* cost stays flat (~30 us/node) all the way to
the 3.4M-node / 9 GB working set, so the evolve-side rise is algorithmic, not
a memory-hierarchy effect.

### Walk anatomy

Per-walk progress reports (verbosity=warn) show the two regimes directly.
Healthy walks evolve ~80% of live nodes. Stalled walks - the merger-sync
cascades - visit everything to evolve almost nothing (samples from one
m_res=1e8 tree: 1,262 live/21 evolved, 1,259/18, ..., 1,242/0). At m_res=1e8,
11-100% of a tree's total visits (tree-to-tree scatter) occur in walks
evolving <10% of live nodes.

### Subhalo census

z=0 subhalo counts per tree: 133 (M/m=1e4), 1,204 (M/m=1e5) -> S ~
(M/m)^0.96. The population is roughly steady from z~1.5 to 0 (accretion
balancing destruction), so the host carries an O(S) satellite list for most of
its evolution.

## 4. Variant experiments (m_res = 1e8, t_evolve per tree)

| variant | t_evolve [s] | walks | vs baseline |
|---|---:|---:|---:|
| baseline (reference config) | 17.10 | 195 | 1.00 |
| `fractionTimestepSatelliteMinimum=0` (code default) | **142.09** | 195 | **8.31x** |
| `backtrackToSatellites=false` | 17.00 | 195 | 0.99 |
| `reuseODEStepSize=true` | 16.43 | 216 | 0.96 |
| merger sync window x10 (0.1 Gyr) | 17.05 | 195 | 1.00 |
| merger sync window /10 (1 Myr) | ~17.1 | 195 | 1.00 |

The no-defer result is the decisive one: identical trees and walk structure,
but satellites take a micro ODE call at every visit instead of deferring until
they can take >= 0.75 of the host-relative step. All the penalty is per-call
fixed overhead on tiny steps. **Any configuration that leaves the parameter at
its default 0.0 pays this everywhere** - and its cost share grows with the
stalled-walk visit product, i.e. superlinearly.

## 5. What is *not* responsible

- **Tree construction temporal refinement**: node count is linear in 1/m_res
  (slope 0.95-0.97); nodes per branch do not grow with resolution.
- **Tree construction cost**: linear, and a constant ~20% share.
- **The merger sync window size**: neutral at x10 either way (at M/m <= 1e5).
- **Backtracking**: neutral.
- **ODE tolerance pathologies**: per-step cost and evaluations/step constant.

## 6. Profiling

gdb sampling during evolution shows the dominant linear-term payload inside
satellite ODE RHS evaluations: tidal-heating specific-energy terms ->
`massDistributionSphericalHeated` initial-radius Brent root-finds -> NFW
velocity dispersion (log/fma), plus mass-distribution object copies
(`memmove` in `darkMatterOnlyGet`).

TODO(batch 3): perf record shares at M/m = 1e4 / 1e5 / 1e6 - quantify the
growth of walker + timeEvolveTo + serialization shares with resolution.

## 7. Defaults-series scaling: the observed 1/m^1.5 reproduced

Same model, but with the evolver's **code-default** settings
(`fractionTimestepSatelliteMinimum=0`, `backtrackToSatellites=false`) instead
of the reference-file values (0.75/true); 10^13 M_sun host:

| m_res | M/m_res | t_evolve [s] | vs reference | slope(t_evolve) |
|---:|---:|---:|---:|---:|
| 1e9    | 1e4    | 3.89   | 2.17x | -    |
| 3.16e8 | 3.16e4 | 20.63  | 3.72x | 1.45 |
| 1e8    | 1e5    | 140.75 | 8.23x | 1.67 |

Global power law over this range: **t_evolve ~ (1/m_res)^1.56** - the observed
production scaling, reproduced at practical resolutions. The mechanism is the
one quantified in section 4: without deferral, every stalled-walk visit to a
satellite pays the full ~130 us ODE-call overhead for a micro-step, so the
walk-driven quadratic term carries a ~100x larger coefficient and its
crossover moves from M/m ~ 1e7 down to M/m ~ 3e4. A "1/m^1.5" measurement is
this crossover in progress; the asymptote in both configurations is
quadratic.

## 8. Strategies to improve the scaling

Ranked by expected impact on the *exponent* first, then on the constant.

**S1. Event-driven evolution (priority queue) - the structural fix.**
Replace "walk everything until nothing evolves" with a queue of evolvable
nodes keyed by allowed evolve-to time. The infrastructure already exists:
`timeEvolveTo` computes, for every blocked node, *which* node blocks it
(`nodeLock`, today used only for deadlock reporting). Invert it: per-node
waiter lists; when a node advances or a tree event fires, wake exactly its
waiters. Complexity O((evolutions + wakes) log N) replacing O(walks x
N_live) - removes the stalled-walk regime entirely (walks that visit 1,300
nodes to evolve 20). This converts the visit product's (1/m)^2.06 into
~(1/m)^1.0 log.

**S2. Keep satellites deferring (config + maybe default change).**
`fractionTimestepSatelliteMinimum=0.75` is worth 8.3x at M/m=1e5 with no
measured accuracy concern in this model class (it is the reference/validation
setting). Consider making ~0.75 the code default, or warning when 0 is used
with orbiting satellites. Cheapest possible win for anyone running defaults.

**S3. Cache the host's minimum-satellite time and satellite-list tail.**
The O(S) hosted-satellite scan runs every time the host is considered; a
cached min-time (invalidated lazily: recompute only when the cached min's
holder advances/leaves) makes it O(1) amortized. Similarly a `lastSatellite`
pointer (or doubly-linked satellite list) makes insertion/removal O(1), and
trunk promotions can splice the satellite list without re-walking it if the
parent pointer is resolved through the host rather than stored per satellite
(or re-parent lazily). Together these remove the intrinsically quadratic
pointer work.

**S4. Reduce per-ODE-call fixed overhead.**
Reuse a per-thread GSL driver/workspace keyed by system dimension instead of
alloc/free per call; skip re-serialization scaffolding when the property
count is unchanged (the common case for a satellite stepping repeatedly);
`reuseODEStepSize=true` is a free ~4%. This scales down the coefficient of
*both* the linear term and the forced-micro-step term.

**S5. Batch merger synchronization.**
Process mergers due within one sync window together: let the host take one
step to the earliest merge time and complete *all* mergers falling within the
window, instead of one walk per merger. Reduces event-driven walk count by
the mean number of mergers per window (grows with S, so this directly attacks
the superlinear term when S1 is too invasive).

**S6. Cheapen the satellite RHS (linear-term constant).**
Memoize/interpolate the heated-profile initial-radius solution (currently a
Brent root-find per evaluation) and NFW dispersion lookups; avoid
mass-distribution deep copies per call (`darkMatterOnlyGet` memmove). Given
tidal physics dominates RHS cost, plausibly a 1.5-3x constant-factor win for
subhalo-heavy models.

**S7. Physics-level coarse-graining (already available).**
The validation configuration subsamples branches below a mass threshold
(`mergerTreeBuildController=subsample`) - the established way to cut S when
full subhalo statistics below a mass scale are not needed.

## 9. Reproduction

```bash
export GALACTICUS_EXEC_PATH=... GALACTICUS_DATA_PATH=... GALACTICUS_TOOLS_PATH=...
cd performanceAnalysis/scripts
./run_one.sh base_m1.0e8 1.0e8 4                 # one point of the base series
python3 collect.py base_m1.0e8                   # harvest treeTiming + logs
python3 analyze_scaling.py base                  # slopes table
./run_one.sh diag_m1.0e8 1.0e8 4 profiler=simple profileOdeEvolver=true
python3 profcompare.py diag_m1.0e8               # ODE-step histograms
./profile_one.sh prof_m1.0e8 1.0e8 2             # perf sampling profile
```

`make_params.py` builds fully-resolved parameter files from the reference
includes (paths are environment-specific - edit `REFDIR`); the E&H transfer
function is substituted for CAMB. Timing instrumentation is
`mergerTreeOperator=treeProcessingTimer` (per-tree construct/evolve CPU and
node counts in `metaData/treeTiming`) plus `verbosityLevel=warn` walk logs.
