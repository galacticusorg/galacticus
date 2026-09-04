# Galacticus performance survey (September 2026)

A code-reading survey of the Galacticus Fortran sources and the Python build
tooling, looking for micro- and macro-optimizations *beyond* those already
proposed in issues [#1431](https://github.com/galacticusorg/galacticus/issues/1431)
(satellite deferral default), [#1432](https://github.com/galacticusorg/galacticus/issues/1432)
(O(1) satellite bookkeeping), [#1433](https://github.com/galacticusorg/galacticus/issues/1433)
(ODE driver reuse / node-evolver fixed overhead),
[#1434](https://github.com/galacticusorg/galacticus/issues/1434) (batched merger
synchronization) and [#1435](https://github.com/galacticusorg/galacticus/issues/1435)
(event-driven evolver). Those five attack the tree-evolver *scheduling* layer;
this survey covers everything else.

Tree surveyed: `master` @ e8d9c46. Line numbers refer to that tree.

## 0. How to read this document

* **Measured [M]** items were timed or profiled in this survey (Python build
  tooling only — no Fortran compiler was available, so no model was run).
  Run-time findings are grounded where possible in the `perf` profiles from the
  earlier scaling analysis
  ([`claude/galacticus-performance-analysis-giadol`](https://github.com/galacticusorg/galacticus/tree/claude/galacticus-performance-analysis-giadol/performanceAnalysis),
  reference dark-matter-only subhalo model), whose self-time shares were:
  satellite/tidal physics RHS ~15%, math library ~10%, ODE core ~9%,
  mass distributions ~8%, memory alloc/copy ~7%. Everything else here is
  **inferred [I]** from reading the code, with frequencies stated (per RHS
  evaluation, per node per step, per tree, once per model).
* "Per RHS" means per ODE right-hand-side (derivative) evaluation. A node sees
  O(10–100) RHS evaluations per timestep, so a per-RHS cost is the most
  expensive kind.
* **Bitwise-neutral** means the change cannot alter model output; such changes
  can be validated by `h5diff` against `master`. Non-neutral changes need
  statistical validation.
* Effort: S (hours), M (days), L (a week or more).
* **Overlap with open issues.** The open issue list was checked after the
  survey. Items already covered by an open issue are marked **[filed: #N]**
  and kept only for context; the rest are new. Related open optimization
  issues: [#701](https://github.com/galacticusorg/galacticus/issues/701) and
  [#913](https://github.com/galacticusorg/galacticus/issues/913) (heated
  profiles, §2.3), [#575](https://github.com/galacticusorg/galacticus/issues/575)
  (adiabatic-contraction root find, §2.2/§5.3),
  [#1417](https://github.com/galacticusorg/galacticus/issues/1417) (per-call GSL
  interpolator in history interpolation, a sibling of §4.2–4.3),
  [#1419](https://github.com/galacticusorg/galacticus/issues/1419) (tabulated
  interpolate tidy-up), [#111](https://github.com/galacticusorg/galacticus/issues/111)
  (`integration2` clean-up, touches §4.3) and
  [#58](https://github.com/galacticusorg/galacticus/issues/58) (lightcone
  luminosities, related to §3.1).

## 1. Executive summary — the ten most valuable items

| # | Item | Kind | Expected gain | Effort | Neutral |
|---|------|------|---------------|--------|---------|
| 1 | `Calculations_Reset` fires **per RHS** and its mass-distribution handler discards *all* cached distributions; combined with 29 non-pooled profile `get` implementations this makes several object allocations + destructor chains per RHS per node (§2.1, §2.2) | run-time, macro | medium (part of the measured ~7% alloc/copy and ~8% mass-distribution shares) | S per class, M for smarter invalidation | yes |
| 2 | Heated-profile `radiusInitial` is a Brent root find with a one-entry memo that thrashes within one RHS (§2.3). **Already filed**: [#701](https://github.com/galacticusorg/galacticus/issues/701) carries a full tabulate-and-invert plus Newton design, and [#913](https://github.com/galacticusorg/galacticus/issues/913) the shared cross-subhalo table; listed here only because it remains the largest per-satellite constant-factor lever | run-time, macro for subhalo models | plausibly 1.5–3× on the satellite RHS (dominant linear term) | M | no |
| 3 | `rootFinder`/`integrator` objects constructed on the stack per call — ~8–10 heap ops plus two OpenMP lock inits each; the freefall-radius path nests *three* levels of them (§4.1–4.3) | run-time | small–medium, macro on decorated/composite profiles | S–M | yes |
| 4 | Equilibrium structure solver always runs ≥2 full iterations per RHS (§5.3) | run-time, macro for baryonic models | medium–large per RHS in baryonic runs | M | partly |
| 5 | HDF5 output: dataset re-probed and re-opened (2× `h5dopen`) for every property × tree × output under the global lock; output buffer grows by +1024 (quadratic); per-tree metadata rebuilt (§5.1, §5.2, §5.7) | run-time I/O, macro for many-small-tree and N-body runs | large where output dominates | S–M | yes |
| 6 | `make` spends ~6 s per invocation trying built-in RCS/SCCS pattern rules on ~2000 prerequisites; `.SUFFIXES:` does not cancel them (§6.1) [M] | build | 6 s per make, 13× faster no-op builds | S | n/a |
| 7 | Quadratic string appends in the preprocessor's serializer and comment embedder: ~36 s on the node-component module alone; whole-tree `deepcopy` in Generics: ~20 s on the HDF5 module (§6.2, §6.3) [M] | build | ~55 s off the critical path of every clean build | S | n/a |
| 8 | FunctionClass parents re-parse **and re-process** every implementation file; catalog XML parsed up to 13 times per invocation; ~0.17 s fixed cost × 2000 invocations (§6.4–6.6) [M] | build | roughly half of the ~680 CPU-s preprocessing total | M | n/a |
| 9 | Object-typed `<prop>Rate` accumulators (`stellarLuminosities`, `abundances`, …) allocate and copy the serialized object three times per call; only `history` uses the in-place path (§3.1) | run-time | 1–5% of RHS time in luminosity-heavy runs | S | yes |
| 10 | Release builds target baseline x86-64 with no `-march`, no PGO; the global `Display_Lock` is taken *before* the verbosity test on every message (§7, §3.2) | build flags / run-time | 5–15% (PGO/`-march`, to be measured); removes a global serialization point | S | flags: no; lock: yes |

## 2. Run time: the satellite / dark-matter-profile hot path

The earlier profiles put the linear (per-node) cost of the reference subhalo
model in satellite tidal physics → heated mass distribution → NFW kinematics,
plus allocation churn. Since that analysis, `master` has gained an object pool
(`source/objects/pool.F90`) and pooled the NFW, heated, tidal-heating and
beta-profile factories, and the disk/spheroid component distributions. The
following remain.

### 2.1 `Calculations_Reset` is per RHS and flushes every cached mass distribution [I, verified]

* `source/merger_trees/node_evolver/standard.F90:1042` — `call Calculations_Reset(node)`
  inside `standardDerivativesCompute`, i.e. **every derivative evaluation**
  (issue #1433 discusses only the per-evolve-call site at `:452`).
* `source/objects/nodes/_class.F90:1863-1881` — `massDistributionCalculationReset`
  destroys all 20 slots of the memoized `massDistributions__` array and its
  `uniqueID` argument is explicitly `unused`.
* Consequence: every `node%massDistribution(...)` chain (dark matter → heated →
  NFW → scaler → kinematics, plus hot halo, disk, spheroid where present) is
  rebuilt on every RHS. Pooled classes pay a re-`initialize`; non-pooled
  classes pay `allocate` + constructor + destructor + `deallocate` for each
  object in the chain.
* Fix (bitwise-neutral): (a) finish pooling (§2.2); (b) in the reset handler,
  keep distributions whose inputs cannot have changed within a step — the
  dark-matter-only distribution depends only on `basic%mass`, the scale radius
  and time, all inactive between steps, so it can be kept when those match the
  cached values for the same `uniqueID`. This is the natural home for the
  generation-counter idea in #1433(D).

### 2.2 Non-pooled `get` factories [I, verified]

Only four factories use `objectPool` (`dark_matter_profiles_DMO/NFW.F90`,
`dark_matter_profiles_DMO/heated/_class.F90`, `dark_matter_profiles/heating/tidal.F90`,
`hot_halo/mass_distribution/beta_profile.F90`). The 29 `get(self,node,weightBy,...)`
implementations without a pool include, in likely order of use:

* `dark_matter_profiles/dark_matter_only.F90:125-196` — allocates a
  `massDistributionSphericalScaler` **and** a `kinematicsDistributionSphericalScaler`
  per call, wrapping a pooled DMO object. This is in the reference subhalo
  model's chain, and is the `darkMatterOnlyGet` `memmove` seen in the profile.
* `dark_matter_profiles/adiabatic_Gnedin2004.F90:229,265,277` — three
  allocations per call; used by `parameters/quickTest.xml` and most baryonic
  models. The assessment on [#575](https://github.com/galacticusorg/galacticus/issues/575)
  measured the adiabatic cluster at ~2% of the milkyWay benchmark, dominated by
  `constructorInternal`/`get`/`initialize` rather than the root find, and
  called repeated construction "the better target" — this is that item.
* `dark_matter_profiles_DMO/{Einasto,Burkert,Zhao1996,isothermal,cusp-NFW,decaying,
  soliton_NFW,soliton_NFW_heated,multiple,accelerator,SIDM_parametric}.F90`,
  `dark_matter_profiles_DMO/{finite_resolution,truncated,heated/monotonic,
  accretion_flow/*,SIDM/*}`, `dark_matter_profiles/{SIDM/isothermal,accelerator}.F90`,
  `hot_halo/mass_distribution/{PatejLoeb2015,Enzo_hydrostatic,Ricotti2000,null}.F90`.

Fix: apply the existing pool pattern (see `NFW.F90:150-200` for the template).
Note the stale comment in `heated/_class.F90:heatedGet` ("the heating object is
currently constructed afresh on each call") — it is now pooled. Effort S per
class. Bitwise-neutral.

### 2.3 Heated-profile `radiusInitial`: root find per new radius [I, verified] **[filed: #701, #913]**

* `source/mass_distributions/spherical/heated/_class.F90:369-437`
  (`sphericalHeatedRadiusInitial`): a one-entry memo (`radiusFinalPrevious`)
  then a Brent solve (`self%finder%find`) with range expansion. Every
  `density`, `massEnclosedBySphere`, `tidalTensor` (which calls both) and the
  heated kinematics Jeans integrand (`kinematic_distributions/heated.F90:163,266`)
  go through it. Within one RHS the profile is queried at several distinct
  radii (orbital radius, virial radius, half-mass radius, tidal radius), so the
  memo thrashes and several Brent solves run per RHS per satellite. The profile
  showed "Brent + root finder ~7%" and "NFW dispersion ~7%" (the latter is the
  root function's inner evaluations).
* The forward map is explicit: `sphericalHeatedRadiusEnclosingMass` (`:505-530`)
  already uses `r_f = 1/(1/r_i − 2ε(r_i)/(G M(r_i)))`. So `r_i(r_f)` is the
  inverse of a cheap monotonic function.
* **Already filed.** [#701](https://github.com/galacticusorg/galacticus/issues/701)
  (with a 2026-08-19 design comment) proposes exactly this: Design A,
  per-object tabulate-and-invert of the explicit forward map, gated by an
  evaluation-count heuristic; Design B, a derivative-based `fdfsolver` solve
  using the analytic `dr_f/dr_i`, with the bracketing solver as fallback;
  recommended A+B. [#913](https://github.com/galacticusorg/galacticus/issues/913)
  goes further for the tidal-heated tabulated-NFW case: fold the heating
  parameter `Q/(G ρ_s)` in as one extra axis of a *shared*
  `massDistributionSphericalTabulated` table, amortized across every subhalo
  in the run. Nothing to add beyond noting that the per-RHS cache flush in
  §2.1 means #701's per-object table is rebuilt every RHS unless the
  invalidation there is also refined — #913's run-wide table avoids that
  coupling entirely.
* Related: `Gnedin1999.F90:195-206` also evaluates `radiusEnclosingMass` and
  `rotationCurve` at the half-mass radius per RHS; with the table in place these
  become cheap.

### 2.4 ODE scale of the integrated tidal tensor [I]

`source/objects/nodes/components/satellite/orbiting.F90:217,248-252` sets the
absolute scale of the six `tidalTensorPathIntegrated` components to 1% of the
virial value. The earlier ODE-step diagnostics found this property limits most
satellite steps. Whether 1% is needed for the tidal-heating accuracy actually
used downstream (`tidalHeatingNormalized`, itself scaled at 1%) is a physics
question; exposing the fraction as a parameter would make the experiment cheap.
Not bitwise-neutral; potentially large step-count savings; validate against the
DMO subhalo statistics.

## 3. Run time: core infrastructure

### 3.1 Object-typed rate accumulators copy the whole object three times per call [I, verified]

`python/Galacticus/Build/Components/Properties/Evolve.py:20-27,372-381`: only
`history` is in `INCREMENT_SERIALIZED_TYPES`. Every other non-intrinsic
property rate (`stellarLuminosities`, `abundances`, `chemicalAbundances`)
does `current=self%<prop>Data; deserialize; increment; serialize`, i.e. one
malloc/free plus three full copies per rate call, per component, per RHS. For
`stellarLuminosities` with hundreds of bands this is a measurable share of the
RHS. Fix: add `incrementSerialized` to those three types (mirror
`History_Increment_Serialized`, `history.F90:1196`) and list them; the
generator already emits the fast path. Effort S. Bitwise-neutral.

### 3.2 `Display_Lock` taken before the verbosity check [I, verified]

`source/display/_module.F90:282-284,312-314` (`displayMessageChar`,
`displayMessageVarStr`, and the indent/unindent routines): the global
`!$omp critical(Display_Lock)` is entered, then `showMessage(verbosity)` is
tested. Every info/debug-level message anywhere serializes all threads even
when nothing prints. Some callers also build the `varying_string` message
eagerly (e.g. `objects/nodes/_class.F90:786-788`, `Tree_Node_Remove_from_Mergee`,
whereas its twin at `:752-756` guards on `displayVerbosity()`). Fix: test
verbosity outside the critical (the level is set once at start-up), and guard
eager string construction at hot call sites. Effort S. Bitwise-neutral.

### 3.3 Per-node memory footprint from the component-array design [I]

`treeNode` holds one `class(nodeComponentX), allocatable, dimension(:)` per
component class (~80 B descriptor each, even when empty);
`treeNodeXGet` allocates a generic placeholder component on first query of a
class the node lacks (`Components/TreeNodes/Classes.py:155-159`); and
`Implementation_Creation` emits `allocate(self%<prop>Data(0))` for every
rank-1 property (`Components/Implementations/CreateDestroy.py:137-149`).
Together these plausibly account for a large fraction of the ~3 kB/node and of
malloc traffic during tree build/destroy. Fix: shared immutable per-class
placeholder instead of per-node allocation; drop zero-size allocations (the
`Count` functions already treat unallocated as zero). Effort M–L (many
`allocated(...)` tests in generated code). Bitwise-neutral. Macro for
memory-bound read-tree runs.

### 3.4 `treeNode%isSatellite()` is an O(#siblings) scan [I, verified]

`objects/nodes/_class.F90:661-680`: walks `parent%firstChild → sibling…`
looking for `self`. 140 call sites, including `basic%timeLastIsolated()` (66
callers) and every tree-walker ascent. Fix: cache a logical on `treeNode`,
maintained at the (few) list-mutation sites (#1432 already enumerates them).
Gain micro on built trees (short child lists), moderate on branchy N-body
trees. Effort M. Bitwise-neutral.

### 3.5 Smaller infrastructure items

* `nodeEventIncrement` uses `!$omp critical` (`_class.F90:475-478,1531-1534,1720-1723`)
  whereas the uniqueID counter was already converted to `!$omp atomic capture`
  (`:395-398`). Same fix. S.
* `reportMemoryUsage()` per tree per thread (`tasks/evolve_forests/_class.F90:788,993`)
  calls glibc `mallinfo2`, which locks and walks every malloc arena. Throttle
  by wall clock. S. Small–medium at high thread counts with tiny trees.
* Whole-tree walks for `earliestTime()/latestTime()` per evolve iteration
  (`_class.F90:1555-1616`; callers `evolve_forests/_class.F90:838,876,966`) —
  compute once per iteration and share. S.
* Memoized `massDistribution` lookup compares five keys over 20 slots with no
  guaranteed short-circuit (`Components/TreeNodes/Utils.py:365-390`); nest the
  `uniqueID` test. S, micro but hot.
* ODE serialization recomputes `serializeCount` per component and uses `pack()`
  with masks (`Components/TreeNodes/ODESolver.py:118-160,264-274`); a per-step
  gather index would remove O(P) redundant work per RHS. M, micro–small.
* Rank-1 meta-property `Get` allocates and copies where `GetReference` exists
  (`Components/Classes/MetaProperties.py:66-77,105-171`). Audit hot callers. S.
* Object-valued `+`/`*` operators on `stellarLuminosities`/`abundances`
  allocate on return (`structure.F90:716-743`, `abundances.F90:414-431`);
  prefer `increment` in per-RHS code. S per site.
* Verified fine: property getters are plain copies; debug guards compile only
  under `-DDEBUGGING`; event-hook dispatch is lock-free; the functionClass
  graph is deep-copied once per thread, not per tree; uniform tables use O(1)
  index arithmetic with last-value caching; no run-time `inputParameters`
  lookups occur in hot paths (keep it that way).

## 4. Run time: numerical kernels

The dominant pattern is **GSL-backed `integrator`/`rootFinder` objects
constructed on the stack per call**. Each costs `gsl_integration_workspace_alloc(1000)`
(~48 kB, five arrays) or `gsl_root_fsolver_alloc`, a `malloc` for the
`gsl_function`, and two `resourceManager` constructions (each allocating a
counter and an OpenMP lock, `utility/resource_manager.F90:109-113`) — roughly
8–10 heap operations plus two lock init/destroy pairs per object lifetime. The
`rootFinder` type already supports persistent, recursion-safe use (threadprivate
`currentFinders` stack, `numerical/root_finder.F90:215-230`), so the fix is
usually just to hoist the object.

### 4.1 Triple-nested integration in the freefall radius [I, verified]

`mass_distributions/spherical/_class.F90:1209-1223`: `rootRadiusFreefall`
constructs an `integrator` **inside the root function**; its integrand calls
`potentialDifference`, which for profiles without an analytic override
(`sphericalPotentialDifferenceNumerical`, `:769-797`) constructs another
`integrator`. Per freefall radius: Brent × QAG61 × (alloc + QAG61) ×
`massEnclosedBySphere` ≈ 10⁴–10⁵ enclosed-mass evaluations and ~10³ workspace
mallocs. Called per node per cooling-rate RHS via
`cooling/freefall_radii/dark_matter_halo.F90:150`. NFW/Burkert/Einasto/isothermal
override `radiusFreefall` with a tabulation; decorated, heated, SIDM,
finite-resolution and composite profiles do not. Fix (in order): hoist both
integrators to `save`d threadprivate objects (neutral, S); compute
`potentialDifference` from the existing potential table (M); tabulate
`t_ff(r)` once per distribution as `nfwTimeFreefallTabulate` does (M).

### 4.2 Per-call `rootFinder` in the generic radius solvers [I, verified]

* `mass_distributions/_class.F90:339-367,393-421,449-467,496-514,541-559,611-630`
  (`radiusEnclosingMassNumerical`, `radiusEnclosingDensityNumerical`,
  `radiusFromSpecificAngularMomentumNumerical`, `radiusRotationCurveMaximumNumerical`, …),
  `mass_distributions/spherical/_class.F90:465-499`, `finite_resolution/NFW.F90:752,1051`,
  `satellites/dynamical_friction/acceleration/Kaur2018.F90:232`,
  `satellites/orbits.F90:116,245`.
* Hot callers: `satellites/tidal_stripping/radius/King1962.F90:266,277` (per
  satellite per RHS), `satellite_merging/radius_trigger.F90:305,317`,
  `galactic_structure/radius_solver/equilibrium.F90:369`.
* Several use `rootGuess=1.0d0` (1 Mpc) for kpc-scale galaxies, forcing 7–10
  bracket-expansion evaluations before Brent starts; only
  `radiusEnclosingDensityNumerical` keeps a previous-solution guess.
* Fix: a threadprivate stack of prebuilt finders parallel to the existing
  `massSolvers` stack (neutral, S–M); previous-solution guesses for the other
  methods (not neutral, S). [#1417](https://github.com/galacticusorg/galacticus/issues/1417)
  addresses the same per-call-construction pattern for the `interpolator` in
  history interpolation; these could share one fix PR.

### 4.3 Other per-call construction sites [I, verified]

* `potentialDifferenceNumerical` (`_class.F90:1345`, `spherical/_class.F90:786`)
  — per satellite per step via `satellites/orbits.F90:376` and black-hole code.
* `concentrationRadius` (`dark_matter_profiles/structure/scale/concentration.F90:244-262,296-297`)
  creates and destroys a work `treeNode` with two components plus a
  `rootFinder` per call; cache them in the existing threadprivate `state_`.
* Zhang & Hui table build constructs a QAGS integrator per (i,j) pair
  (`Zhang_Hui.F90:422-432`), once per model but O(N²) times.
* `multiVectorizedCompositeGaussKronrod1DEvaluate` (`numerical/integration2.F90:1637-1671,1720-1811`)
  allocates interval nodes (8 allocations per bisection) and sorts through a
  C-callback heapsort on every accepted ODE step of every node when
  `useJacobian` is on (via `latentIntegrator`, `ODE_solver/solver.F90:640-668`).
  Pre-allocate an interval arena. M, neutral.
* ODE `solve` allocates `z0` per call (`solver.F90:559-564`). S.

### 4.4 Finite-difference Jacobian recomputes the base derivative [I, verified]

`node_evolver/standard.F90:975-980`: `standardJacobian` calls
`standardDerivativesCompute` at the unperturbed state before the N perturbed
calls. The vendored `msbdfactive_corrector` (`external/gslODEInitVal2/msbdfactive.c:780-798`)
has just evaluated `f(t+h, ytmp)` at the same predicted point before the
Jacobian request. Passing `dydt_in` through (the stepper is vendored, so the
signature is ours to change) removes one of N+1 full node-operator chains per
Jacobian (~3–9% of Jacobian cost for N≈10–30); writing columns directly into
`derivativeRatesValues` removes the `reshape` copy at `:1015`. S–M. Neutral.

### 4.5 Vendored GSL compiled without `-DHAVE_INLINE -DGSL_RANGE_CHECK_OFF` [I, verified]

`Makefile:149` sets no GSL macros, so the Newton loops in
`msbdfactive.c:834-882` call out-of-line, range-checked
`gsl_vector_get/set`/`gsl_matrix_get/set` O(dim) times per iteration per step
per node. LTO cannot inline across the shared `libgsl` boundary, but these
are in *our* translation units, so the macros apply. Also
`lu.c:82,122` allocates `ipiv` per LU decomposition. S. Neutral.

### 4.6 Table and interpolation micro-items

* **Logarithmic-table range check is both slow and inconsistent** [verified]:
  `objects/tables/_module.F90:2534` tests `x > self%x(self%xCount)` — a
  virtual call that evaluates `exp()` for logarithmic tables — while the lower
  bound at `:2517` uses `self%xv(1)`. Because `Table_Logarithmic_1D_Interpolate`
  (`:1542`) passes `log(x)` down, the upper test compares a logarithm with a
  linear bound and is effectively never true for tables with `x_max > 1`;
  `..._Interpolate_Gradient` (`:1561`) passes linear `x` instead. Use
  `self%xv(self%xCount)`; decide deliberately what upper-range behavior log
  tables should have (this is a latent correctness issue as much as a
  performance one).
* Linear `interpolator%interpolate` (`numerical/interpolation/1D.F90:661-707`)
  goes through the GSL FFI (`gsl_interp_eval_e`) although a cached pure-Fortran
  `locate` exists; `interpolatorLinearFactors` asserts twice. Evaluate linear
  interpolation in Fortran with GSL's exact formula. Micro × 10⁷–10⁹ calls.
* `sphericalPotentialNumerical` calls `Range_Pinned` before testing whether the
  table needs rebuilding (`spherical/_class.F90:599-609`). S.
* Non-integer `**` pairs where one `exp(log)` would do:
  `atomic/cross_sections/ionization/photo/Verner.F90:2050-2052,2078-2080`,
  `mass_distributions/spherical/Zhao1996.F90:363-368`. Not neutral. Micro.
* `hubbleParameterRateOfChange` (`cosmology/functions/matter_lambda.F90:638-650`)
  recomputes E(a) three times; `equalityEpoch*` take `**(1/3)` per call of a
  parameter-only quantity (`:820,855`). S.
* `sortIndex` and `searchArrayDouble` are FFI calls per comparison/lookup
  (`numerical/sort/_module.F90:70-90`, `numerical/search.F90:60-75`); fine for
  large arrays, slow for the tiny lists in §4.3.
* Verified fine: `rootFinder` bracket reuse and threadprivate stack; 1-D
  `interpolator` index-hint cache; all `table1D*` memoization; cosmology, σ(M),
  δ_c, growth, Farahi/Zhang-Hui, loss-cone, Kummer2018, Sersic tabulations on
  absolute lattices with HDF5 caching; `Values_Agree/Differ` elemental; the
  disk/spheroid component distributions already pooled.

## 5. Run time: tree pipeline, tasks, output

### 5.1 HDF5 output: probe + open + write + close per property per tree per output [I, verified]

`utility/IO/HDF5/_module.F90:4123-4190` (`writeDataset` with `appendTo`):
`hasDataset` (`:3944-3966`, implemented as `h5dopen_f` + `h5dclose_f`, not
`h5lexists_f`) → `openDataset` (`h5dopen_f`, `h5dget_create_plist_f`,
`h5pget_chunk_f`, two `h5eset_auto_f`) → extent + hyperslab + write → close.
Callers: `merger_trees/outputter/standard.F90:743-830` (one append per property
per tree per output) and `:333-337` (five scalar appends per tree per output),
all inside the process-wide `hdf5Access` lock, so every thread serializes here.
For 10⁵ trees × 10 outputs × 80 properties that is ~10⁸ metadata calls. Fix:
keep dataset handles open per output group; replace the probe with
`h5lexists_f`; buffer the per-tree scalars; track `referenceStart` as a counter
instead of re-reading the dataset extent (`:304-311`). M. Neutral. Macro for
many-small-tree runs.

### 5.2 Output buffer grows linearly (+1024) [I, verified]

`outputter/standard.F90:56,681-682,867-960`: the buffer is extended by
`standardBufferSizeIncrement=1024` and copied each time, and only dumped at the
end of a tree, so a tree with N output nodes costs O(N²·P/1024) element copies.
For a 10⁶-node N-body forest with 80 properties that is ~4×10¹⁰ copies per
output. Fix: dump when full (the append path already exists) or grow
geometrically. S. Neutral. Macro for large trees.

### 5.3 Equilibrium structure solver always iterates twice per RHS [I, verified]

`galactic_structure/radius_solver/equilibrium.F90:278`:
`do while (countIterations <= 1 .or. (fitMeasure > tolerance .and. ...))`,
hooked to `preDerivativeEvent` (`:171`), i.e. per RHS for every node with a
galaxy. Iteration 1 seeds from stored radii; iteration 2 always does the full
enclosed-mass (adiabatic contraction → numerical) and rotation-curve
evaluation per component even when the seed already satisfies the tolerance.
Fix: (a) key the converged solution on the exact structural input vector for
the node and skip the solve when unchanged (neutral) — the node evolver
already short-circuits identical `(time, y)` RHS calls at
`node_evolver/standard.F90:868-877`, so repeated-input calls do occur; (b)
make iteration 2 conditional on the seed failing the tolerance (changes
results at tolerance level). M. This is likely the largest per-RHS cost in
baryonic models with adiabatic contraction. The adiabatic root find *inside*
the enclosed-mass evaluation is [#575](https://github.com/galacticusorg/galacticus/issues/575)
(assessed there at ~0.4% of the milkyWay benchmark — small); the lever here is
avoiding the solve altogether when inputs are unchanged.

### 5.4 N-body reader: per-forest probes and opens under the global lock [I]

`merger_trees/construct/read/importer/galacticus.F90:1024-1398`: seven
`hasDataset` probes whose answers never change between forests, and up to ~40
`readDatasetStatic` calls each opening and closing its dataset, per forest,
under `hdf5Access`. Per-node `allocate(nodes(iNode)%reals(...))` at
`:1270-1291`. N-body runs have 10⁵–10⁷ mostly tiny forests, so metadata
dominates data volume. Fix: resolve presence once in `galacticusOpen`; keep
handles open; read blocks of consecutive small forests in one hyperslab; use a
single 2-D array for per-node values. M. Neutral. Macro for read-path runs.

### 5.5 Cole2000 build: root find plus two σ(M) evaluations per branch step [I]

`merger_trees/construct/builder/Cole2000/_class.F90:563-565,699-702`: every
branch step calls `criticalOverdensity_%timeOfCollapse` (a root find whose
memo never hits because δ_c and mass change every step) and two `rootVariance`
lookups. When `.not.timeParameterIsMassDependent` (the common ΛCDM case,
already computed at `:321`) the growth-factor ratio is `δ_c(t₀)/δ_c(t)` with no
inversion, and `time` is needed only when a branch actually occurs. Fix:
compute `time` lazily; cache the mass-only factors of `deltaWEarliestTime`.
M. Not neutral in general (different arithmetic path) unless restricted to
skipping unused evaluations. Plausibly 20–40% of Cole2000 build time, which
is ~20% of total in the reference model.

### 5.6 `nodeOperatorMulti` dispatches every phase to every operator [I]

`nodes/operators/multi.F90:246-425`: sixteen fan-out routines each loop over
all operators (20 in `quickTest.xml`) through polymorphic dispatch, although
most operators override only one or two phases (e.g. `Post`, `Scales`,
`Inactives`, `PostStep`: 0/20 override; `differentialEvolution`: 14/20). ~100
no-op virtual calls per node per step and ~20 per RHS, plus a virtual
`multiIsActive` that always returns `.true.`. Fix: at construction build one
compact list per phase of operators that override it (compare `c_funloc`
against the generated default, or add a generated `provides(method)` query).
S–M. Neutral. Micro–small (a few % of the derivative path).

### 5.7 Per-tree per-output metadata churn in the outputter [I]

`outputter/standard.F90:266-275` re-runs `propertiesCount`, `buffersAllocate`
and `propertyNamesEstablish` (`:961-1233`) for every tree at every output,
reallocating metadata containers and re-fetching `varying_string` names and
descriptions (~5–10 allocations per property) although nothing changes between
trees at a given output. Cache the layout per `indexOutput`. S–M. Neutral.

### 5.8 Other pipeline items

* Default `hdf5ChunkSize=1024` elements (8 kB) with a 1 MB chunk cache
  (`output/HDF5/open.F90:137-138,196-197`): for many datasets the cache holds
  fewer chunks than there are datasets being appended, so every dump thrashes.
  Raise the defaults (8k–64k elements, 32–64 MB cache) or size the cache from
  the property count. S. Data-neutral.
* Per-thread startup is fully serialized (`evolve_forests/_class.F90:645-649,709-719`):
  parameter re-parse + graph deep copy under two critical sections, so
  start-up scales as #threads. Matters on 64–128-thread nodes. M–L.
* Per-tree `!$omp critical` around the RNG deep copy in the tree constructors
  (`construct/build/_class.F90:290-296`, `construct/read/_class.F90:1010-1016,1879-1885`);
  the source object is already thread-private, so a plain `clone()` suffices. S.
* `hostTidalMassLoss` timestep walks the host ancestry per satellite per step
  (`evolve/timesteps/host_tidal_mass_loss.F90`); memoize per `uniqueID`,
  invalidate on promotion. S.
* `evolveForests` counts nodes with a full walk per tree (`:769-776`) although
  the read path already knows `nodeCounts`. S, micro.
* Extractor chain allocates per output node per extractor
  (`outputter/standard.F90:436-680`, `nodes/property_extractor/multi.F90:233-282`);
  preallocate per-thread scratch. M, small.
* Verified fine: forest scheduling (descending mass, FCFS ≈ LPT); output-time
  binary search; tree walkers iterative and allocation-free; reader index maps
  sorted + bisected; `darkMatterHaloScale`, cooling radius, `timeOfCollapse`
  memoization; no `displayMessage`/`allocate`/string compares inside any
  `*DifferentialEvolution` routine under `nodes/operators`.

## 6. Build system and Python tooling [all M unless noted]

Measured in this survey on a 4-core sandbox (Python 3.11, GNU make 4.3) with a
temporary `BUILDPATH` populated by running the real catalog scripts; no Fortran
compiler was available, so compile times are not included. Profiles are
reproducible with `python3 -m cProfile -o out.prof scripts/build/preprocess.py <src> <dst>`.

| Item | Measured |
|---|---|
| `preprocess.py` import overhead (warm pycache) | 65–85 ms (bare Python ≈ 10 ms); 99 `re.compile` at import |
| `preprocess.py` on an 86-line file with no directives | 0.17–0.20 s, ~0.10 s of it catalog XML loading |
| `preprocess.py` on a 74-line functionClass parent with 33 implementations | 3.7 s |
| 12 sampled `_class.F90` parents | mean 0.96 s → ×271 parents ≈ 260 CPU-s |
| 30 sampled ordinary files | mean 0.21 s → ×1715 files ≈ 360 CPU-s (~90% fixed overhead) |
| `source/utility/IO/HDF5/_module.F90` | 22.7 s, ~88% in `copy.deepcopy` |
| `source/objects/nodes/_class.F90` (component builder) | 33–36 s, 13 s of it system time (24k `brk` calls) |
| Extrapolated preprocessing per clean build | ≈ 680 CPU-s (~6 min wall at CI's `-j2`), > 500 s avoidable |
| `make -n Galacticus.exe` (nothing to do) | 6.3 s; with `-r`: 0.46 s |
| `moduleDependencies.py` cold / warm | 32 s CPU / 0.13 s |
| `sourceDigests.py` for `Galacticus.exe`, cold / warm | 34 s CPU (9.8 s wall) / 0.19 s |

### 6.1 `make` tries built-in RCS/SCCS pattern rules on every prerequisite

`Makefile:773` has `.SUFFIXES:`, which cancels suffix rules but **not** GNU
make's built-in *pattern* rules (`%: %,v`, `%: RCS/%,v`, `%: s.%`, …). During
makefile remaking, make searches implicit rules for every prerequisite of the
tree-scanning rules (lines 795, 860, 870, 898, 904, 908 list all ~2000 sources
and Python files), each probing the 555-directory `vpath`. `make -d` shows
14,311 "Trying pattern rule" attempts vs 8 with `-r`. Cost: ~6 s on every
make invocation, paid 2–3× per clean build because of makefile restarts, and
the floor of every incremental build. Fix: `MAKEFLAGS += --no-builtin-rules`
near the top of the Makefile (not `-R`). **Caveat:** the `LINK_METADATA`
recipe (`Makefile:697-701`) greps `MAKEFLAGS` for ` -j` patterns to decide on
MD5 locking; confirm that detection still works with the added flag. Effort S.

### 6.2 Quadratic string appends in the preprocessor

`python/Galacticus/Build/SourceTree/__init__.py:247` (`code_commented = code_commented + line`),
`:1008-1016` (`mappings +=`, `stripped +=`, `serialization +=`), and
`python/Galacticus/Build/Components/__init__.py:350-430` (23 `build['content'] +=`
sites). On `objects/nodes/_class.F90` (4.7 MB / 86k generated lines):
`_serialize` 23.7 s, `_comment_embedded` 17.7 s, `functions_serialize` 17.3 s
under profile; a micro-benchmark of the two loops on the real output goes from
18.5 s → 0.03 s and 19.0 s → 0.04 s with list + `''.join`. This file is on the
critical path of every clean build and is regenerated whenever any of 74
component-related files change. Effort S; output byte-identical.

### 6.3 Whole-tree `deepcopy` in Generics

`Process/Generics.py:142,281`: `copy.deepcopy(code_node)` on a node that
carries a `parent` back-pointer copies the entire tree (~3000 nodes) each time;
2768 calls → 33.9 M `deepcopy` frames on the HDF5 module (22.7 s wall). Fix:
build the temporary node from `content`/`source`/`line` directly, or detach
`parent`/`sibling` before copying. Effort S.

### 6.4 FunctionClass parents re-parse and re-process every implementation

`Process/FunctionClass/__init__.py:319-338` (`_load_and_sort_classes`) runs
`parse_file` **and** `process_tree` on each implementation file; each
implementation is also preprocessed by its own make rule, so every
implementation is processed twice per clean build, and again whenever the
parent is re-preprocessed (fan-out up to 163 files for
`nodes/property_extractor/_class`). ≈ 200 of the 260 CPU-s spent on parents.
Fix: cache the parsed+processed tree per file on disk (keyed by mtime +
`preprocessorSources.digest`), or run only the hook subset the downstream
builders need. Effort M–L; CI's determinism `cmp` on
`mass_distributions/_class.p.F90` guards the output.

### 6.5 Catalog XML re-parsed by 13 independent loaders

`directiveLocations.xml` (876 kB, 8,798 entries) and `stateStorables.xml`
(150 kB) are loaded by separate code in `MetaPropertyDatabase.py:25-30,68-69`
(uncached, every `process_tree`), `InputParameter.py:59-72` (parsed **and
discarded**, plus every `<exe>.d` read, on every invocation),
`InputParametersValidate.py:25-36`, `EventHooksStatic.py`, `FunctionsGlobal.py`,
`ObjectBuilder.py` (cached but loaded before checking relevance),
`Constructors.py`, `StateStorable.py`, `FunctionClass/__init__.py:129-141`,
`EventHooks.py`, `ComponentBuilder.py`, `SourceDigest.py`, `Parse/Directives.py`.
On a plain 86-line file ~0.10 s of the 0.17 s total is catalog parsing (10
parses) — half the per-file cost, ×1700 files ≈ 170 CPU-s. Fix: one shared
lazily-loaded catalog module, loaded only when a hook finds a relevant node;
have the catalog writers also emit a pickle sibling (~10× faster to load);
delete the dead reads. Effort S–M.

### 6.6 Fixed per-invocation overhead × ~2000 invocations

~70 ms import (`Fortran/Utils.py` compiles 99 regexes at import; `lxml.etree`
imported eagerly by `Parse/Directives.py`) + ~100 ms catalogs + 10 ms
interpreter ≈ 0.18 s × 1715 files ≈ 300 CPU-s against ~0.02 s of actual
parsing per file. Fix: §6.5 first; then a batch mode (`preprocess.py --batch`
or N files per recipe) that amortizes imports and catalogs; lazy-import `lxml`.
Effort S (lazy imports) to M (batch mode).

### 6.7 Parser hot spots

* `_parse_units` tries all nine unit-opener regexes on every logical line
  (`SourceTree/__init__.py:412`; patterns in `Fortran/Utils.py:95-178`):
  3.16 M `re.match` = 13.9 s on the component module; a keyword prefilter gives
  7.1 s → 0.10 s with an identical match set. S.
* `_function_class_names` / `function_class_entry`
  (`StateStorables.py:78-96`) rebuild a set / linear-scan per node:
  40k calls = 5.9 s on the component module. Memoize per `process_*` call. S.
* `walk_tree` is a recursive generator (`__init__.py:169-178`): 5.5 M frames
  for 1.26 M nodes; ~35 hooks each walk the full tree. Iterative walker plus a
  once-per-tree directive index. S–M.
* Uncompiled regexes per line in `FortranUtils.py:48,70-79`,
  `Parse/Declarations.py:38`, `Parse/Directives.py:371-374`. S.

### 6.8 Dependency scripts

* `moduleDependencies.py:113-160` does a full `parse_file` of every
  implementation just to find `type … extends`; 32 s CPU cold, invalidated on
  any source add/remove. A regex scan (as `_scan_file` already does for
  modules) or the data already in `stateStorables.blob` would do. S–M.
* `sourceDigests.py:247` re-parses every functionClass instance at **link
  time**, per executable: 34 s CPU / 9.8 s wall on the serial tail before the
  LTO link; the blob is per target, so `make all` (156 executables) redoes it
  per executable. Share the per-file `class_dependencies` blob. M.
* `enumerateOpenMPCriticalSections.py`, `includeDependencies.py`,
  `findExecutables.py` have no incremental cache and re-run (with a make
  restart) on any source edit (~3.5 s + 6 s from §6.1 before any compile).
  Give them per-file blobs like the others. S each.

### 6.9 CI

Seven Linux jobs each do a clean build at `-j2` (`cicd.yml:125,180,223,255,286,323,363`)
because the cap protects the -O3 `_class` compiles; the ~680 CPU-s of Python
preprocessing (≈100 MB RSS) runs at that parallelism too. Fix: two-phase build
(`make -j$(nproc)` on a phony `preprocess` target listing every `.p.F90`, then
`make -j2 Galacticus.exe`), and `actions/cache` for `work/build/*.p.F90`,
`.up`, `.blob` and catalogs keyed on `preprocessorSources.digest` + source
hash (the only-if-changed semantics make restored files safe). S–M.

## 7. Compiler and link flags [I]

* **No `-march`/`-mtune`**: `Makefile:97-108` compiles for baseline x86-64
  (SSE2). Galacticus is pointer-chasing heavy so gains are modest, but the
  math-library share (~10%) and the vectorizable table/serialization loops
  would benefit from `-march=x86-64-v3` (AVX2/FMA) in the static release
  binary (documented as a `GALACTICUS_FCFLAGS` option; `x86-64-v3` is safe on
  any CPU from 2013 on).
* **No profile-guided optimization**: the code is dominated by polymorphic
  dispatch and branchy generated code, exactly where `-fprofile-use` helps
  most (typically 5–15%). The `profile.yml` workflow already runs
  representative models and could produce the `.gcda` set; PGO composes with
  LTO. Worth one experiment on the reference DMO and a baryonic model.
* **`-ffpe-trap=invalid,zero,overflow`** in release flags implies
  `-ftrapping-math` semantics (GCC default anyway), which inhibits some
  if-conversion and vectorization of conditional FP code. Removing the traps
  would be a debugging-vs-speed policy change; at minimum, measure the delta
  once so the cost is known.
* `-fstack-arrays` (part of `-Ofast`) is not set; automatic arrays above
  `-fmax-stack-var-size` go to the heap, which matters for the many small
  temporaries in generated code. Cheap experiment.
* `-DHAVE_INLINE -DGSL_RANGE_CHECK_OFF` for the vendored GSL sources (§4.5).

## 8. Suggested order of attack

1. **Build (one afternoon, all bitwise-identical output):** `MAKEFLAGS += -r`
   (§6.1); list/join in the three serializers (§6.2); detach before
   `deepcopy` (§6.3); shared catalog loader and dead-read removal (§6.5);
   opener prefilter and regex hoisting (§6.7). Expected: ~6 s off every make,
   ~55 s off the clean-build critical path, and roughly half of the per-file
   fixed cost.
2. **Run time, neutral, small:** pool the remaining `get` factories starting
   with `dark_matter_only.F90` and `adiabatic_Gnedin2004.F90` (§2.2); hoist
   the per-call `rootFinder`/`integrator` objects (§4.1–4.3); `Display_Lock`
   ordering (§3.2); `incrementSerialized` for the three object types (§3.1);
   `h5lexists_f`, geometric buffer growth and cached output layout (§5.1,
   §5.2, §5.7); GSL inline macros (§4.5). Validate each with `h5diff` against
   `master` on `quickTest.xml` and the DMO reference.
3. **Run time, neutral, medium:** keep dark-matter-only distributions across
   `Calculations_Reset` when inputs are unchanged (§2.1); input-keyed skip of
   the equilibrium solver (§5.3a); open-handle HDF5 output and reader (§5.1,
   §5.4); Jacobian `dydt_in` pass-through (§4.4); operator-phase dispatch lists
   (§5.6).
4. **Run time, not neutral, measured on the benchmark grid of the earlier
   analysis:** heated-profile inversion per the existing plans on #701/#913
   (§2.3); Cole2000 lazy
   `timeOfCollapse` (§5.5); tidal-tensor scale fraction as a parameter (§2.4);
   `-march`/PGO/`-fstack-arrays` experiments (§7).
5. **Larger refactors:** implementation-tree cache for functionClass parents
   (§6.4); per-node component-array footprint (§3.3); parallel per-thread
   start-up (§5.8).

## 9. Caveats

* No model was run and no Fortran was compiled in this survey; every run-time
  estimate is from code reading, cross-checked against the earlier `perf`
  profiles where the symbols overlap. Items marked as candidates for large
  gains (§2.3, §5.3, §5.1) should be confirmed with a `perf` profile of a
  representative model before implementation.
* `master` has been actively optimized recently (object pools, pinned
  tabulations, cache files, lockless event hooks); several items the survey
  set out to check were already fixed and are listed under "verified fine" so
  they are not re-done.
* The logarithmic-table upper-range check (§4.6) is a behavioral
  inconsistency as well as a micro-cost and deserves a deliberate decision
  rather than a silent fix.
