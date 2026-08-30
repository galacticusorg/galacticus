# Avoid per-call ODE driver allocation and reduce fixed node-evolver overhead

*Filed as [galacticusorg/galacticus#1433](https://github.com/galacticusorg/galacticus/issues/1433).*

## Summary

Every call to `mergerTreeNodeEvolverStandard::evolve` pays a large fixed cost
regardless of step size: a freshly allocated GSL ODE driver (system + stepper
+ evolve + control objects, with internal arrays malloc'd and freed each
call), three serialization-offset passes, value/scale serialization with saved
copies, and an initial-step-size search. Measured: ~130 us per ODE step,
constant across step sizes; perf shows malloc/free at ~4-6% self-time plus
kernel page-zeroing rising to 2.5% at high resolution (analysis:
`performanceAnalysis/` on branch
`claude/galacticus-performance-analysis-giadol`). This constant multiplies
*every* term of the run-time scaling - including the forced micro-steps that
dominate high-resolution subhalo runs when satellite deferral is off. This
issue proposes reusing driver workspaces across calls, plus smaller
reductions of the per-call scaffolding.

## Current per-call allocation chain

`node_evolver/standard.F90:585-643` constructs an `odeSolver` per evolve call
(and again per retry/fallback). The constructor
(`source/numerical/ODE_solver/solver.F90:335-426`) performs:

- `allocate(self%gsl_odeiv2_step_type)`, `allocate(self%system)`,
  `allocate(self%driver)`, `allocate(is_non_negative(...))`;
- `gsl_odeiv2_system_init` - mallocs a system struct
  (`source/interface/GSL/C/odeiv2.c:73-87`);
- `gsl_odeiv2_driver_alloc_scaled2_new`
  (`source/external/gslODEInitVal2/driver2.c:104`) - mallocs the driver,
  stepper state (several dimension-sized arrays), evolve workspace, control
  object, and *copies* of `scale_abs` and `is_non_negative`
  (`cscal2.c:217-`).

The destructor frees it all. That is ~10+ malloc/free pairs per node-evolve
call, at ~1.3M+ calls per 10^13 M_sun tree at m_res=10^7.

Because `gslODEInitVal2` is vendored in-repo (it already carries the custom
`scaled2` control and the `msbdfactive` stepper), a reuse API can be added
without any upstream-GSL constraint.

## Design

### A. Driver reuse across calls (expected main win)

1. New C entry point in `source/external/gslODEInitVal2/driver2.c`:

       int gsl_odeiv2_driver2_reuse(gsl_odeiv2_driver *d,
                                    const double hstart,
                                    const double scale_abs[],
                                    const int is_non_negative[]);

   For a driver of unchanged dimension and tolerances: call the existing
   `gsl_odeiv2_driver_reset` (resets stepper and evolve state, sets
   `h = hstart`), `memcpy` the new `scale_abs` / `is_non_negative` into the
   `sc2_control_state_t` arrays, and reinitialize the error vector
   (`gsl_odeiv2_driver_init_errors`). No allocation.
2. Fortran side (`solver.F90`): give `odeSolver` a `reset` method wrapping the
   above; record `dim`, stepper type, and tolerances so a mismatched reuse
   request is refused (assertion).
3. `mergerTreeNodeEvolverStandard` holds a small per-instance (the evolver is
   already per-thread via deep copies) cache of `odeSolver` objects keyed by
   (dimension, algorithm, jacobian mode). In practice only a handful of
   distinct dimensions occur per model (satellite vs isolated-node property
   counts), so a linear-searched array of ~8 entries suffices. The evolve
   call replaces construct-with-`conditionalCall` (l.613-636) by
   lookup + `reset`; the callback targets (`standardODEs`,
   `standardPostStepProcessing`, etc.) are module procedures, identical for
   every call, so the cached system object remains valid.
4. Tolerances vary only via constructor parameters (fixed per evolver);
   `hStart` varies per call and is handled by `reset`. The non-Jacobian
   fallback path constructs with a different algorithm - covered by the cache
   key.

Results should be bitwise-identical: a reset driver is in the same state a
freshly allocated one would be (same control parameters, cleared stepper
state, `h = hstart`). Assert this in validation.

### B. Serialization-offset memoization (phase 2, measurement-gated)

`node%serializationOffsets` runs three times per call
(`node_evolver/standard.F90:455,518-519`), walking all components. The
`propertyTypeAll` pass depends only on the node's component layout; the
active/inactive passes depend on flags that operators set per call, so only
the first is safely memoizable by (node uniqueID, layout revision). Requires
adding a layout-revision counter to the node serialization machinery -
inspect `python/Galacticus/Build/Components/TreeNodes/Serialization.py`
before committing to this; gate on a profile showing the offset passes still
matter after A.

### C. Step-size reuse in shipped configurations

`reuseODEStepSize=true` measured a free ~4% (skips the initial-step-size
search by seeding from the node's previous step). Not bitwise-neutral
(stepping changes within tolerance). Propose enabling in the reference
configurations; a default change can ride the same release as other
behavior-affecting defaults.

### D. (Future, separate design) `Calculations_Reset` fan-out

The `calculationReset` event hook fires for every subscribed object on every
evolve call (`node_evolver/calculations_reset.F90`). A generation-counter
scheme (objects compare a global epoch instead of being called) would turn
this O(#subscribers) into O(1), but touches the memoization pattern used
across many classes - out of scope here, noted for completeness.

## Validation and benchmarking

- A: bitwise-identical outputs vs master on `quickTest.xml` and the DMO
  subhalo reference configuration (fixed seed); any drift is a bug.
- Run under valgrind massif / heaptrack on a short model to confirm the
  allocation-rate drop.
- Benchmark (`performanceAnalysis/scripts/run_one.sh`) at m_res = 1e8 and
  3.16e7 (10^13 M_sun): expect ~5-15% evolve-time reduction in the reference
  configuration (malloc/free + page-zeroing + init shares), and substantially
  more in configurations that micro-step (deferral off), where per-call
  overhead dominates. Re-profile to confirm malloc/free/page-zeroing shares
  drop.
- Full testSuite pass.

## Effort

Moderate: ~40 lines of C in the vendored library, a `reset` method and
constructor bookkeeping in `solver.F90`, cache plumbing in
`node_evolver/standard.F90`. B and D are follow-ups behind measurement gates.
