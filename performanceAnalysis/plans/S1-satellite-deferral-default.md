# Avoid satellite micro-stepping by default (`fractionTimestepSatelliteMinimum`)

*Filed as [galacticusorg/galacticus#1431](https://github.com/galacticusorg/galacticus/issues/1431).*

## Summary

`mergerTreeEvolverStandard` supports deferring a satellite's evolution when the
step it is allowed to take is only a small fraction of its host-relative
timestep (`fractionTimestepSatelliteMinimum`). The reference configuration
files set this to `0.75`, but the **code default is `0.0`** (no deferral), and
any configuration not built on the reference includes runs without it. The
measured cost of running with the default is large and grows with resolution:
8.3x at M_host/m_res = 10^5, 15.2x at 3.16x10^5, and a run-time scaling of
(1/m_res)^1.57 in place of the reference configuration's (1/m_res)^0.98-1.16
(see the analysis under `performanceAnalysis/` on branch
`claude/galacticus-performance-analysis-giadol`). This issue proposes making
deferral the default behavior, in phases.

## Mechanism (why the default is expensive)

Tree evolution proceeds by repeated full-tree walks
(`source/merger_trees/evolver/standard.F90:359`). Merger events synchronize the
merge target to within `timeOffsetMaximumAbsolute` (10 Myr), fragmenting the
host's natural ~0.1/H stepping into event-sized pieces, each costing a full
walk; walk counts grow superlinearly with M/m_res (measured 8 at 10^3 up to
10,468 at 3.16x10^6). During these "stalled" walks, satellites are permitted
steps of only the host's small advance. With deferral off, every such visit
becomes a micro ODE call carrying the full fixed per-call overhead
(~130 us: node-operator pre/post hooks, three serialization passes, scale
computation, a freshly allocated GSL driver, initial-step-size search) for a
~10 Myr or smaller step. With deferral on (0.75), those visits cost ~1 us and
the walk machinery stays subdominant until M/m_res ~ 10^7. The variant
experiment (identical trees, identical walk counts, `0.75` vs `0.0` only)
isolates this cleanly: 17.1 s vs 142.1 s evolve time per tree at m_res=10^8,
M_host=10^13 M_sun.

## Proposed change (phased)

Phase 1 - shipped configurations (no behavior change for users' own files):

- Set `<fractionTimestepSatelliteMinimum value="0.75"/>` in the shipped
  configurations that currently omit it, notably `parameters/quickTest.xml`
  and `parameters/parametersProfile.xml` (the reference `evolution*.xml`
  includes and `baryonicPhysicsConstrained.xml` already set it). Currently 42
  of 235 `testSuite/parameters/*.xml` set the parameter.

Phase 2 - visibility:

- Extend the parameter's `<description>` (in
  `source/merger_trees/evolver/standard.F90`, the `inputParameter` block for
  `fractionTimestepSatelliteMinimum` at ~l.199-206) to state the measured
  performance consequences of `0`, and reference the scaling analysis.
- Optionally: emit a one-time `displayMessage` at `verbosityLevelWarn` from
  the evolver constructor when the value is `0`, noting that satellite
  micro-stepping may dominate run time at high resolution.

Phase 3 - default change:

- Change `defaultValue` from `0.0d0` to `0.75d0`. This is a
  behavior-affecting default: models that relied on the old default will
  produce slightly different results (identical physics; different ODE step
  boundaries, so differences are at the integration-tolerance level).
  Announce in release notes; refresh any validation baselines that shift
  outside tolerance (the DMO subhalo validation configuration already runs
  with `0.75`, so its baselines are unaffected).

## Correctness considerations

- Deferral delays a satellite's evolution by at most one host-relative step
  (<= min(0.1/H, 1 Gyr) x (1 - f)); the satellite then takes one larger step
  over the same interval. All rate integrals span the same total time; the
  approximation is the same one already accepted in the reference and
  validation models.
- Merger synchronization is unaffected: the satellite timestep class
  (`source/merger_trees/evolve/timesteps/satellite.F90`) still caps a merging
  satellite at its merge time, and the deferral logic applies only to the
  "satellite in host" criterion (`evolver/standard.F90:925-931`).
- `backtrackToSatellites` measured neutral (17.0 vs 17.1 s) - no change
  proposed there.

## Validation and benchmarking

- Rerun the resolution grid with the scripts under `performanceAnalysis/scripts`
  (`run_one.sh def_m<r> <r> 4` with and without the new default): acceptance is
  the defaults-series evolve-time slope returning to ~1.0 over
  M/m_res = 10^4-3x10^5.
- `testSuite/validate-darkMatterOnlySubhalos.py` and one baryonic validation
  unchanged within their tolerances.
- Confirm no change at all for configurations that already set the parameter.

## Effort

Small: parameter default + description text + two shipped configuration files
+ release-note entry. The main cost is the validation-baseline audit for
Phase 3.
