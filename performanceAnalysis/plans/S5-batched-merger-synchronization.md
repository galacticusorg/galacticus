# Batch merger synchronization to reduce event-driven tree-walk fragmentation

## Summary

Subhalo mergers synchronize the merging satellite with its merge target to
within a tolerance window (`timeOffsetMaximumAbsolute` = 10 Myr /
`timeOffsetMaximumRelative` = 10^-3 in the satellite timestep class). The
evolver, however, caps the target's step at *exactly* the earliest mergee's
merge time, so each merger forces its own host step - and each host step in
the stalled regime costs a full O(live nodes) tree walk. Walk counts grow
superlinearly with resolution (measured 8 at M/m_res = 10^3 to 10,468 at
3.16x10^6; analysis under `performanceAnalysis/` on branch
`claude/galacticus-performance-analysis-giadol`). This issue proposes letting
the target *overshoot* each merge time by up to the already-accepted
synchronization tolerance, so that all mergers falling within one window
complete against a single host step instead of one step (and walk) per
merger.

## Current behavior

1. A merger triggers as a mid-ODE interrupt when a satellite falls inside
   0.01 r_vir; `timeOfMerging` is set to the current time
   (`source/nodes/operators/physics/satellite_merging/radius_trigger.F90:222`).
2. The satellite's own step is capped at `timeOfMerging`, and it may not
   complete the merger until the target has reached
   `timeOfMerging - delta`, where
   `delta = min(timeOffsetMaximumAbsolute, timeOffsetMaximumRelative x t)`
   (`source/merger_trees/evolve/timesteps/satellite.F90:178-207`). There is no
   *upper* bound on the target's time in that trigger condition - the merger
   fires whenever the target is at or beyond `timeOfMerging - delta`.
3. The upper bound is imposed by the evolver instead, in two places in
   `standardTimeEvolveTo` (`source/merger_trees/evolver/standard.F90`):
   - the mergee scan caps the target at `max(timeOfMerging, timeNode)`
     (l.990-1004); note the existing handling of a target *already past* a
     mergee's merge time - it is frozen at its current time until the mergee
     completes, i.e. "target ahead of merge time" is already a legal,
     handled state;
   - the hosted-satellite scan caps a host at `min(basic%time())` over its
     satellites (l.963-984), and a satellite parked at its merge time holds
     the host there just the same.

So with k mergers inside one 10 Myr window the target advances k times, each
advance requiring at least one additional full tree walk.

## Design

Introduce an opt-in overshoot in the evolver,
`timeOffsetMergeeMaximumAbsolute` / `...Relative` (defaults `0` - current
behavior preserved bit-for-bit), applied as
`overshoot = min(absolute, relative x t)`:

1. **Mergee scan** (l.990-1004): cap becomes
   `max(timeOfMerging + overshoot, timeNode)`. The existing freeze logic then
   holds the target at its arrival time until the parked mergees complete -
   a state the code already handles today.
2. **Hosted-satellite scan** (l.963-984): a satellite parked for merging
   (`satellite%timeOfMerging() <= basic%time()`) caps its host at
   `basic%time() + overshoot` instead of `basic%time()`. Reading the
   satellite component in the scan is only needed on the binding path (when
   `timeSatellite < evolveToTime`), and this is the single point that must be
   coordinated with the earliest-satellite-time cache proposed in the
   satellite-bookkeeping issue (cache the *effective* limit).
3. **No change** to `timesteps/satellite.F90`: its trigger condition already
   tolerates a target ahead of the merge time; each parked mergee completes
   its merger via its own end-of-step task as soon as it is visited with the
   target inside the window. With the target landed at
   `t_earliest + overshoot`, every mergee with
   `timeOfMerging in [t_earliest, t_earliest + overshoot]` completes without
   any further target advance - the batching.
4. Recommended wiring: set the overshoot equal to the satellite timestep
   class's `timeOffsetMaximumAbsolute/Relative` (0.01 Gyr / 10^-3) in the
   reference configurations, keeping the total asynchrony between a merging
   pair within the tolerance the model already accepts (the offset changes
   sign: the target is up to one window *ahead* at the merger instead of up
   to one window behind).

Parameter duplication note: the window values then exist in two classes. An
alternative is a shared parameter object, but given both are plain numbers
with the same meaning ("acceptable merger asynchrony"), documenting the
recommended pairing in both descriptions seems sufficient; opinions welcome.

## Expected gain and its limits

The gain is bounded by the share of walks that merger synchronization forces:
batching collapses runs of mergers separated by less than the window into
one target step. Destruction events (the other event class) do not cap the
host and are unaffected. This is the least-certain of the four proposed
optimizations - it should be implemented behind its default-off parameter and
judged on the benchmark: acceptance is a measurable walk-count and
evolve-time reduction at M/m_res >= 10^6 (10^12 M_sun host at m_res <= 10^6
is the cheap probe: 3,325 walks, 178.7 s/tree at 10^6 on the analysis
machine), with unchanged validation statistics; if the gain is <5-10% there,
close as not-worth-the-complexity.

## Correctness considerations

- The overshoot never exceeds the merger-asynchrony tolerance the model
  already accepts by construction; convergence: overshoot -> 0 recovers the
  current behavior exactly.
- Mergers within one window are processed as effectively simultaneous - which
  is precisely the semantics the existing window already assigns them
  (targets today may lag each merge time by up to delta).
- Galaxy-merger physics (`galaxiesMerge` node operators, satellite removal,
  sub-subhalo promotion) sees a target at most one window ahead of the merge
  time; today it sees one at most one window behind. Any operator that
  assumed exact time equality would already be wrong today.
- Deadlock safety: parked mergees still advance to their merge times
  independently (their caps are unchanged), so the trigger condition is
  reachable; the target frozen at `t + overshoot` is released as mergees
  complete (tree mutations reset the deadlock status).

## Validation and benchmarking

- Default-off: bitwise-identical outputs to master with the parameter unset.
- With overshoot = 0.01 Gyr: `validate-darkMatterOnlySubhalos.py` statistics
  unchanged within tolerance; subhalo mass function / radial distribution
  compared explicitly between overshoot = 0 and 0.01 at m_res = 1e8.
- Benchmark grid as in the other issues; report walk counts (verbosity=warn
  logs) alongside evolve times.

## Effort

Small-to-moderate: two scan-site changes plus two parameters in
`mergerTreeEvolverStandard`, documentation, and the benchmark/validation
runs. Interacts constructively with the satellite-bookkeeping issue (shared
effective-limit accessor).
