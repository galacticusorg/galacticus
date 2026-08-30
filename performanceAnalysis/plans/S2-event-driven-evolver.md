# Event-driven merger-tree evolver (priority queue scheduling)

*Filed as [galacticusorg/galacticus#1435](https://github.com/galacticusorg/galacticus/issues/1435).*

## Summary

The standard evolver schedules work by *sweeping*: it repeatedly walks every
node of a tree, attempting to evolve each, until a complete walk evolves
nothing (`source/merger_trees/evolver/standard.F90:359`). At high resolution
the sweep count is driven by event synchronization and grows superlinearly
(measured: 8 walks per tree at M_host/m_res = 10^3, 195 at 10^5, 10,468 at
3.16x10^6), while the cost of each sweep grows with the live node count, so
total scheduling work - node visits - scales as (1/m_res)^2.06 (measured
2.4k visits/tree at M/m = 10^4, 245k at 10^5, 1.4M at 3.16x10^5). In the
stalled regime a sweep visits ~1,250 live nodes to advance ~20, sometimes 0.
This term is the structural cause of the run-time scaling steepening from
1/m_res toward 1/m_res^2 (through the observed ~1/m_res^1.5 crossover).

This issue proposes an event-driven evolver: a priority queue of evolvable
nodes ordered by current time, plus per-node waiter lists built from
dependency information the evolver already computes (`nodeLock` - which node
blocks each blocked node - today used only for deadlock reporting). Popping a
node evolves it to its allowed limit and wakes exactly the nodes waiting on
it. Scheduling cost becomes O((evolutions + wakes) log N) - proportional to
work actually performed - and the stalled-sweep regime disappears. The
scheduler changes; the physics stack (timestep classes, node operators, node
evolver, interrupts, end-of-step tasks) is untouched.

Analysis: `performanceAnalysis/` on branch
`claude/galacticus-performance-analysis-giadol`.

## Why this is possible with existing machinery

- `standardTimeEvolveTo` already identifies, for every limit it imposes, the
  node responsible, via its optional `nodeLock`/`lockType` arguments
  (`evolver/standard.F90:762-1035`), and every timestep class propagates the
  same through `mergerTreeEvolveTimestepClass::timeEvolveTo`. Inverting this
  relation ("who blocks me" -> "whom do I wake when I advance") is the entire
  scheduling insight.
- Every blocking constraint is relaxed only by an advance or mutation of its
  reported lock node. Enumerating the lock types: "sibling" (sibling advance
  or its merger), "hosted satellite" (satellite advance/removal), "satellite
  in host" (host advance), "mergee" and "satellite (host)" merger
  synchronization (target advance / mergee merger), "promotion" (a target,
  not a wait), tree/node events (fire when the node itself reaches the event
  time). So registering a blocked node against its single reported lock node,
  and re-computing on wake, misses no wakeups - see also the safety valve
  below.
- The evolver is a `functionClass` (`merger_trees/evolver/_class.F90`), and
  `mergerTreeEvolverThreaded` already demonstrates the extension pattern:
  `type, extends(mergerTreeEvolverStandard)`, overriding `evolve` while
  reusing `nodeIsEvolvable`, `timeEvolveTo`, and the deadlock machinery
  (`evolver/threaded.F90:251-`). The event-driven evolver follows the same
  pattern: `<mergerTreeEvolver value="eventDriven"/>`, fully opt-in.

## Design

### Data structures

- **Ready heap**: array-backed binary min-heap of (time, insertion counter,
  node pointer), keyed on the node's current `basic%time()`; the insertion
  counter breaks ties deterministically (reproducible runs). No heap utility
  exists yet in `source/utility/` - add a small one (~120 lines) local to the
  evolver or as a `Priority_Queues` utility module.
- **Waiter lists**: per-node intrusive singly-linked list. New `treeNode`
  members via the generator (`python/Galacticus/Build/Components/TreeNodes/
  __init__.py`, `_BASE_DATA_CONTENT`): `waiterFirst`, `waiterNext`
  (pointers), plus an `indexHeap` integer for O(log N) removal of destroyed
  nodes. (~24 bytes/node against ~3 kB/node.)

### The loop

Seed: one `mergerTreeWalkerAllNodes` pass at `evolve()` entry (spanning the
forest's linked trees) enqueues every node passing `nodeIsEvolvable` - the
cost of a single sweep, paid once per evolve call instead of thousands of
times. Then repeat until the heap is empty:

1. Pop the earliest node. Compute `timeEvolveTo` (the inherited routine,
   always requesting `nodeLock`).
2. Zero step (evolveToTime == current time): unlink from heap and append to
   `nodeLock`'s waiter list. (Zero step with no lock node - defensively -
   goes to a global retry list woken by any mutation.)
3. Positive step: call `mergerTreeNodeEvolver_%evolve` exactly as the
   standard evolver's per-node body does (same interrupt loop, same
   end-of-step `timestepTask_` handling, same promotion/merger logic,
   `evolver/standard.F90:424-507` transplanted). Then:
   - **advance**: wake the node's waiters (move list to heap), re-enqueue the
     node at its new time (unless at `timeEnd` / end of tree - terminal,
     dropped; rediscovered by the next evolve call's seeding walk);
   - **promotion**: the surviving node inherits and wakes both waiter lists
     and is re-enqueued (this is also how a parent becomes evolvable: its
     last child promotes into it);
   - **merger / destruction** (node destroyed inside the task): the popped
     node is not in the heap while being processed, so no dangling entry;
     wake its waiters, its pre-mutation host/parent's waiters, and its
     merge target's waiters, plus the global retry list (mutations are O(S)
     total, so coarse waking is cheap and removes missed-wakeup risk).
4. Empty heap with non-empty waiter lists -> **safety valve** (below); empty
   heap with empty waiters -> done (`treeDidEvolve` per the usual
   bookkeeping).

Tree events (`tree%event`) already cap node steps inside `timeEvolveTo`;
perform them (the inherited `standardTreeEventsPerform`) whenever the heap
drains, matching the standard evolver's end-of-sweep semantics.
`systemClockMaximum`, `status` passthrough, `suspendTree`
(`deadlockStatusIsSuspendable`), progress reporting
(`metaTreeProcessingTime`), and the evolve profiler hooks all keep their
existing behavior; `deadlockReporting` mode simply runs the inherited
reporting sweep.

### The safety valve (de-risking)

If the heap empties while waiters remain, fall back to one inherited standard
sweep over the waiters (recompute `timeEvolveTo` for each; re-enqueue any
that can advance) and count the event. A missed wakeup - the failure mode of
any event-driven scheduler - thus degrades to one extra sweep and a loud log
counter instead of a hang; a sweep that frees nothing is a true deadlock and
flows into the existing deadlock report. During burn-in the valve counter
measures how complete the wake rules are; the target is zero valve
activations on the full testSuite.

## What this does *not* fix (complementarity with #1431-#1434)

- The host's O(S) satellite/mergee scans inside `timeEvolveTo` still run once
  per *host consideration*; with ~S host steps that is still O(S^2) -
  #1432's caches remain necessary and complementary.
- Forced micro-steps from satellite non-deferral remain (they are
  schedule-independent ODE calls) - #1431 remains the first-order lever.
- Merger-synchronization still fragments the host's stepping (fewer full
  sweeps, but each fragment still costs heap traffic and a host ODE call) -
  #1434 still helps.
- The threaded evolver's intra-tree parallelism is orthogonal; a parallel
  event-driven scheduler (workers popping the heap with dependency guards) is
  a possible future unification, deliberately out of scope.

## Correctness considerations

- **Same rules, different order.** Each scheduling decision uses the same
  `timeEvolveTo` on the same current state; what changes is the interleaving
  of decisions. Node ODE results depend (at integration-tolerance level) on
  their step boundaries, so results are *not* bitwise-identical to the
  standard evolver - the same situation as any walk-order change. Validation
  is therefore statistical (below), unlike the pure data-structure issues.
- **Progress**: every pop either advances time, mutates the tree, or moves a
  node to a waiter list; nodes enter waiter lists only with a recorded waker.
  Earliest-first ordering means a popped node is never blocked by a node
  behind it in time.
- **Determinism**: heap with tie-break counter + single-threaded per-tree
  processing gives run-to-run reproducibility, as now.
- **Destroyed nodes**: only the in-flight (popped) node is ever destroyed;
  `indexHeap` covers any future case where structural code must remove a
  queued node.

## Phasing

1. Implement `mergerTreeEvolverEventDriven` (opt-in), sharing the per-node
   evolve body with the standard evolver by extracting it into a (private)
   method of `mergerTreeEvolverStandard` - the threaded evolver may adopt the
   same extraction, reducing its duplication too.
2. Burn-in: full testSuite with the valve counter required to read zero;
   statistical validation (below).
3. Benchmark; if clean, consider making it the default in a later release
   (the standard evolver remains available and is the reference
   implementation).

## Validation and benchmarking

- `validate-darkMatterOnlySubhalos.py` and one baryonic validation:
  statistics unchanged within tolerance vs the standard evolver; explicit
  comparison of subhalo mass function / radial distribution at m_res = 1e8.
- Tolerance-convergence check: halving `odeToleranceAbsolute/Relative` must
  shrink the standard-vs-eventDriven differences, demonstrating they are
  integration-tolerance-level, not systematic.
- Full testSuite under both evolvers; valve counter = 0.
- Benchmark (`performanceAnalysis/scripts/run_one.sh`) on the measured grid:
  (10^13 M_sun; m_res = 1e8, 3.16e7, 1e7) and (10^12 M_sun; m_res = 1e6,
  3.16e5). Acceptance: the walk/visit term is gone - evolve-time slope ~1.0
  through M/m = 3.16e6 (vs measured 1.16), with the residual excess
  attributable to the O(S) scans (#1432). Also benchmark with deferral OFF,
  where sweep overhead is largest today.
- Report scheduling statistics (pops, wakes, valve activations) at
  verbosity=warn, replacing the walk-count reports.

## Effort

The largest of the five proposals: a new ~700-900 line evolver class plus the
shared-body extraction from the standard evolver, two-to-three generated
`treeNode` members, and a small heap. High-confidence design (the dependency
information and the extension pattern both already exist), with the safety
valve turning residual wake-rule bugs into measurable performance blips
rather than hangs.
