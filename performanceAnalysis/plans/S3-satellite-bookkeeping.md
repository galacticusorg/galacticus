# O(1) satellite-list bookkeeping in tree evolution

## Summary

Several hot paths in tree evolution perform O(S) walks over a host's
satellite linked list, where S is the subhalo count (S ~ M_host/m_res;
measured 1,204 per 10^13 M_sun tree at m_res=10^8, z=0). Because these walks
are themselves executed O(S) or O(walks) times, they contribute quadratic
terms to the run-time scaling with resolution. perf profiling shows their
symbols (`basicStandardTimeGet` in scan loops 4.9%, mergee-scan
`timeOfMergingGet` 2.5%, `Satellite_Move_To_New_Host`+`standardPromote` 2.0%,
`simpleTimeEvolveTo` 1.0%) rising from <0.5% each at M/m_res <= 10^5 to a
combined ~10% at 10^6 - matching the measured excess of evolve time over
linear scaling there (analysis: `performanceAnalysis/` on branch
`claude/galacticus-performance-analysis-giadol`). This issue proposes making
these operations O(1): a satellite-list tail pointer, an earliest-satellite
time cache, an earliest-mergee cache, and O(1) satellite-list transfer on
promotion.

## Current O(S) operations and their sites

1. **Append to a host's satellite list** via `lastSatellite()` (linear walk;
   `source/objects/nodes/_class.F90:703-716`). Callers:
   - `source/merger_trees/node_merger/single_level_hierarchy.F90:102`
   - `source/merger_trees/node_merger/multi_level_hierarchy.F90:102`
   - `source/satellites/promotion.F90` (`Satellite_Move_To_New_Host`)
   - `source/merger_trees/node_evolver/standard.F90:1293` (`standardPromote`)
   Every satellite accretion, sub-subhalo move, and trunk promotion pays this.
2. **Hosted-satellite minimum-time scan** in `standardTimeEvolveTo`
   (`source/merger_trees/evolver/standard.F90:963-984`): O(S) per host
   consideration, at least twice per tree walk.
3. **Mergee scan** (`evolver/standard.F90:990-1004`, and the satellite branch
   at l.899-906): O(#mergees) per consideration. (The stray development
   comment at l.989 - "need to also figure out how to speed up satellite
   evolve checks" - refers to exactly this.)
4. **Satellite re-parenting on promotion**
   (`node_evolver/standard.F90:1284-1305`): every trunk promotion re-points
   `parent` on the full transferred satellite list - O(S) per promotion, with
   promotions ~ 1/m_res.
5. **Removal from a host's list** (`objects/nodes/_class.F90`,
   `Tree_Node_Remove_From_Host`): O(position); called on every merger and
   destruction event.

## Design

The `treeNode` type is generated from
`python/Galacticus/Build/Components/TreeNodes/__init__.py` (pointer members at
`_BASE_DATA_CONTENT`, currently `parent, firstChild, sibling, firstSatellite,
mergeTarget, firstMergee, siblingMergee, formationNode`); new members are
auto-nullified on creation via `CreateDestroy.py`. Add:

- `lastSatellite` (`type(treeNode), pointer`): tail of the satellite list.
- `nodeSatelliteEarliest` (`type(treeNode), pointer`) and
  `timeSatelliteEarliest` (`double precision`): cache of the
  earliest-in-time satellite and the time at which it was cached.
- (Phase 2) `nodeMergeeEarliest` / `timeMergeeEarliest`: same pattern for the
  mergee list, caching min over `timeOfMerging()`.

Memory cost: 3-4 pointers + 2 doubles ~ 40 bytes/node against a measured
~3 kB/node - ~1%.

### A. Tail pointer (exact, bitwise-neutral)

Maintain `lastSatellite` at every list mutation:

- Append sites (the four callers above): `host%lastSatellite%sibling => node`
  (or `firstSatellite => node` when empty), then `host%lastSatellite => node`.
- `Tree_Node_Remove_From_Host`: when the removed node is the tail, walk is
  already in progress there - set `lastSatellite` to the predecessor found by
  the existing loop (no extra cost).
- List transfer in `standardPromote` and destruction/merger sub-subhalo moves:
  splice using both ends.
- Sites that null or rebuild `firstSatellite` wholesale must null/rebuild the
  tail too: `merger_trees/construct/read/_class.F90:3083,3091`,
  `merger_trees/pruning/utilities.F90:93`,
  `nodes/operators/satellite/subsampling.F90:204`,
  `nodes/operators/analyses/node_formation_time/Cole2000.F90:213`.
- Reimplement `Tree_Node_Get_Last_Satellite` to return the stored pointer;
  keep the walk as a consistency assertion under a debug build flag so any
  missed maintenance site is caught loudly.

### B. Earliest-satellite cache (exact, bitwise-neutral)

Invariant exploited: a newly attached satellite enters at a time <= the
current minimum (it becomes a satellite at its host's current time), and
satellite times only increase.

- On attach (same four sites): if cache empty or `time(new) <=
  timeSatelliteEarliest`, set cache to the new node.
- On removal: if the removed node is `nodeSatelliteEarliest`, invalidate
  (null the pointer).
- On query (`standardTimeEvolveTo` l.963): if the cache is set and the cached
  node's current `basic%time()` equals `timeSatelliteEarliest`, the cached
  value is the minimum - O(1). Otherwise (the argmin advanced, or cache
  invalid) rescan once, storing the new argmin. Rescans therefore occur only
  when the binding constraint actually changed - amortized O(1) per host
  advance instead of O(S) per consideration.
- The scan fallback stays as the single implementation point
  (`timeSatelliteEarliestGet(node)`), so behavior is identical by
  construction; results must be bitwise-identical to master.

Note for compatibility with the batched-merger-synchronization proposal
(separate issue): implement the cache over an *effective* limit time
(currently just `basic%time()`), so a later change to how parked mergees cap
the host slots into one place.

### C. Promotion splicing

With A in place, `standardPromote`'s transfer becomes an O(1) splice plus the
existing O(S) re-parent loop (`satelliteNode%parent => hostNode`,
l.1300-1304). Keep the re-parent loop for now (pointer writes only; it is the
smaller share of the measured cost). A follow-up investigation - out of scope
here - could remove it via identity-preserving promotion (the child node takes
over the parent's tree links instead of moving component contents into the
parent), but that changes node index/uniqueID continuity semantics visible to
outputs and event hooks, so it needs its own design discussion.

### D. (Phase 2) Earliest-mergee cache

Same pattern as B over `firstMergee`/`siblingMergee`, keyed on
`timeOfMerging()`. Maintenance sites are more scattered - mergee transfer
loops in `satelliteMergerProcess`
(`evolve/timesteps/satellite.F90`), `destructionProcess`
(`evolve/timesteps/satellite_destruction.F90`), `standardPromote`, mergee
add/remove in `objects/nodes/_class.F90`, and `timeOfMergingSet` calls (e.g.
`nodes/operators/physics/satellite_merging/radius_trigger.F90:222`) - so gate
this phase on profiling after A-C (the mergee-scan getter was 2.5% at
M/m=10^6 and growing).

## Validation and benchmarking

- **Bitwise-identical outputs** to master for A+B+C on: `quickTest.xml`, the
  DMO subhalo reference configuration at m_res=1e8 (fixed seed), and one
  read-from-file (N-body) tree test - this is a pure data-structure change
  with no numerical-path change, so any difference is a bug.
- Debug-flag assertion runs of the testSuite to catch missed maintenance
  sites.
- Benchmark with `performanceAnalysis/scripts/run_one.sh` at
  (10^13 M_sun; m_res = 1e8, 3.16e7) and (10^12 M_sun; m_res = 1e6, 3.16e5):
  expect the ~10% overhead at M/m=10^6 to largely vanish and the slope at
  3.16e6 to drop toward 1.0-1.05. Re-run the perf profile and confirm
  `basicStandardTimeGet` / promote/move shares return to <1%.

## Effort

Moderate: one generator-file change, ~8 maintenance sites, evolver query
rewrite, debug assertions. The enumerated-site list above is intended to be
exhaustive; the debug assertion is the safety net for anything missed.
