# Continuous integration and pull requests

*Read this before diagnosing a CI failure, judging a PR's checks, or editing a workflow.* It covers how the GitHub
Actions pipeline behaves and the ways it misleads. These entries record traps found in practice; each names when it was
last verified, so check that a named file, flag or command still exists before relying on it.

## CI "unrecognized parameter" warnings: where they show, and how to find them all
*Added 2026-11-11; verified 2026-11-11*

`testModel.yml` writes these warnings only to the job **summary**, which no GitHub API returns. A job log
(`gh api --allow-escape-sequences repos/.../actions/jobs/<id>/logs`; without the flag `gh` refuses the output) holds the
tee'd `test.log`, but not the `galacticus.log` files of models a test script launches itself, so scanning every log of a
run misses those. Before PR #1773 (2026-11-10), `parameterValidate.py` also skipped files with a `<parameterGrid>` root; it
now validates them, including `parameterLevel`.

**Why:** a log scan of all 503 jobs found none of the 50 warnings visible in one job summary.

**How to apply:** combine a log scan with the validator, run over every `<parameters>` block (regressions, grids, and XML
embedded in test scripts), and confirm with a run on a current binary.

## A stacked PR's MERGED label does not mean it reached master
*Added 2026-09-17; verified 2026-11-05*

A PR stacked on another targets its *base branch*, not master. If the base merges upward first, merging the stacked PR
puts its commits on a branch that master has already taken, so GitHub shows MERGED while master has nothing. On 2026-09-17
this stranded four of the seven dust PRs (#1489, #1491, #1493, #1497): #1487 merged at 15:04:30 and #1489 merged into
#1487's branch at 15:04:32, two seconds too late. Merging the stacked PR *first* does not protect it either: master uses a
merge queue, which merges the head the base PR had when it was queued. On 2026-10-22 #1643 merged into #1641's branch at
15:35:49, #1641 reached master at 15:36:09 with its older head (`cd6236e794`), and #1643's three commits were stranded.
**Prefer unstacked PRs.**

A stacked PR also gets no PR-Checks (spell check, test markers, formatting, ...): `prChecks.yml` fires only for PRs whose
base is `master`, and its event types omit `edited`, so retargeting to master doesn't start it either. Only CI/CD runs. So
#1742 reached master, and #1745 sat green, without it (2026-11-05). After retargeting, close and reopen the PR to run it.

**Why:** the whole stack was reported as merged and the work looked done; master had none of it. Recovery needed a further
PR (#1515) from the one branch that still held every commit.

**How to apply:** after a stack merges, verify on master, not on labels: `git merge-base --is-ancestor <commit> origin/master`
for one commit per PR, and check that a file each PR added exists. To find a branch holding the remainder, test candidate
heads the same way, and tell the maintainer. A recovery PR must be merged with a merge commit when any commit in it carries a
parameter-migration hash.

## A Milky Way benchmark "Performance Alert" on a PR is usually noise
*Added 2026-10-30; verified 2026-10-30*

`github-action-benchmark` posts a review comment when a PR's milkyWay wall time exceeds the previous *single* recorded
result by 1.40. That series is bimodal and flips between a ~40 s and a ~61 s state from one master commit to the next, so
the alert fires on unrelated PRs. Measured over the 128 recorded runs: median 63.3 s, min 22.8, max 94.5, between-commit
stdev 17.2 s, while each commit's own `range` is at most 1.6 s (median 0.111) - a given binary is reproducible, different
binaries land in different states. **19 of 127 commit-to-commit steps (15%) would trip the 1.40 threshold**, 9 slower and
10 faster.

The cause was diagnosed during the investigation of an apparent slowdown blamed on `cf4a1d0cf0` (PR #1625): the cost is
dominated by chance bursts of tiny ODE steps limited by `disk:standard:massGas`, and which build gets them turns on tiny
floating-point perturbations - that commit added an *always-zero* meta-property and moved the benchmark as much. The
root-cause fix `146103d19c` (`fractionMassImplausibleNegligible`) reached master 2026-10-22 and reduced the severity but
not the bimodality.

**Why:** it reads as a 1.5x regression caused by the PR, and invites hunting for a cost that is not there; on PR #1708 the
flagged 60.5 s was at the 45th percentile of the history, below the median, while the 40.0 s baseline was at the 12th.

**How to apply:** before investigating, get the history and place both numbers in it -
`gh api repos/galacticusorg/galacticus/contents/dev/bench/milkyWayModel/data.js?ref=gh-pages --jq .download_url`, then
parse `window.BENCHMARK_DATA`. Also check `git diff --name-status <base>..HEAD -- source/`: a PR which only *adds* files
changes no executed code path, and can move this benchmark only by perturbing the trajectory.

## A lost CI runner looks exactly like a test failure
*Added 2026-09-18; verified 2026-09-18*

A GitHub hosted runner which dies mid-job marks the job `failure`, with the annotation "The hosted runner lost
communication with the server". There is no test output to read, and the log blob is unfetchable while the run is still in
progress, so it presents as an unexplained failure in one test.

**Three signatures distinguish it from a real failure**, and all are visible without the log:
- **Duration.** Compare against the same job on an earlier run (`gh api .../actions/jobs/<id>` gives
  `startedAt`/`completedAt`). A numerical failure takes as long as a pass. On 2026-09-18 a job which normally ran 7 minutes
  took 47 and then "failed".
- **Step states.** `gh api repos/galacticusorg/galacticus/actions/jobs/<id> -q '.steps[]|...'`. If the job concluded
  `failure` but its steps are still `running`, it was killed mid-step. A test which really failed leaves `Run test` as
  `success`, because these tests exit 0 and report in their output.
- **Local reproduction.** The model simply runs and passes.

**Why:** a cross-PR comparison (passes on other PRs, fails on this one) and a before/after comparison on the same branch
both look like strong evidence that the change caused it. Neither can separate "our change" from "a bad runner that day";
on 2026-09-18 that reasoning wrongly blamed the change.

**How to apply:** check duration and step states *before* attributing a single-job CI failure to the change under test,
and prefer a re-run to a day of investigation. Re-running needs the run to be complete first - `gh run rerun <id> --failed`
refuses while it is still in progress.

## Master's CI/CD runs are almost always cancelled, by design
*Added 2026-09-17; verified 2026-11-23*

`cicd.yml` sets `concurrency: {group: cicd-${{ github.ref }}, cancel-in-progress: true}`, on master as well as on PR
branches. On master this means each merge cancels the run validating the previous one, and merges arrive faster than the
suite completes (five landed inside one minute on 2026-09-17). Of 30 consecutive master CI/CD runs checked then, **27 were
cancelled and 2 succeeded**. This is deliberate, maintainer-confirmed behavior (2026-11-23): changes are validated by their PR
checks, which run to completion.

**Why:** it looks as though master is verified after each merge, or as though cancelled runs signal a problem. Neither is
so, and waiting for a post-merge run can burn hours for nothing: a 12-hour watch on 2026-09-17 never saw one finish.

**How to apply:** judge a change by its *PR* checks; do not wait on a master run. If asked to confirm master is green after
a merge, explain this rather than watching. Do not propose changing the concurrency setting.

## galacticusorg Actions jobs can queue for hours
*Added 2026-09-15; verified 2026-09-15*

GitHub Actions runners are shared across the galacticusorg organization. When many Galacticus PRs have CI running, jobs in
other repos (e.g. the dendros release's "Publish to PyPI") can sit in `queued` for over an hour, even after environment
approval and with githubstatus.com green.

**Why:** on 2026-09-15 a queued dendros publish job looked stuck, and a cancel and re-run was nearly proposed; it was just
waiting its turn.

**How to apply:** don't treat a long `queued` state as stuck, and don't propose cancelling or re-running it. Poll with a
long, lightweight background watcher (curl every 60 s, hours of timeout).

## GitHub Actions caches: exclusions, scope and eviction
*Added 2026-10-13; verified 2026-10-13*

- **`!` exclusions in `actions/cache` paths only filter what the include patterns resolve to.** It resolves with
  `implicitDescendants: false` (`cacheUtils.resolvePaths` in actions/toolkit) and archives each resolved path whole, so
  `path: dir` plus `!dir/big*` caches everything. Use `dir/*` plus `!dir/big*`. A 3.2 GB "slim" cache proved it.
- **Scope:** an entry is readable only from the branch (or PR merge ref) that saved it, and from the default branch's.
  Galacticus master's CI is nearly always cancelled (see above), so there are almost no master-scoped caches to share;
  every `testModel.yml` cache is PR-scoped.
- **Eviction:** the galacticusorg/galacticus 10 GB budget is usually full of 2-3 GB PR caches, and least-recently-used
  eviction removed new entries within hours (2026-10-12/13).

**Why:** a cache designed to carry tables between weekly PGO runs saved 3 GB entries that were evicted before any later run
could use them.

**How to apply:** list a directory's entries when excluding; check entry sizes with `gh cache list` and usage with
`gh api repos/<repo>/actions/cache/usage`; don't count on a cache surviving across a weekly schedule. Use an artifact if
data must persist.
