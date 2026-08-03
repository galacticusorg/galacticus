"""Run an independent per-source-file scan concurrently across a process pool.

Several build scripts walk the `source/` tree and parse every file. The per-file
cost is dominated by blocking file I/O (open/read over NFS); on a loaded node
each open stalls in close-to-open metadata revalidation, so doing them one at a
time serialises thousands of millisecond waits. Scanning the files concurrently
overlaps those waits (and the parsing CPU work).

A *fork*-based process pool is used so workers are fully isolated -- callers need
make no assumptions about the thread-safety of the parser -- and the read-only
inputs a worker needs are inherited via copy-on-write rather than pickled per
task. Results are returned in the SAME order as the input tasks, so callers that
depend on deterministic ordering (e.g. Makefile emission) are unaffected.

Andrew Benson (2026).
"""

import concurrent.futures
import multiprocessing
import os
import re
import sys

from Galacticus.Build.Jobserver import slots as jobserver_slots

# Worker count when run standalone (no `make`, no explicit override). Capped at
# the task count by `resolve_jobs`. Oversubscribing cores is intentional: the
# work is I/O-latency bound, so workers blocked on NFS still overlap waits.
_DEFAULT_STANDALONE = 16

_MAKEFLAGS_J_RE = re.compile(r'(?:^|\s)-j(\d*)')


def _cpu_count():
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def resolve_jobs(n_tasks):
    """Decide how many workers to use.

    Precedence:
      * `GALACTICUS_BUILD_JOBS` (explicit override; 0/"unlimited" -> cpu count);
      * the `-j` value `make` passed down via `MAKEFLAGS` (so the build's own
        `-jN` flows through -- bare `-j` means unlimited -> cpu count, and a
        `make` invocation with no `-j` is serial -> 1);
      * when not run under `make` at all, `_DEFAULT_STANDALONE`.
    The result is clamped to `[1, n_tasks]`.
    """
    override = os.environ.get('GALACTICUS_BUILD_JOBS')
    if override is not None and override.strip() != '':
        try:
            jobs = int(override)
        except ValueError:
            jobs = 0
    else:
        makeflags = os.environ.get('MAKEFLAGS')
        if makeflags is None:
            jobs = _DEFAULT_STANDALONE          # standalone (no make)
        else:
            m = _MAKEFLAGS_J_RE.search(makeflags)
            if m is None:
                jobs = 1                         # under make, no -j -> serial
            elif m.group(1) == '':
                jobs = 0                         # bare -j -> unlimited
            else:
                jobs = int(m.group(1))           # make -jN -> N

    if jobs <= 0:
        jobs = _cpu_count()
    return max(1, min(jobs, n_tasks))


def scan(tasks, worker, script_name, jobs=None):
    """Apply `worker(task)` to every task, returning results in task order.

    `worker` must be a module-level (picklable) callable. `script_name` is used
    only to make failures legible. Falls back to a serial loop when one worker
    is requested or there is a single task. On a worker failure the underlying
    error is surfaced with a hint to re-run serially for a full traceback,
    rather than a raw pool traceback.

    When run under a `make` jobserver (from a `+`-marked recipe) the requested
    worker count is further capped by the job slots make actually has free, so a
    scan running alongside other recipes cannot oversubscribe the build.
    """
    tasks = list(tasks)
    if not tasks:
        return []
    if jobs is None:
        jobs = resolve_jobs(len(tasks))
    if jobs <= 1:
        return [worker(task) for task in tasks]

    with jobserver_slots(jobs) as granted:
        if granted <= 1:
            return [worker(task) for task in tasks]
        return _scan_pool(tasks, worker, script_name, granted)


def _scan_pool(tasks, worker, script_name, jobs):
    """Run `tasks` across a fork-based pool of `jobs` workers, in task order."""
    results = [None] * len(tasks)
    ctx = multiprocessing.get_context('fork')
    try:
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=jobs, mp_context=ctx) as pool:
            futures = {pool.submit(worker, task): i
                       for i, task in enumerate(tasks)}
            for future in concurrent.futures.as_completed(futures):
                results[futures[future]] = future.result()
    except concurrent.futures.process.BrokenProcessPool as exc:
        sys.exit(f"{script_name}: a parallel worker was killed before it could "
                 f"return (out of memory, or a crash in the parser). Re-run "
                 f"with GALACTICUS_BUILD_JOBS=1 to reproduce serially.")
    except BaseException as exc:  # noqa: B036  (re-raised cleanly below)
        sys.exit(f"{script_name}: a parallel worker failed: {exc}\n"
                 f"Re-run with GALACTICUS_BUILD_JOBS=1 for a full serial "
                 f"traceback.")
    return results
