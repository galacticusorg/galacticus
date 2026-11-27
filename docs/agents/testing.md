# Testing and verifying

*Read this before comparing model outputs, writing or trusting a test, or benchmarking.* It collects traps that make a
correct change look broken, or a broken one look correct. These entries record traps found in practice; each names when
it was last verified, so check that a named file, flag or command still exists before relying on it.

## A Python test failing only in CI: suspect OpenBLAS kernels and thread count
*Added 2026-10-31; verified 2026-10-31*

OpenBLAS picks its kernels per CPU, and with the Haswell and Zen kernels its rounding depends on the BLAS thread count;
with SkylakeX (the default on a Zen4 machine) it did not. `test_train_parallel_matches_serial` (emulation) failed only in
CI for this reason, and passed locally at every thread count until reproduced with
`OPENBLAS_CORETYPE=Haswell OPENBLAS_NUM_THREADS=4 python -m pytest ...`.

**Why:** a test demanding bit-identical results across processes with different BLAS thread counts passes locally and
fails on the CI runners (about 4 cores).

**How to apply:** reproduce with `OPENBLAS_CORETYPE` set to Haswell or Zen and a thread count of 2–64. Fix such a test by
running both sides under `threadpoolctl.threadpool_limits(1)` rather than by loosening its tolerance.

## Nondeterminism: what differs between runs, and how to compare
*Added 2026-06-23; verified 2026-10-15*

### Comparing outputs

- **lightcone model:** runs of the same binary with 8 threads differ in ~85 of 123 datasets, because RNG draw order (halo
  spins) depends on thread scheduling. Geometry and N-body quantities are stable. Compare `OMP_NUM_THREADS=1` runs after
  sorting by `nodeIndex`; they should be bit-identical. Row order also varies, and SFH *mass* columns differ even
  single-threaded, so exclude them.
- **validate_darkMatterOnlySubHalos:** bins differ by up to ~88% between runs, and the output has only analyses. Compare
  each analysis's `logLikelihood` against the spread between two runs of the same binary.
- **Model outputs are not row-comparable across runs.** With OpenMP, galaxies are written in thread-dependent order, so two
  runs of the same executable differ row by row. For before/after "no change" checks run both with `OMP_NUM_THREADS=1`
  (bit-identical on repeat), or match rows on a unique key such as `basicMass`.
- **Separate `-O3 -flto` builds:** these differ at ~1e-9 from FP codegen. For bit-exact comparison, rebuild both with
  `-ffp-contract=off`. **The seed is that small; the output difference need not be.** Measured 2026-10-13 on `quickTest`:
  relocating one procedure within `source/utility/IO/HDF5/_module.F90`, body byte-for-byte identical, moved 25 of 36
  `nodeData` datasets, by up to **1160%** in `spheroidAngularMomentum`, from single-ulp seeds (`darkMatterProfileScale`
  3.3e-16) amplified through the baryonic evolution. Doing the same in `merger_trees/outputter/standard.F90` moved 22 of 36
  by up to 2.9%. Whether a given edit shifts codegen at all is arbitrary - a *larger* edit to that outputter changed
  nothing, a smaller one did - so a clean default-build comparison is luck, not evidence. Integer datasets and tree
  structure (`nodeIndex`, `mergerTreeIndex`, `mergerTreeStartIndex`, `mergerTreeCount`) stayed identical throughout.
- **Round-off changes in columns you didn't touch** have two known causes:
  - The 20-slot `node%massDistribution` memo cache: changing one extractor's requests changes cache eviction, and iterative
    profiles reconverge slightly differently. Test whether the diff disappears when the intended change is removed.
  - Adding `target` to a hot array such as `nodeRates`, which changes aliasing for all generated code. Put
    `target, contiguous` on the *dummy argument* only.
- **Parameter-file copies:** `xi:include` paths are relative to the parameter file, so copies placed elsewhere break.
- **Output file names:** never run two models writing the same `outputFileName` at the same time.

### Dynamic-data cache traps

- **First run after a cache miss:** a run that retabulates a quantity (e.g. `-> retabulating σ(M)` /
  `writing σ(M) data to: …/cache/…` in the log) uses the fresh in-memory table, while every later run reads it back from
  `GALACTICUS_DYNAMIC_DATA_PATH` and differs at ~1e-7, even in `basicMass` and `mergerTreeWeight`. Seen on 2026-09-18
  (`molecularCooling.xml`): 154 of 272 datasets "differed" until each binary was run twice. Discard a reference run whose
  log shows retabulation, or run each binary twice.
- **Cache keys include source digests, so builds do not share tables.** Tabulations whose file names use
  `hashedDescriptor(includeSourceDigest=.true.[, includeFileModificationTimes=.true.])` (spherical collapse, σ(M), adaptive
  SFH, …) get new keys whenever a rebuild changes the digests - e.g. any edit reaching `Galacticus_Nodes`. So do the SED
  templates of `nodePropertyExtractorSED`
  (`$GALACTICUS_DYNAMIC_DATA_PATH/stellarPopulations/nodePropertyExtractorSED_<md5>_1.hdf5`), which append
  `:sourceDigest{…}` for `SED.F90` itself, so **any** edit to that file renames every template and the cached files of two
  builds can never be compared directly - compare the model output instead (verified 2026-10-23).
  The new binary then tabulates afresh while the old one reads tables extended by earlier runs, so outputs "differ" (even
  `basicMass`), and running each binary twice does *not* fix it: a binary's own runs can differ too, as the shared cache is
  rewritten between runs. Confirmed 2026-09-21: two binaries 280/341 and 37/102 "different" became bitwise identical once
  each run got a **private copy of a frozen cache holding only files whose names both binaries open** (list them with
  `strace -f -e trace=openat`, take the intersection, `cp -a` a fresh copy per run, set `GALACTICUS_DYNAMIC_DATA_PATH` to
  it). Probably also the cause of past "incremental build gave different output" episodes.
  It can also flip the **pass/fail of a borderline test**: on 2026-09-22 the fuzzing suite's "failed evolution" case
  (`rateMaximumExpulsion=1e30`) ran to completion on master with the shared cache but failed with an empty one. For
  `quickTest`, an **empty** private cache per run is enough - master and branch were then bit-identical.
  Confirmed again 2026-10-12: editing a function the model *never executes* (a scaler moment method) renamed the
  σ(M)/spherical-collapse tables, so even `-ffp-contract=off` builds diverged (node counts differed) until each got a frozen
  private cache. Filter the strace list to regular files (`[ -f ]`): it also lists directories, and `cp -a` of the cache
  root copied all 29 GB. An **empty** cache is not always safe either: there, a run hung forever at a RecFast "Enter the
  choice of modification for HeI" prompt (master too); always run comparisons with `< /dev/null` and a `timeout`.
- **The converse trap (added 2026-10-15):** a type's composite digest covers only its own source, its parent types and
  same-family members - *not* `class(xClass)` members of other families (`Utils.py:class_dependencies`). So editing, e.g.,
  `transferFunctionFile` does **not** rename the σ(M) cache, and a stale table (with its stored `sigmaNormalization`) is
  reused. A change to what a member object returns needs explicit invalidation (a `fileFormat` attribute was added to the
  σ(M) file for this reason).

**Why:** these differences look like regressions when they aren't.

**How to apply:** pick the comparison method before diffing outputs.

## Reading Galacticus HDF5 outputs in test scripts
*Added 2026-09-21; verified 2026-09-21*

Moved to the user guide: *Reading outputs from a script* in `docs/manuals/user-guide/running.rst`. Read it
before writing a test which reads model outputs.

## Trees evolve from the highest tree number down
*Added 2026-08-10; verified 2026-09-13*

"Evolving tree number 128" is the *start* of a 128-tree run.

**Why:** reading the numbers as ascending inverts how far a run has progressed.

**How to apply:** judge progress, or whether a run passed a crash point, by the *lowest* tree number reached.

## Runs with nearby seeds share their trees
*Added 2026-11-17; verified 2026-11-17*

Moved to the user guide: *Random seeds and independent runs* in `docs/manuals/user-guide/running.rst`. Space base
seeds by more than the tree count, or use `<mergerTreeSeeds value="random"/>`.

## Benchmarking and training runs: what the models really do
*Added 2026-09-29; verified 2026-10-11*

- **`validate_milkyWay.xml` evolves only ~1 tree in 4.** `filterLMC.xml`'s `anyNode` `mergerTreeFilter` skips every built
  tree without an LMC analog *before* evolution. With few trees a seed can evolve **none** (2 trees, seed 7349: 16 s, no
  "Evolving tree" lines). Count "Evolving tree" lines before trusting a run as a benchmark or a PGO training input.
- **`-fprofile-generate` binaries are ~9x slower** on the Milky Way model (3 trees in ~2 h single-threaded), ~7x on the DMO
  reference.
- **`-fstack-arrays` segfaults on large forests** (the UNIT lightcone trees, in `mergerTreeConstructorRead` building node
  arrays) at the default 8 MB stack; it runs under `ulimit -s unlimited`.
- **PGO or `-march` change outputs** (last-bit changes amplified by the ODE solver), so their timings partly measure a
  different trajectory; `-ffp-contract=off` restores bitwise-identical DMO under `-march`. A PGO build can change even DMO
  output, via the random tree builder (2026-10-09: 0.6% of halos differ) - compare cycles per instruction.
- **An empty `GALACTICUS_DYNAMIC_DATA_PATH` is not a cold start.** Tables are also found under
  `$GALACTICUS_DATA_PATH/dynamic` (e.g. `stellarPopulations/simpleStellarPopulationsFSPS*`) and FSPS keeps its outputs in
  the tools path (`fsps-*/OUTPUTS`), so locally a "cold" run started evolving in 4 min, while on CI's fresh checkout the
  galaxy formation model spent hours building FSPS SSPs serially first. To test a cold start, point `GALACTICUS_DATA_PATH`
  and `GALACTICUS_TOOLS_PATH` at fresh locations too.

**Why:** each of these silently invalidated a benchmark or training run in the #1449 flag study.

**How to apply:** check tree counts and exit codes of every training/benchmark run, and pair each variant with a baseline
run in the same round.

## testSuite scripts using `launch.py` skip models whose output directory exists
*Added 2026-10-30; verified 2026-10-30*

`scripts/aux/launch.py` (local method) skips any job whose `testSuite/outputs/<name>/` directory already exists, silently.
So rerunning such a test (e.g. `test-metal-conservation.py`) after swapping in another `Galacticus.exe` just re-checks the
previous binary's outputs.

**Why:** on 2026-10-30 a "without the fix" mutation run passed, and looked like proof the test couldn't detect the bug.
With the outputs moved aside, it failed in five of its fixed-metallicity cases.

**How to apply:** before every run of a launch.py-based test, move or delete its `testSuite/outputs/<name>*` directories,
and check the log shows the models actually ran.

## Failures that aren't your change
*Added 2026-08-01; verified 2026-09-15*

- **Poisoned Blitz2006 cache.** `Fatal error: attribute 'coefficientFactorBoostMinimum' does not exist` (or a similar
  missing HDF5 attribute) means an interrupted run left the star-formation cache without attributes. Check with `h5ls -v` on
  `<cache>/starFormation/starFormationRateSurfaceDensityDisksBlitz2006_*.hdf5`, then move the file aside and rerun.
- **Corrupted CAMB Makefile.** `failed to build CAMB code` (`Error 127`) means a missing forutils stamp re-ran a sed that
  isn't idempotent, leaving `gfortran-16-16` in `<tools>/CAMB-*/fortran/Makefile`. Fix the compiler name there rather than
  rerunning.
- **Network downloads.** Download failures (e.g. recfast from astro.ubc.ca, Zenodo) are usually transient.
- **CI artifact download.** A `Test-Model` job failing with
  `Unable to download artifact(s): Failed to ListArtifacts: Unable to make request: ENOTFOUND` never ran the model. It is
  GitHub infrastructure, not a test result, and needs only a re-run (seen on PR #1496, 2026-09-15).
- **GitHub 5xx in a post-test step.** A `Test-Code` job whose log ends `Tests passed: N (100%)` and then
  `HttpError: <!DOCTYPE html>` (a GitHub "Unicorn!" error page, from the `qoomon/actions--context` step) passed its tests;
  only the reporting step failed. Re-run it (seen on PR #1497, 2026-09-16).

See also [`ci.md`](ci.md) for a CI runner that is lost mid-job.

**Why:** these look like regressions caused by whatever you just edited.

**How to apply:** rule these out before debugging your own change.

## First use downloads and builds tools; prebuild them where nodes lack internet
*Added 2026-10-18; verified 2026-10-18*

With an empty `GALACTICUS_TOOLS_PATH`, a model with a CAMB transfer function (e.g. the MDPL and Symphony PMF models)
downloads CAMB 2.0.3 and forutils from GitHub and compiles them, then spends ~15 min per cosmology on the transfer
function and σ(M) tables. A model whose transfer function is read from a file (COZMIC) fetches nothing. A
`<task value="buildToolCAMB"/>` parameter file (all tools: `parameters/buildTools.xml`) downloads and builds CAMB alone in
a few minutes.

**Why:** cluster compute nodes may have no internet access, so the first jobs would fail with `unable to download CAMB`.
And many jobs starting cold all wait on the same table locks.

**How to apply:** on a new machine, run the `buildTool*` task on a login node first, then warm the tables with a small run
per cosmology before fanning out jobs.

## CI Python environments are minimal; yours is not
*Added 2026-09-16; verified 2026-10-22*

A test or notebook that imports a package indirectly can pass locally and fail in CI. Two distinct traps, both hit on PR
#1497 (2026-09-16), one after the other:

- **Missing.** CI installs only what its workflow step names, and optional *extras* are not dependencies: dendros renders
  `list_outputs(format='tabulate')` through pandas, an extra, and the local environment happened to have pandas.
  `ModuleNotFoundError: No module named 'pandas'`.
- **Present but broken.** The CI container mixes Pythons: pip-installed packages live in
  `/usr/local/lib/python3.14/dist-packages`, the distribution's in `/usr/lib/python3/dist-packages`. Importing the
  distribution's astropy there raises `AttributeError`, so the astropy fallback of the same call failed too.

Two more, in the Test-Preprocessor job:

- **Test-Preprocessor's pytest** (whole-repository `python -m pytest`) has only the packages its "Install pytest" step in
  `cicd.yml` lists (no numpy until PR #1554 added numpy and h5py). A test module importing anything else fails
  *collection* (pytest exit status 2) and takes the whole job down (2026-09-23). Check that list, and reproduce with a venv
  holding exactly it.
- **No editable install in CI.** There `python/` is importable only because `python/conftest.py` prepends it to
  `sys.path`, inside the pytest process. A test that runs a build script as a **subprocess** must put `<repo>/python` on
  the child's `PYTHONPATH` itself, as the Makefile does, or it fails with `No module named 'Galacticus'`. Locally,
  `pip install -e .` hides this (PR #1666, 2026-10-22). Reproduce by running the *whole* suite (`python -m pytest` from the
  repo root) in the venv. Running one file skips `python/conftest.py` and fails even at collection.

**Why:** each cost a 1.5-2 hour CI job, and the second was masked because the test truncated the traceback before its
final line.

**How to apply:** prefer an API path needing no optional renderer at all — print what the notebook needs directly. To
check, re-run with stub modules on `PYTHONPATH` (`pandas.py`, `astropy.py` … each containing `raise ImportError(...)`) and
confirm the old call fails under the same block. When printing a captured traceback, print its end, not just its start.

## An independent reference must match Galacticus' unit definitions and the model's own inputs
*Added 2026-09-17; verified 2026-11-02*

Moved to the developer guide: *Writing an independent reference for a test* in
`docs/manuals/developer-guide/continuous-integration.rst`. Read it before writing a reference calculation that a test
compares against.

## Test coverage traps
*Added 2026-06-18; verified 2026-09-17*

- **`test-python-utils.py`** is outside `pytest python/`. It deep-imports private `Build/SourceTree/Process` symbols, so a
  rename can pass pytest and break CI. Grep it for any symbol you remove. It also defines shared helpers between test
  functions, so don't remove tests by slicing between `def` lines.
- **Whole-repository `pytest`.** CI's Test-Preprocessor job runs `python -m pytest` from the repo root with `BUILDPATH`
  unset, collecting `scripts/*/tests/` too, not just `python/`. A build-system change can break a synthetic fixture far
  from the code touched: enforcing directive schemas broke a `<functionClass>` fixture in `scripts/doc/tests/`, and the
  docs extractor reported it only as a warning. Run that exact command before pushing a build-system change.
- **Unit-test executables.** Each `source/tests/*.F90` is its own executable (`tests.<name>.exe`), not part of
  `Galacticus.exe`. Build the ones touched by an interface change. `libgalacticus` is a third target with its own generated
  Fortran. Unit-test programs must call `Error_Handler_Register()`, or GSL errors abort regardless of `status=`. Those that
  read HDF5 (directly or through a class constructor) must call `ioHDF5AccessInitialize()` first, or they segfault in
  `mutexset` at the first `hdf5Access%set()`.
- **Spell-check.** A clean local docs spell-check doesn't mean CI passes:
  - contributor names need a full git history;
  - words naming importable modules are silently ignored locally;
  - `sphinx -b spelling` exits 0 even when it flags words (check for non-empty `*.spelling` files);
  - stale generated pages cause local-only warnings;
  - a parameter name written `[likeThis]` (the `Error_Report` convention) in an `!![ ]!!`/`!!{ }!!` description passes the
    commit hook but fails CI (PR #1693, 2026-10-27). In RST docs write ``` ``likeThis`` ```.

  Never whitelist a British spelling in `aux/words.dict`.

  Where `sphinxcontrib-spelling` and `pyenchant` are installed, *Spell-Check-RST* can be reproduced exactly:
  `SPHINX_SPELLING=1 python -m sphinx -b spelling docs docs/_spelling`, judged by a non-empty `*.spelling` file rather than
  the exit status. Do that before and after a fix, so a clean run means something. `docs/conf.py` enables the extension
  only under `SPHINX_SPELLING=1` and reads `aux/words.dict`.
- **Unicode homoglyphs defeat the dictionary.** `aux/words.dict` whitelists `μm` with U+03BC GREEK SMALL LETTER MU, and the
  tree uses that form (15 files against 2), but U+00B5 MICRO SIGN renders identically and is a *different* word to enchant,
  so it is flagged. Fix by normalizing the character to U+03BC, never by adding the twin to the dictionary - that would
  leave both acceptable and let the failure recur. A `grep` for one form silently misses the other, so search by
  codepoint. Cost a red *Spell-Check-RST* on PR #1518 (2026-09-17), inherited from the dust tutorial merged in from master.
- **Notebook tutorials: edit the source, not the path CI names.** The tracked notebooks are `tutorials/library/*.ipynb` and
  `tutorials/models/*.ipynb`; `docs/conf.py` copies them to `docs/tutorials/<group>/`, which is what warnings cite and what
  `.gitignore`'s `docs/tutorials/*/*.ipynb` covers.
- **A test existing is not coverage: read what it asserts.** Several unit tests assert *internal consistency* rather than
  values, and look like coverage until read. `tests.initial_mass_functions` integrated each IMF between its own mass limits
  and required unit mass - satisfied by any shape the constructor normalizes, at only the default parameters; a
  normalization defect in the Chabrier IMF survived beneath it for years and was found the moment the test compared against
  an independently derived shape (2026-09-15). `tests.spin_distribution_nbody_errors` asserts only that its tabulation is
  reproducible. `tests.stellar_populations` compares against two numbers its own comment calls "a previously computed
  value". Before concluding a subsystem is covered, check whether its test compares against anything external.
- **A node memoizes the mass distributions built from its components.** In a unit test, resetting a component property
  (`disk%radiusSet`, …) and then calling the consuming class' own `calculationReset` leaves the *old* mass distribution in
  place, silently and with no error - in a model it is the `calculationResetEvent` that flushes it. Build one node per
  geometry instead. Cost a "fully molecular" limit that was evaluated on the uncompressed disk and wrong by 1,570x
  (2026-09-15).
- **A test that aborts loses its output, so a crash carries no failure marker.** Unit tests are built with
  `-ffpe-trap=invalid,zero,overflow`, so an expression such as `f/(1-f)` used to recover a quantity raises SIGFPE the moment
  the code under test returns the degenerate value; a NaN reaching GSL trips its error handler instead. Either way the
  executable dies *before* any buffered `MM:` output is flushed, so the log shows a crash with no marker and CI's grep sees
  no failure. Guard any such recovery and assert a value that cannot match instead. Seen twice on 2026-09-15/16.

**Why:** each of these let a real failure through, or produced a false one.

**How to apply:** check this list before trusting a green test result.

## Exercising an output analysis in a test model
*Added 2026-10-06; verified 2026-10-06*

- `<outputAnalysis>` is only constructed under the `analyzer` merger tree outputter. With `quickTest.xml`'s `standard`
  outputter it is silently ignored: no crash, no `analyses` group, even under `--dry-run`. Wrap both in
  `<mergerTreeOutputter value="multi">` and add an `outputTimes` list.
- Array parameters take one space-separated `value="1.5 0.5"`. Repeating the element keeps only the first, and the sole
  sign is `WARNING: problems found with input parameters` in the log.
- `quickTest.xml` forms no galaxies above ~10¹⁰ M☉. To populate high-mass bins, use the analysis's own
  `systematicErrorPolynomialCoefficient` (offset + slope) rather than raising tree masses.

**Why:** on 2026-10-06 (issue #1565) the first two cost a false "can't reproduce" and a false "the polynomial slope is
ignored" lead.

**How to apply:** check for the `analyses` group in the output before concluding an analysis ran, and grep the log for that
warning after adding parameters.

## `perf report` and `valgrind` hang at 0% CPU when `DEBUGINFOD_URLS` is set
*Added 2026-07; verified 2026-10-20*

Where `DEBUGINFOD_URLS` is set system-wide (common on recent Linux distributions), `perf report`/`script`/`annotate`
block forever trying to fetch debug symbols, and so does `valgrind`, which sits waiting on a child `debuginfod-find`
(`wchan=poll_schedule_timeout`; seen 2026-10-16 on a memory-leak run). An agent's sandbox often cannot reach the
server, so the fetch never completes.

**How to apply:** clear it for the command: `DEBUGINFOD_URLS= perf report ...`, or `export DEBUGINFOD_URLS=` before
`valgrind` (including `testSuite/test-memory-leaks.py`).

## Profiling without `perf`: sample stacks with `eu-stack`
*Added 2026-10-23; verified 2026-10-23*

Where `perf` is not installed (common on clusters), a poor-man's profiler needs no root as long as
`/proc/sys/kernel/yama/ptrace_scope` is `0`. Run the model in the background, then loop
`eu-stack -p $pid >> stacks.txt` with `sleep 1` between samples (`eu-stack`, `gstack` and `pstack` are in elfutils).
`eu-stack` fully symbolizes a `-flto` Galacticus binary, including the OpenMP outlined frames (`..._omp_fn.0`) and GSL,
at well under a second per sample. Aggregate by counting, over samples, the innermost frame (self cost) and the samples
whose stack contains a symbol anywhere (inclusive cost). About 180 samples was ample to separate an 80%-inclusive symbol
from a 17% one.

**Why:** without `perf` the run looks unprofilable, which invites a `-pg` rebuild (slow, and it perturbs the hot path)
or guessing.

**How to apply:** use this to find where a model's time goes before optimizing it. Prefer `OMP_NUM_THREADS=1`, so each
sample is one stack.

## On a shared node, an oversubscribed OpenMP run looks deadlocked
*Added 2026-11-10; verified 2026-11-10*

A run with more threads than the node's free CPUs spins at its OpenMP barriers: a 2 s solve at 96-128 threads ran for
over 9 minutes (8.5 CPU-hours), and `eu-stack` showed threads waiting at different, consecutive barriers. With
`OMP_WAIT_POLICY=passive` it finished in 27 s.

**Why:** threads at different barriers look like a worksharing mismatch in the code, which invites a hunt for a bug.

**How to apply:** check the load (`top`, or `ps -eo pcpu,user,comm --sort=-pcpu | head`) before a thread-scaling test,
and stay below the free CPUs or set `OMP_WAIT_POLICY=passive`. A real deadlock still hangs with passive waiting.
