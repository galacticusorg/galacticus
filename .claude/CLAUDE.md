# CLAUDE.md

Guidance for [Claude Code](https://docs.claude.com/en/docs/claude-code) (and
other AI assistants) working in this repository. This is a short operational
quick-reference; it does **not** restate the full documentation. For anything
beyond the essentials below, defer to:

- [`README.md`](../README.md) — quickstart, prerequisites, troubleshooting.
- [`CONTRIBUTING.md`](../CONTRIBUTING.md) — contribution workflow, commit format,
  attribution, and the AI-assisted-contribution policy (**all AI-generated work
  must be verified by a human before submission**).
- The [developer guide](https://galacticus.readthedocs.io/en/latest/manuals/developer-guide/index.html)
  (`docs/manuals/developer-guide/`) — the authoritative reference for coding
  style ([`coding.rst`](../docs/manuals/developer-guide/coding.rst)), the
  preprocessor-directive / `functionClass` code-generation system, and the
  workflow for
  [adding a new class](../docs/manuals/developer-guide/creating-a-new-class.rst).

Galacticus is a semi-analytic model of galaxy formation, written primarily in
Fortran with a Python-based build/code-generation system.

## Environment setup

Run these from the repository root before building or running:

```bash
export GALACTICUS_EXEC_PATH=`pwd`   # the repo root; required — see below
pip install -e .                     # puts python/ on the import path + installs deps
```

- **`GALACTICUS_EXEC_PATH`** must point at the repository root. It is how the
  executable and the build/analysis scripts locate the source tree and its data.
  If it is unset the code silently falls back to `./` (the current directory),
  which usually causes hard-to-diagnose failures when you are not sitting in the
  repo root — so always export it. (Resolved in `source/input/paths.F90`.)
- **`pip install -e .`** makes the modules under `python/` importable and pulls
  in numpy/scipy/h5py/lxml/matplotlib/… Alternatively,
  `export PYTHONPATH=$(pwd)/python:$PYTHONPATH`. Use `pip install -e '.[test]'`
  to also get `pytest`.
- **`GALACTICUS_DATA_PATH`** points at the separate
  [`galacticusorg/datasets`](https://github.com/galacticusorg/datasets) repo.
  It is not needed to *compile*, but is required to *run* real models and most
  of the model regression tests.

## Building

From the repository root:

```bash
make -jN Galacticus.exe   # N = number of parallel cores
```

Gotchas:

- **Builds are slow and memory-hungry.** Concurrent `-O3` compiles of the large
  generated `*_class` units are the peak. On a ~16 GB machine, `-j4` can OOM
  (exit 137); CI caps at **`-j2`**. On memory-constrained machines use a lower
  `-jN` and/or `LTO=disabled` (e.g. `make -j2 LTO=disabled Galacticus.exe`).
- **One `make` per tree.** Never run two `make` invocations concurrently in the
  same working tree — they race over generated files.
- Because a full build takes a long time, **prefer launching it in the
  background** and continuing other work, rather than blocking on it.

## Waiting on long builds & models (avoid the `pgrep` trap)

Builds and models can run for a long time, so it is tempting to poll with
`pgrep`. A naive `pgrep -f <pattern>` frequently **matches the poller's own
command line** (the shell/command running the `pgrep`), so it looks like the job
is "still running" forever, or matches unrelated processes. Prefer, in order:

1. **Let the harness manage it.** Claude Code can run a command in the
   background and notify you when it exits — use that instead of busy-polling.
2. **Capture the PID and wait on it:**
   ```bash
   make -j2 Galacticus.exe > build.log 2>&1 &
   build_pid=$!
   wait "$build_pid"; echo "build exited $?"
   ```
3. **Use a sentinel** — have the job `touch done.flag` (or write a final log
   line) on completion, and check for that.

If you must use `pgrep`, make it exact and exclude yourself:

```bash
pgrep -x Galacticus.exe                 # exact process-name match, not full cmdline
pgrep -f Galacticus.exe | grep -vw "$$" # or drop the current shell's PID
```

Never rely on a bare `pgrep -f <pattern>` whose pattern can appear in the
polling command itself.

## Running a model & the smoke test

```bash
./Galacticus.exe parameters/quickTest.xml
```

`quickTest.xml` is a deliberately small model that completes quickly — the
standard smoke test after a build. A successful run exits `0` and writes an
HDF5 file; since `quickTest.xml` sets no output name, it defaults to
`galacticus.hdf5` in the working directory. Set an `outputFileName` parameter to
change that.

## Testing

- **Python unit tests** (build-system / tooling modules — no model runs):
  ```bash
  export GALACTICUS_EXEC_PATH=`pwd`
  python -m pytest -v
  ```
- **Model regression tests** (require a built `Galacticus.exe` and
  `GALACTICUS_DATA_PATH`). The suite is a collection of independent
  `testSuite/test-*.py` (and `validate-*.py`) scripts. **Run only the handful
  that exercise what you changed** — pick the scripts whose names match the
  physics/subsystem you touched and run them directly:
  ```bash
  python3 testSuite/test-Python-interface.py   # example — choose ones relevant to your change
  ```
  Logs are written to `testSuite/outputs/`. Avoid `python3 testSuite/test-all.py`
  locally: it runs *every* script and takes a very long time (CI only manages it
  by sharding across many parallel jobs).
- **How pass/fail is judged:** tests and CI run a model, then grep the run log
  for failure markers (`FAIL`, `FAILED`, `fatal`, `aborted`,
  `ODE integration failed`, `unrecognized parameter`). A test "passes" when none
  of those appear — so when triaging a test, read its log, don't just trust the
  exit code.
- **Test scripts must exit `0` even when the test fails.** Failure is signaled
  *in the output*, not in the exit status: print a line containing `FAILED`
  (conventionally `FAILED: <what went wrong>`, with `SUCCESS: <what passed>` on
  the happy path) and then `sys.exit(0)`. The harness
  (`testSuite/test-all.py:42`) and the CI workflows judge each test with
  `grep -q -e FAIL -e FAILED` over the captured log and ignore the return code
  entirely, so a script that exits non-zero on failure is still recorded as a
  *pass* if it never printed the marker. Two corollaries:
  - Guard early-exit paths the same way — a missing input file or a failed model
    run must `print("FAILED: …")` before exiting, or the test silently vanishes.
  - The grep matches `FAIL` as a *substring anywhere* in the log, so never emit
    it in a passing context (no `"0 FAILures"`, no echoing a parameter named
    `...FAIL...`); that alone turns a green test red.

## Editing Fortran source

- **Never run a whole-file Fortran formatter** (`findent`, `fprettify`, or an
  editor "format document" action) on Galacticus `*.F90` sources. The build's
  code generation reads directive blocks embedded in comments — `!![ … !!]` (XML
  directives) and `!!{ … !!}` (reStructuredText docs) — and relies on `!$`
  OpenMP sentinels; whole-file formatters corrupt those blocks and relocate the
  sentinels, silently breaking the generated code. Edit by hand and match the
  surrounding style. (See
  [`docs/manuals/developer-guide/editor-setup.rst`](../docs/manuals/developer-guide/editor-setup.rst)
  and [`coding.rst`](../docs/manuals/developer-guide/coding.rst).)

## Spelling

- **Use US English spelling throughout** — source code (identifiers, variable
  and procedure names), comments, embedded `!![ … !!]`/`!!{ … !!}` directive
  blocks, documentation, and commit messages. So `color`, `normalization`,
  `analyze`, `center`, `behavior` — not `colour`, `normalisation`, `analyse`,
  `centre`, `behaviour`. This matches the existing API (e.g. `...Normalization`
  class names and parameters), so a British spelling in an identifier is not
  merely a style slip — it breaks consistency with names users write in their
  parameter files. The documentation build enforces this: `docs/conf.py` sets
  `spelling_lang = 'en_US'`, and the *Spell-Check-RST* CI job reports
  misspellings on PRs (add legitimate new words to `aux/words.dict`).

## Developer tools (`galacticusDevTools`)

A companion repository,
[`galacticusorg/galacticusDevTools`](https://github.com/galacticusorg/galacticusDevTools),
holds maintenance/development utilities that are *not* part of the model. It is
a separate checkout and is **optional** — never assume it is present; check
first, and either clone it or fall back to doing the job by hand if it is not.
Its `ReadMe.md` is the authoritative description of each tool; the ones most
likely to be useful when working in this repo:

- **`GPLerize.py <sourceDir>`** — adds or refreshes the standard Galacticus GPL
  header on every Fortran (`.f`, `.f90`, `.inc`) and C/C++ (`.c`, `.cpp`, `.h`)
  file in a directory, with the copyright year range generated dynamically. Use
  it instead of hand-copying headers into new source files (see *Commits &
  attribution* below). It rewrites files in place and leaves `~`-suffixed
  backups — so run it on a clean tree, and remember to delete the backups.
- **`migrateAllParameterFiles.py`** — run from the root of a Galacticus checkout
  (with `GALACTICUS_EXEC_PATH` set); walks `parameters/`, `constraints/`, and
  `testSuite/` and runs `scripts/aux/parametersMigrate.py` in place on every XML
  parameter file. Use this after a change that renames/restructures parameters,
  rather than editing the bundled parameter files one at a time. It also resets
  the `lastModified revision` in `testSuite/.../strictOutdated.xml` and
  `unstrictOutdated.xml`, which intentionally exercise the "outdated parameter
  file" paths.
- **`deltaTestCaseReducer/delta.sh`** — wrapper around the
  [Delta](https://github.com/dsw/delta) debugging tool; reduces a source file to
  a minimal case that still reproduces an error (compiler ICE, runtime crash,
  …) by successively deleting lines:
  ```
  ./delta.sh -test=testScript.sh -suffix=.F90 file.F90
  ```
  where `testScript.sh` exits `0` **iff** the error still occurs (see
  `testScriptExample.sh`). Reach for this when a compiler bug or crash needs a
  minimal reproducer for an upstream report — not for routine debugging, as each
  iteration recompiles and the reduction is slow.
- **`runBenchmarks.sh -e <exe> [-e <exe> …] [-r <repeats>]`** — runs executables
  alternately, `taskset`-pinned to a single CPU, and (given `sudo`) pins the
  governor to `performance` with turbo disabled, to make micro-benchmark timings
  comparable. Use it for before/after timings of a performance change; do not
  use it for ordinary test runs. Note it hard-codes CPU 3 and touches system
  CPU-frequency settings via `sudo`, so confirm with the user before running it.
- **`retrieveGHPagesArtifacts.sh <runID>`** — pulls validation, benchmark, and
  build-profile artifacts from a CI/CD run into a local `gh-pages` checkout so a
  PR's metric pages can be inspected before merge. Requires the `gh` CLI, and
  must be run from a directory with Galacticus' `gh-pages` branch checked out
  (not this working tree). It opens pages in a browser at the end.

Two further tools are for data/model maintenance rather than day-to-day work:
`extractSDSSBPTData.py` (rebuilds the SDSS DR8 emission-line constraint HDF5
datasets under `GALACTICUS_DATA_PATH`) and `promptCusps.py` (independent Python
reference values for validating `source/tests.prompt_cusps.F90`).

## Commits & attribution

- **Conventional Commits are enforced** by the `commit-msg` git hook (from the
  separate [`gitHooks`](https://github.com/galacticusorg/gitHooks) repo). Use
  `type(scope): summary`; allowed types: `fix`, `feat`, `build`, `chore`, `ci`,
  `docs`, `style`, `test`, `refactor`, `perf`, `revert`, `clean`. For a breaking
  change append `!` (e.g. `feat!:`) and add a `BREAKING CHANGE:` footer.
- **New source files need the GPL v3 header** (copy the block from any existing
  `source/**/*.F90`). Galacticus is licensed under GPL-3.0 (`COPYING`).
- **Attribution:** add inline `!+` markers to record contributors (auto-extracted
  by `scripts/doc/extractContributors.py`), and note AI assistance where used —
  see the *Contributor Attribution* and *AI-Assisted Contributions* sections of
  [`CONTRIBUTING.md`](../CONTRIBUTING.md).
