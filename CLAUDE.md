# CLAUDE.md

Guidance for [Claude Code](https://docs.claude.com/en/docs/claude-code) (and
other AI assistants) working in this repository. This is a short operational
quick-reference; it does **not** restate the full documentation. For anything
beyond the essentials below, defer to:

- [`README.md`](README.md) — quickstart, prerequisites, troubleshooting.
- [`CONTRIBUTING.md`](CONTRIBUTING.md) — contribution workflow, commit format,
  attribution, and the AI-assisted-contribution policy (**all AI-generated work
  must be verified by a human before submission**).
- The [developer guide](https://galacticus.readthedocs.io/en/latest/manuals/developer-guide/index.html)
  (`docs/manuals/developer-guide/`) — the authoritative reference for coding
  style ([`coding.rst`](docs/manuals/developer-guide/coding.rst)), the
  preprocessor-directive / `functionClass` code-generation system, and the
  workflow for
  [adding a new class](docs/manuals/developer-guide/creating-a-new-class.rst).

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

## Editing Fortran source

- **Never run a whole-file Fortran formatter** (`findent`, `fprettify`, or an
  editor "format document" action) on Galacticus `*.F90` sources. The build's
  code generation reads directive blocks embedded in comments — `!![ … !!]` (XML
  directives) and `!!{ … !!}` (reStructuredText docs) — and relies on `!$`
  OpenMP sentinels; whole-file formatters corrupt those blocks and relocate the
  sentinels, silently breaking the generated code. Edit by hand and match the
  surrounding style. (See
  [`docs/manuals/developer-guide/editor-setup.rst`](docs/manuals/developer-guide/editor-setup.rst)
  and [`coding.rst`](docs/manuals/developer-guide/coding.rst).)

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
  [`CONTRIBUTING.md`](CONTRIBUTING.md).
