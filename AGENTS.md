# AGENTS.md

Guidance for AI coding agents (Claude Code, Codex, Copilot, Gemini CLI, Cursor,
…) working in this repository; humans may find it useful too. This is a short
operational quick-reference; it does **not** restate the full documentation.
For anything beyond the essentials below, defer to:

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

This file is always loaded, so it is kept short. Longer, situational knowledge
is in [`docs/agents/`](docs/agents/): read the matching file **before** the
task it names (see [Knowledge files](#knowledge-files) below). Claude Code
loads this file through `.claude/CLAUDE.md`, which imports it.

## Core rules

These cost the most when forgotten. The detail behind each is in the knowledge
files.

1. **Export `GALACTICUS_EXEC_PATH` in the same command** as every `make`,
   commit, build script or test (``export GALACTICUS_EXEC_PATH=`pwd` && …``).
   Shell state may not persist between an agent's commands, and a value left
   pointing at another checkout gives racing builds, silently stale compiled
   code, and audits of the wrong tree (see
   [`docs/agents/building.md`](docs/agents/building.md)).
2. **Update a branch by merging master in, never by rebasing.**
   `scripts/aux/migrations.xml` and the `lastModified` stamps in parameter files
   are keyed to commit hashes, and a rebase silently disables them. Give the
   merge a Conventional Commits message (`chore(merge): bring master into …`).
3. **Never write a full commit hash by hand.** Substitute it from
   `git rev-parse`; a hand-expanded prefix invents the rest of the hash.
4. **A commit that renames or removes a parameter needs a follow-up commit**
   registering the migration in `scripts/aux/migrations.xml` (see
   [`docs/agents/committing.md`](docs/agents/committing.md)).
5. **Commits and pushes can take minutes** because of the git hooks: run them
   in the background and check `git log` before retrying.
6. **When an error appears in code you did not touch, suspect the build tree
   first** (stale module lists, dependency files or objects; see
   [`docs/agents/building.md`](docs/agents/building.md)), and include
   `mv: cannot stat` in every build-failure grep.
7. **Don't run `scripts/aux/formatModuleUses.py` or
   `scripts/aux/formatDeclarations.py`.** This is current project policy: the
   formatters are not yet considered ready, and the maintainer runs them when
   needed. Match the surrounding alignment by hand. The hook and CI formatting
   checks are advisory, and follow-up `style(...)` alignment commits are normal.

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
  generated `*_class` units are the peak. USE consolidation (`CONSOLIDATE_USES=yes`,
  the default) keeps each compile under ~3 GB, so `-j4` peaks at ~6 GB, and CI
  builds at **`-j4`**. With `CONSOLIDATE_USES=no` the largest compile needs
  ~14 GB, and `-j4` can OOM (exit 137) on a ~16 GB machine. On memory-constrained
  machines use a lower `-jN` and/or `LTO=disabled` (e.g.
  `make -j2 LTO=disabled Galacticus.exe`).
- **One `make` per tree.** Never run two `make` invocations concurrently in the
  same working tree — they race over generated files.
- Because a full build takes a long time, **prefer launching it in the
  background** and continuing other work, rather than blocking on it.

## Waiting on long builds & models (avoid the `pgrep` trap)

Builds and models can run for a long time, so it is tempting to poll with
`pgrep`. A naive `pgrep -f <pattern>` frequently **matches the poller's own
command line** (the shell/command running the `pgrep`), so it looks like the job
is "still running" forever, or matches unrelated processes. Prefer, in order:

1. **Let the harness manage it.** If your agent can run a command in the
   background and be notified when it exits, use that instead of busy-polling.
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
- **Test scripts must exit `0` even when the test fails.** The outcome is
  signaled *in the output*, by a marker line at the start of the line:
  - `FAILED: <what went wrong>` — the test failed;
  - `SUCCESS: <what passed>` — the test passed;
  - `SKIPPED: <why>` — a prerequisite could not be met (too few MPI processes, a
    non-debugging build, a missing executable), so the test ran nothing.

  Then `sys.exit(0)`. Both harnesses — `testSuite/test-all.py` and the CI
  workflows (`testModel.yml`, and the inline steps in `cicd.yml`, all of which
  run under `set -o pipefail`) — judge a test by grepping the captured log for
  `FAIL`.
  **A non-zero exit is also a failure to both**, so an uncaught exception *is*
  detected — but as a traceback rather than as a readable message. The point of
  exiting `0` is that the message, not the exit code, is what a reader sees; the
  one thing that genuinely disappears is a failure path that exits `0` without
  printing a marker. Four corollaries:
  - Guard early-exit paths the same way — a missing input file or a failed model
    run must `print("FAILED: …")` before exiting, or the test silently vanishes.
  - The grep matches `FAIL` as a *substring anywhere* in the log, so never emit
    it in a passing context (no `"0 FAILures"`, no echoing a parameter named
    `...FAIL...`); that alone turns a green test red. Keep `FAILED` at the start
    of the line for the same reason.
  - Use the canonical spellings. `FAIL:`, `FAILURE` and `PASS:` all happen to
    work today, but `scripts/aux/auditTestMarkers.py` and any future per-test
    reporting key on one spelling.
  - Reserve `SKIPPED` for a *whole-test* skip. A sub-check that is skipped (an
    optional validator that is not installed, say) should say so in plain words
    without the marker word, or `test-all.py` will report the whole test as
    skipped.
- **Check test scripts with `scripts/aux/auditTestMarkers.py`** after adding or
  editing one. It audits every `testSuite/test-*.py` and `testSuite/validate-*.py`
  against the convention above in under a second, and the *Validate-Test-Markers*
  PR check runs it with `--check`. Unlike the source formatters this check is
  blocking: the suite is fully compliant, so any report is new drift.
  ```bash
  ./scripts/aux/auditTestMarkers.py           # report
  ./scripts/aux/auditTestMarkers.py --check   # exit 1 if anything is reported
  ```
  A string literal that must contain a marker word without being one (a `grep`
  pattern, say) is exempted by ending its line with `# markers: exempt`.
- **Three scripts are not run by any CI workflow** and are run manually (or via
  `test-all.py`) only: `testSuite/test-regressions.py` (superseded by the
  `testSuite/regressions/` matrix), `testSuite/validate-concentrationConvergence.py`
  and `testSuite/validate-mergerTreeConvergence.py` (both are long validation
  runs that submit jobs through a queue manager).

## Before you push

Run the checks that CI enforces but the git hooks miss (each has cost a red
PR; see [`docs/agents/committing.md`](docs/agents/committing.md)):

```bash
export GALACTICUS_EXEC_PATH=`pwd`
python3 scripts/aux/docBlockAudit.py --check     # unknown :cite: keys, doc blocks
python3 scripts/aux/auditRevisions.py --check    # lastModified stamps are real ancestors
git diff --name-only --diff-filter=AM origin/master -- '*.xml' | xargs -r grep -nP '[^\x00-\x7F]'  # non-ASCII in XML
```

- For a **new source file**, copy it into an empty `source/` directory somewhere
  temporary and run `python3 scripts/doc/convertDocstringsToRST.py --check` on
  that directory (every `<methods>`-style directive needs `docformat="rst"`).
  Never run it without `--check` on the checkout: it rewrites files.
- Grep the build log for `Warning:` as well as errors: CI fails on unsuppressed
  compiler warnings.
- For a **build-system change**, run the whole-repository `python -m pytest`
  from the repository root, as CI does.

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
- Before adding a class, parameter or directive, or changing code that copies,
  assigns or privatizes objects or uses OpenMP, read
  [`docs/agents/fortran.md`](docs/agents/fortran.md).

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
Its tools are grouped into directories by purpose (`maintenance/`,
`debugging/`, `codeAudits/`, `dataPreparation/`, `physicsReferences/`). The
top-level `ReadMe.md` indexes them, and each directory's `ReadMe.md` is the
authoritative description of its tools. Paths below are relative to the
repository root. The ones most likely to be useful when working in this repo:

- **`maintenance/GPLerize.py <sourceDir>`** — adds or refreshes the standard Galacticus GPL
  header on every Fortran (`.f`, `.f90`, `.inc`) and C/C++ (`.c`, `.cpp`, `.h`)
  file in a directory, with the copyright year range generated dynamically. Use
  it instead of hand-copying headers into new source files (see *Commits &
  attribution* below). It rewrites files in place and leaves `~`-suffixed
  backups — so run it on a clean tree, and remember to delete the backups.
- **`maintenance/migrateAllParameterFiles.py`** — run from the root of a Galacticus checkout
  (with `GALACTICUS_EXEC_PATH` set); runs `scripts/aux/parametersMigrate.py` in
  place on the root `parameters.xml` and on every XML parameter file under
  `parameters/`, `constraints/`, and `testSuite/`. Use this after a change that
  renames/restructures parameters, rather than editing the bundled parameter
  files one at a time. It skips `testSuite/.../strictOutdated.xml`,
  `unstrictOutdated.xml`, and `parameterMigration.xml`, which must stay out of
  date to exercise the "outdated parameter file" and migration paths - do not
  migrate them by hand either.
- **`debugging/deltaTestCaseReducer/delta.sh`** — wrapper around the
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
- **`debugging/runBenchmarks.sh -e <exe> [-e <exe> …] [-r <repeats>]`** — runs executables
  alternately, `taskset`-pinned to a single CPU, and (given `sudo`) pins the
  governor to `performance` with turbo disabled, to make micro-benchmark timings
  comparable. Use it for before/after timings of a performance change; do not
  use it for ordinary test runs. Note it hard-codes CPU 3 and touches system
  CPU-frequency settings via `sudo`, so confirm with the user before running it.
- **`maintenance/retrieveGHPagesArtifacts.sh <runID>`** — pulls validation, benchmark, and
  build-profile artifacts from a CI/CD run into a local `gh-pages` checkout so a
  PR's metric pages can be inspected before merge. Requires the `gh` CLI, and
  must be run from a directory with Galacticus' `gh-pages` branch checked out
  (not this working tree). It opens pages in a browser at the end.
- **`codeAudits/auditLoop/`** — tooling for a scheduled, unattended code audit.
  It provides a `galacticus-audit` Claude Code skill, which audits prioritized
  source files and stages verified fixes on local `audit/*` branches. It also
  provides scripts that track findings and review coverage in a state
  directory, rank files for the next audit, and gate runs on plan usage. If a
  checkout has `audit/*` branches, they belong to that loop: leave them alone
  unless the user asks. Its `ReadMe.md` describes the setup and the formats.

Two further tools are for data/model maintenance rather than day-to-day work:
`dataPreparation/extractSDSSBPTData.py` (rebuilds the SDSS DR8 emission-line
constraint HDF5 datasets under `GALACTICUS_DATA_PATH`) and
`physicsReferences/promptCusps.py` (independent Python reference values for
validating `source/tests.prompt_cusps.F90`). `physicsReferences/` also holds
the independent reference scripts that many unit tests in `source/tests/` cite
by name.

## Knowledge files

[`docs/agents/`](docs/agents/) holds traps found in practice, each with the
date it was last verified. Read the matching file before the task:

| File | Read before |
|---|---|
| [`building.md`](docs/agents/building.md) | debugging a build failure; build options (MPI, profilers, library); editing the `Makefile` |
| [`fortran.md`](docs/agents/fortran.md) | adding a class, parameter or directive; changing object copying, assignment, OpenMP or error handling |
| [`committing.md`](docs/agents/committing.md) | committing, updating a branch, changing a parameter, pushing |
| [`testing.md`](docs/agents/testing.md) | writing or trusting a test; comparing model outputs; benchmarking or profiling |
| [`ci.md`](docs/agents/ci.md) | diagnosing a CI failure; editing a workflow; working with stacked PRs |

A knowledge entry records what was true when it was written: if it names a
file, flag or command, check that it still exists before relying on it.

### Skills

[`.agents/skills/`](.agents/skills/) holds shared agent skills in the Agent
Skills format (`SKILL.md` plus supporting files); `.claude/skills` links to it
for Claude Code. `galacticus-review` reviews a PR, a branch, or everything
merged to master since the last such review, for physics, numerical
correctness and bugs.

### Proposing a change to this guidance

These files, and this one, are shared by everyone's agents, so changes go
through a pull request that a maintainer reviews: never commit a change to them
directly to master. When you learn something that would have saved you real
time and that holds in any checkout, propose it as an entry in the matching
knowledge file (or, if it is short and needed in most sessions, as a rule
here), in this form:

```markdown
## <Short title stating the trap>
*Added YYYY-MM-DD; verified YYYY-MM-DD*

<The fact, with the evidence: dates, PR numbers, commit hashes.>

**Why:** <what goes wrong if you don't know it>

**How to apply:** <what to do>
```

Search the existing files first, and update an entry rather than adding a
near-duplicate. Leave out anything specific to one person or machine (names,
host names, home-directory paths, local checkout names). If an entry is wrong
or stale, propose fixing or deleting it. If a trap could be removed by fixing
the code, a hook or a workflow, say so in the PR: a fix is better than an
entry.

## Commits & attribution

- **Conventional Commits are enforced** by the `commit-msg` git hook (from the
  separate [`gitHooks`](https://github.com/galacticusorg/gitHooks) repo). Use
  `type(scope): summary`; allowed types: `fix`, `feat`, `build`, `chore`, `ci`,
  `docs`, `style`, `test`, `refactor`, `perf`, `revert`, `clean`. For a breaking
  change append `!` (e.g. `feat!:`) and add a `BREAKING CHANGE:` footer.
- **Install the git hooks** as `CONTRIBUTING.md` describes. Two of them (an
  automated review of the staged diff and a suggested commit message) call the
  `claude` CLI and are **opt-in**, since they spend the committer's own usage:
  enable them with `GALACTICUS_CLAUDE_HOOKS=1` or
  `git config galacticus.claudeHooks true`. See the
  [`gitHooks` README](https://github.com/galacticusorg/gitHooks).
- **New source files need the GPL v3 header** (copy the block from any existing
  `source/**/*.F90`). Galacticus is licensed under GPL-3.0 (`COPYING`).
- **Attribution:** add inline `!+` markers to record contributors (auto-extracted
  by `scripts/doc/extractContributors.py`) — see the *Contributor Attribution*
  and *AI-Assisted Contributions* sections of
  [`CONTRIBUTING.md`](CONTRIBUTING.md). A marker records **names only**, in
  the canonical form
  ```fortran
  !+    Contributions to this file made by: Andrew Benson, Claude.
  ```
  Where an AI tool contributed, add its bare name (`Claude`, `Copilot`, …) to
  that list — **never** a description of what it contributed, and never a
  sentence such as "drafted with assistance from Claude and verified by …".
  The marker is parsed as a comma-separated list of names, so prose leaks into
  the generated contributor list in the documentation's Acknowledgments, and
  makes the in-source line needlessly long. Put the description of *what* the AI
  contributed, and how it was verified, in the **commit message** instead.
  `python3 scripts/doc/extractContributors.py --check` enforces this (as does
  the *Validate-Contributor-Markers* CI job); run it after editing a marker.
