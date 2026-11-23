# Building

*Read this before debugging a build failure, changing build options, or editing the `Makefile`.* It covers how a
build tree goes stale or inconsistent and how to recover it, circular module dependencies, the build-option
variables, and the GNU make 3.81 constraint. These entries record traps found in practice; each names when it was
last verified, so check that a named file, flag or command still exists before relying on it.

## What goes wrong when `GALACTICUS_EXEC_PATH` points at another checkout
*Added 2026-06-24; verified 2026-10-30*

`AGENTS.md` says to export `GALACTICUS_EXEC_PATH` in the same command as every `make`, commit, build script or test.
Code generation, the pre-commit hooks, and `test-python-utils.py` all read it, and shell state does not persist
between an agent's separate shell calls. A shell profile that sets it to one checkout silently applies it to every
other checkout and worktree too. When it points at the wrong tree:

- parallel targets each run `rm _class.p.F90; make _class.p.F90` on the same file, causing "multiple recursive
  makes" and `No such file` fatal errors;
- the pre-commit test-marker audit checks the *other* tree's `testSuite/` and blocks the commit;
- XML schema validation uses a stale schema and gives false errors;
- Python tooling imports the other tree's generator;
- `docBlockAudit.py` and `auditTestMarkers.py` audit `$GALACTICUS_EXEC_PATH`, not the current directory, so from a
  worktree they silently report on another tree (2026-10-29: eleven pre-push checks in a worktree audited the main
  checkout; one real problem surfaced only once the path was set);
- `make` can succeed and still compile **stale** code. The `functionClass` preprocessor places each member's
  `.p.F90` relative to `$GALACTICUS_EXEC_PATH/source`. For a source outside that tree it writes
  `work/build/<basename>.p.F90` instead, and the real one stays stale (2026-10-30: a changed `universe/operators`
  member compiled as master's code, and its test then did not pass). If stray `work/build/*.p.F90` files other than
  `Galacticus.p.F90` and `state.p.F90` appear, this has happened: delete them, touch the changed members, and
  rebuild with the path set.

**Why:** each of these looks like a problem with the change being built or committed.

**How to apply:** prefix every `make`, `git commit`, build script and test with ``export GALACTICUS_EXEC_PATH=`pwd` ``
(plus `PYTHONPATH=$(pwd)/python` for manual script runs). After a race like this, delete the partial build
directory and rebuild.

## Worktree builds of older branches print a revision-check warning
*Added 2026-10-21; verified 2026-10-28*

A branch that predates PR #1686 (merge `c81df3a6d4`, 2026-10-28) records `gitHash="unknown"` when built in a
`git worktree`, so every run prints `WARNING: parameter file revision check failed (#2; error code: 2)`: ignore it,
or merge master in when the recorded hash matters.

## Recovering a stale or inconsistent build tree
*Added 2026-07-02; verified 2026-10-23*

- **After merging master (the stale `.m`/`.smod` trap):** a commit which renames modules - e.g. `d8a194faa0`, which
  expanded vowel-stripped implementation names - leaves `work/build/**/*.m` module lists naming `.smod` files the
  compiler no longer produces. The `%.o` recipe cats the `.m` and `mv`s each name it lists out of `moduleBuild/`,
  so the build dies with `mv: cannot stat ...smod` and *no* compiler diagnostic. The renaming commit may be an
  ancestor of the checkout's **old** head, so the residue can predate the merge that exposes it - the merge merely
  forces the first recompile of an affected unit (2026-09-17: residue from before 2026-08-17).
  **The test is only valid on a quiescent tree, and `cannot stat` must be in every build-failure grep.** A `.m`
  regenerated *during* the build (its source changed in the merged commits) can develop the mismatch after a
  pre-build check reported zero: on 2026-09-17 one checkout checked clean, then failed on exactly the unit another
  checkout had failed on. So run the check again *after* a failure, not only before. Three times that day a failure
  grep of `Error:|Fatal|undefined reference` printed nothing because the real marker was `mv: cannot stat`.
  **The exact test is each `.m` against its own recipe.** `work/build/Makefile_Module_Dependencies` holds an
  explicit rule per module-producing unit - `<path>.m:` followed by `@echo <name> > <path>.m` (1843 of them on
  2026-09-17) - so read the expected names from those recipes and delete every `.m` whose contents differ. Those
  regenerate *correctly*; the "regenerates empty" risk applies only to units with **no** explicit rule, which fall
  through to the touch-only `$(BUILDPATH)/%.m : ./source/%.F90` pattern rule in the `Makefile`. Do *not* instead
  compare against "modules the rules produce": that heuristic reported 2 files in one checkout and missed the very
  unit whose `mv` had failed. Then delete the `.mod`/`.smod` on disk that no recipe or target names **and** whose
  name is vowel-stripped (a run of 5+ consonants) - 17 in one checkout, matching its 17 mismatched `.m`, and
  matching another checkout's 17 on 2026-09-16. Leave every other unknown file alone: 40 in that checkout
  (`galacticus_nodes.smod`, `events_hooks.smod`, node-component smods) are needed yet named by no rule, and 226 `.m`
  files list a parent-module `.smod` that is no target. Both are false positives, so a bare "not a target" sweep
  destroys working files.
- **After moving a module out of a file into a new one** (e.g. splitting a test's helper module from its program):
  the old file now has no explicit `.m` rule, and the touch-only fallback never empties its `.m`, so the old unit's
  compile still `mv`s the moved module and dies with `mv: cannot stat ...mod`. Delete that one `.m`; it is recreated
  empty (2026-10-23).
- **After switching branches:** leftover `work/build` artifacts make make fail with "No rule to make target".
  Cached directive and dependency files also embed the file list. Either run
  `grep -rl <name> work/build | xargs rm -f`, `rm -rf` the affected subdirectory and `mkdir -p` it back, or run
  `rm -rf work/build && make`.
  **Building an older base in a tree last built on newer master, then moving forward again:** the old-base build
  rewrites `.d` files without master's new `use`s, and they then look newer than their prerequisites, so the
  forward build links with `undefined reference` to a module that compiled fine (2026-11-09). Delete the `*.d`
  files that build wrote (`find work/build -name '*.d' -newermt <start> ! -newermt <end>`), or don't reuse the tree.
  The same link failure follows simply **adding a `use` of an already-built module** to a program whose `.d` is
  newer than that module's `.mod.d`: make never regenerates the `.d` (2026-11-09,
  `tests.dust.emission_spectra`). Delete that program's `.d` and rebuild.
- **After adding a new functionClass implementation** - including an *abstract intermediate* - the dependency
  blobs are keyed on `_class.F90`'s mtime and miss the new file, so linking fails with `undefined reference`. Run
  `rm work/build/Makefile_{Module,Use}_Dependencies{,.blob}`.
  Second symptom, **fixed by PR #1666 (merged 2026-10-23, `_BLOB_VERSION` 6)**: editing a `use` in an existing
  implementation, or a `functionGlobal`/static-hook file, left stale rules ("No rule to make target
  `iso_c_binding.mod`"). In a tree whose `useDependencies.py` predates 4f48b5b7f7, `touch <dir>/_class.F90` still
  works around it. The first symptom (link-time `undefined reference`) involves the Module_Dependencies blob and is
  not covered by that fix. Also run `make parameters-schema` and stage `schema/parameters.xsd` (the pre-commit hook
  enforces this): it *does* change, gaining one `<xs:enumeration>` entry, because every abstract class is
  registered too - all 32 `abstract="yes"` classes appear there (2026-09-17).
  The hook generates the schema from the *working tree*, not the index. So while such a file is still untracked,
  every earlier commit in a split series is rejected with "Parameter schema is out of date". Move the file aside
  while making those commits (2026-09-22). General fix for any split series (confirmed 2026-10-14): stage a
  commit's files, then `git stash push --keep-index --include-untracked`, commit, `git stash pop`.
- **After adding a `<method>` to a functionClass:** you get "Mismatch in components of derived type" (a stale
  vtable). Delete the stale `.o` files but *not* the `.mod` files:
  `find work/build -name "*.o" ! -newermt "<change date>" -delete`. Then repeat `make -k` until it exits 0.
- **After editing source:** make may not relink. Run `rm -f Galacticus.exe` or `rm -f libgalacticus.so` first, then
  confirm the binary's mtime advanced. A compile that does not succeed leaves the old binary in place.
  Restoring a source file after a temporary edit (a mutation test, say) with `cp -p` gives it an mtime *older* than
  the object built from the edit, so make silently keeps the edited object. `touch` the restored file (2026-10-20).
- **Never hand-compile objects into `work/build`** with different flags (`-O0`, `-flto=auto`). The resulting
  `.mod`/`.smod` mismatch makes later submodule compiles crash gfortran. Use `make work/build/<path>.o`.
- **Never delete `*.p.F90` to force regeneration.** Since `405868dcf5` (2026-08-17) the `%.p.F90.up` sentinel
  depends on `preprocessorSources.digest`, so editing a generator `.py` re-preprocesses every file by itself.
  `%.p.F90`'s own recipe is just `@true`: a deleted `.p.F90` whose `.up` is current is never recreated, and the
  compile then aborts with `f951: Fatal Error: ...p.F90: No such file or directory` (2026-09-19). To force one file,
  `touch` its source instead. Deleting the `.up` too doesn't help: a `functionClass` implementation's `.p.F90` is
  written by preprocessing its `_class.F90` (hand-running `preprocess.py` on the implementation gives an invalid
  non-module file), and a missing `.up` can be satisfied by the `%.up: touch` fallback. Once a tree has been mangled
  this way, `rm -rf work/build<TAG>` and rebuild (2026-10-12).
- Some bugs reproduce only with production `-O3 -flto`.

**Why:** each of these symptoms looks like a bug in the code being changed.

**How to apply:** suspect the build tree first when an error appears in code you didn't touch.

## gfortran crashes and circular module dependencies
*Added 2026-07-22; verified 2026-09-14*

A gfortran internal error while compiling a *submodule* (a segfault in `gfc_get_derived_type`), hitting a different
file on each build, means make dropped a circular dependency. Grep the build log for
`Circular ... dependency dropped`; a clean tree emits none.

Every functionClass implementation is a submodule of its class module. So if any implementation of class A `use`s
class B, module A depends on all of B, and two classes can never reference each other.

**Why:** it looks like a compiler bug, but it's an ordering problem.

**How to apply:**
- Check the direction of use before adding a dependency. When a high-level object is needed, pass it as an
  optional argument, or use the `functionGlobal` pattern in a `utilities.F90` module (see
  `source/nodes/operators/utilities.F90`).
- Module-scope declarations in an implementation file are hoisted into the parent class module, so sibling
  implementations see them without a `use`.
- Procedures after `contains` in a `_class.F90` that implementations call must be `public`. Implementations compile
  to separate objects, and a private module procedure is not exported, so the link stops with
  `undefined reference to __<module>_MOD_<name>` although compiling succeeds (2026-09-14).

## Build options (MPI, profilers, debugging, library)
*Added 2026-07-10; verified 2026-10-17*

Since PR #1561 (merged 2026-10-09) each build variant is its own make option, and they combine: `MPI=yes`,
`LIBRARY=yes`, `PROFILER=gprof|perf`, `ODEPROFILE=yes`, `DEBUGGING=yes`, `BUILDPROFILE=yes`. Options that change the
code add a tag to both the build directory and the executable suffix: `MPI=yes PROFILER=perf` builds in
`work/buildMPIPerf` and produces `Galacticus.exe_MPIPerf`; `MPI=yes` alone gives `Galacticus.exe_MPI` and no longer
overwrites the serial binary. `SUFFIX=` (empty) restores plain `Galacticus.exe`, e.g. for testSuite scripts, which
call `./Galacticus.exe`. `GALACTICUS_BUILD_OPTION` still works with a deprecation warning and its old names (`MPI`
then still overwrites `Galacticus.exe`).

- **MPI:** `MPI=yes` alone still uses the exported `FCCOMPILER` (e.g. `gfortran-16`), so pass the MPI wrappers
  explicitly: `MPIFCCOMPILER=mpif90 MPICCOMPILER=mpicc MPICPPCOMPILER=mpic++` (with whatever `PATH` your machine's
  MPI installation needs). Run with `mpirun --oversubscribe -np N`.
- **The make target stays `Galacticus.exe`** (`make MPI=yes Galacticus.exe`); the suffix is added for you.
  `make MPI=yes Galacticus.exe_MPI` finds the existing file with no rule, prints "Nothing to be done" and exits 0,
  leaving a stale binary (2026-10-17). Check the binary's timestamp after any build.
- **Any compiler or flag change rebuilds everything** in that build directory: `$(BUILDPATH)/buildConfiguration.stamp`
  records them. Expect a full rebuild after changing `LTO`, `GALACTICUS_FCFLAGS`, etc. `make -n` always lists a full
  rebuild (the only-if-changed `.d` rules), so it can't be used to preview one.

**Why:** the old single-value option gave misleading builds (silently non-MPI, overwritten executables, stale
objects after flag changes).

**How to apply:** use the new option names, and pass the `MPI*COMPILER` variables for MPI builds.

## The Makefile must work with GNU make 3.81, which Linux cannot reproduce
*Added 2026-10-08; verified 2026-10-09*

The macOS CI jobs run Apple's GNU make 3.81 (make's ``` `quoted' ``` message style in the log gives it away), so
make 4 features such as `$(file ...)` must not be used: in 3.81 an unknown function silently expands to nothing. A
make 3.81 built on Linux - stock (needs `glob/glob.c` `==`→`>=` on `_GNU_GLOB_INTERFACE_VERSION`,
`CC="gcc-16 -std=gnu89 -fcommon"`) or Apple's `apple-oss-distributions/gnumake` tag `gnumake-136` (also needs
`minus_c_len` declared off macOS) - hangs under `-j` at the same early point on *master* as well, in its own update
loop with no child processes. So it cannot reproduce a macOS-only `-j` hang. At `-j1` it builds normally, which is
enough to check that 3.81 parses a change.

**Why:** on 2026-10-08 hours went into "reproducing" a PR's macOS CI hang with Linux-built 3.81 before a master
control showed the same hang.

**Confirmed cause of one such hang (PR #1561, 2026-10-09):** a stamp written by a rule run on every build (`FORCE`
prerequisite, with a target-specific `export`), which thousands of targets depended on, hung macOS make 3.81 under
`-j`. Writing the stamp at parse time instead, so that it is a plain file prerequisite with no rule, fixed it.

**How to apply:** always run master as the control first. For a macOS-only parallel problem, compare CI logs (PR
against master) and let macOS CI judge the fix.
