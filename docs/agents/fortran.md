# Fortran and code-generation pitfalls

*Read this before adding a class, parameter or directive, or editing Fortran that copies, assigns or privatizes
objects, opens an OpenMP region, or reports a numerical failure.* It covers the directive parser's traps, `<generic>`
templates, abstract intermediate classes, reference counting and OpenMP under gfortran, smaller code traps, and
warnings and GSL status handling. These entries record traps found in practice; each names when it was last
verified, so check that a named file, flag or command still exists before relying on it.

## Parser and naming traps
*Added 2026-08-09; verified 2026-09-13*

- **No `::` inside strings in declarations.** A declaration whose string initializer contains `::` breaks the
  declaration parser with `contains '::'`. Use a braced suffix such as `{defaulted}` instead.
- **Escape `<` and `>` in directives.** Inside a `!![ … !!]` directive, including RST roles with an explicit
  `<target>`, write `&lt;`/`&gt;`. Otherwise the XML parse breaks silently and the class doesn't exist.
- **String parameter defaults.** A `varying_string` input parameter defaulting to `var_str('')` aborts when written
  to output. Use `var_str('none')`.
- **Intrinsic modules.** `<modules><name>ISO_C_Binding</name></modules>` in a method emits a non-intrinsic `use`, so
  make wants `iso_c_binding.mod`. Use a module-level `use, intrinsic ::` instead. `use` isn't allowed inside
  `block`.
- **Registered names.** A functionClass registers under its type name minus the class prefix, and keeps a leading
  acronym: `NGenHalofit`, not `nGenHalofit`. Read the name from `schema/parameters.xsd` rather than guessing.
- **Sub-parameters inherit: object builders walk up the parameter tree.** Every `<objectBuilder>` expands to a loop
  `do while (.not.parametersCurrent%isPresent('<class>').and.associated(parametersCurrent%parent))`
  (see `python/Galacticus/Build/SourceTree/Process/ObjectBuilder.py`), so a class declared once at the top level is
  found from any nested scope, and `parameterNode%objectCreated()`/`objectGet()` makes those scopes share one
  instance. A default is built only when the parameter is absent all the way to the root. Do *not* conclude
  otherwise from `inputParameters%isPresent`, whose `searchInParents` argument defaults to `.false.` and is passed
  in just one place in the tree: the walk lives in the generated builder wrapped around it, not inside
  `isPresent`. So `id`/`idRef` is for reaching a declaration the walk would *not* find (a different branch, or a
  differently named parameter), not for making nesting work. The opposite was once concluded (2026-09-15) from
  reading `isPresent` alone; removing every `idRef` from a test's parameter file left the test passing unchanged.
- **`public ::` lines.** Keep all module `public ::` names on one line; the build silently drops earlier `public ::`
  lines.
- **A module-scope constant in a functionClass implementation file can't be exported.** Its declarations go into the
  implementation's *submodule*, but a `public :: name` line goes into the class module, which then exports an
  implicitly typed, uninitialized variable of that name (no `implicit none` there). It compiles and links silently;
  the importer reads 0 (2026-10-28: a `log(0)` SIGFPE). Duplicate the constant (with a unit test tying the copies)
  or put it in a plain module. Public *procedures* in implementation files (e.g. NFW's state store) do work.
- **63-character name limit, prefix included.** A functionClass method's generated default procedure is named
  `<className><MethodName>` (e.g. `kinematicsDistribution` + 44 characters = 66), and gfortran rejects it with a
  cascade beginning `Name at (1) is too long` in `_class.o`. Implementation function names
  (`kinematicsSphericalScaler...`) hit it too. Count before naming (2026-10-14).
- **functionClass `<argument>` on one line.** If a method's `<argument>` declaration continues onto a second line
  (`&amp;`), the generated default function statement does too. An argument on that continuation line escapes both
  the generator's unused-argument list and `postprocess.py`, which matches a warning to its procedure from the
  echoed source line. So `-Wunused-dummy-argument` can't be suppressed, even with `!$GLC attributes unused`, and CI
  rejects the build (#1493, 2026-09-15).
- **A wrong argument list in a `referenceConstruct` directive** surfaces as
  `Component '<x>' at (1) is a PRIVATE component of '<class>class'`, pointing at the generated
  `object_=constructor(...)`. It means no specific constructor matched, so the intrinsic structure constructor was
  tried - not a visibility problem. After changing a constructor's signature, find call sites with an unfiltered
  `grep -rn "<typeName> *("`.
- **The descriptor generator must match every `objectBuilder` target to a member.** Objects built into a linked list
  the class does not declare with `<linkedList>` abort with "could not find a matching internal object for object
  [...]", and a class may declare only one `<linkedList>`. Take a second list of children as a single object
  parameter (e.g. a `nodePropertyExtractorMulti`) instead.
- **functionClass base `<data scope="self">` members** (2026-10-11): the name must not clash with a `functionClass`
  binding (`isDefault`, `descriptor`, `objectType`, ...; see `source/objects/function_class.F90`), or the generated
  deepCopy/stateStore are rejected with "Expected argument list". And `<data scope="module">` is emitted *after* the
  generated type, so a parameter used in a self member's declaration (e.g. `dimension(n)`) must be declared before
  the `!![ <functionClass>` block. For array indices, an `<enumeration>` directive placed there (`<validator>yes`
  gives `<name>Count`; index with `%ID`) works. The schema fixes element order (`validator`, `visibility`,
  `indexing`), and a *private* enumeration warns that its generated `Describe` function is unused.
  An enumeration taken by a class's constructor must be `<visibility>public`: the `LIBRARY=yes` build compiles each
  implementation separately with a C wrapper that imports the type from the class module ("Symbol
  'enumeration...type' ... not found in module"). Only CI's Build-Library jobs catch it, 1.5-2.5 h in (PR #1713,
  2026-10-30).
- **Check build logs for warnings, not just errors.** CI rejects unsuppressed warnings that a local build only
  prints. Grep `Warning:` too, and remember gfortran quotes names with curly quotes (`‘name’`).

**Why:** most of these only show up at run time, or after a full build.

**How to apply:** check this list when adding a class, parameter or directive.

## `<generic>` conversions: mechanics, and how to verify one
*Added 2026-09-16; verified 2026-09-17*

Collapsing parallel families of procedures onto one `<generic>` template (done in `objects/history.F90`,
`utility/IO/HDF5/_module.F90`, `utility/MPI.F90`):

- **The identifier and every modifier name must be all-lowercase.** Declaration text is lowercased when the file is
  first parsed, *before* Generics runs, so a mixed-case tag such as `{Rank¦rank}` inside `dimension(...)` stops
  matching its identifier and is emitted verbatim: `Error: Expected expression in array specification`. Mixed case
  survives only in subunit openers and other non-declaration text - which is why `history.F90`'s
  `identifier="Type"` works.
- **`schema/generic.xsd` requires `label` *and* `intrinsic`** on every `<instance>`; any other attributes are allowed
  (`anyAttribute processContents="lax"`). The pre-commit hook validates the fragment, so give `intrinsic` its real
  value and use it, rather than adding a dead attribute.
- **`label=""` is legal**, and reproduces un-suffixed names byte-for-byte (keeps `History_*`).
- **Several `<generic>` roots may share one `!![ … !!]` block** (four in `IO/HDF5/_module.F90`), and a `<generic>`
  expands *every following sibling* at its level - so placing it just before `contains` reaches all the module's
  procedures.
- **`<instance>` attributes may be split across lines.** Values carrying alignment padding otherwise give
  400-character lines; the tree's pristine maximum is 256.
- **Padding inside attribute values is the only way to align the generated declarations**, since generated code
  cannot be hand-aligned afterwards. A declaration with no `dimension(...)` (`self`, `iError`) cannot be aligned
  across instances of differing rank at all.

**Verify by preprocessing the one file, not by rebuilding:**
`BUILDPATH=work/build GALACTICUS_EXEC_PATH=$(pwd) ./scripts/build/preprocess.py source/X.F90 out.p.F90` takes
seconds (it needs `BUILDPATH`, which the Makefile normally supplies). Compare each procedure against a copy saved
before the edit, collapsing whitespace *and* punctuation spacing and rewriting `[line NNNN]` to a constant - that
suffix is the only line-number-dependent text in the output. Never blanket-strip digits: the introspection
payload's function and module names, `dim=1`, and `MPI_Integer8` all matter. `cmp` is the wrong test, because
changing the directive block's line count shifts every `[line NNNN]`.

**This verification works only for plain modules.** For a **functionClass implementation** the build's
`work/build/**/*.p.F90` is the *post-split* form - `submodule (Nodes_Operators) …` with `module procedure` bodies and
declaration identifiers lowercased - while a standalone `preprocess.py` run emits the unsplit file with its full
type definition and directives. The two are different representations, so comparing them reports hundreds of bogus
token differences (`-submodule`/`+use`, `-module procedure`/`+subroutine`, `-darkmatterprofile`/`+darkMatterProfile`).
This happened on 2026-09-17 while trying to prove an edit was whitespace-only; for such files verify behaviorally
instead (rebuild, run a covering test, compare `readelf -sW` symbols).

**Why:** each of the traps above cost a broken build or a wrong verdict, and the full-`make` loop costs ~15 minutes
against seconds.

**How to apply:** write the template, preprocess, compare per procedure, and expect every surviving difference to
be one predicted in advance. Then build both configurations: when the module's public interface is unchanged, only
that one object recompiles (every other `.p.F90` regenerates identically, even though adding a directive changes
`directiveLocations.xml` and so re-runs `preprocess.py` on the whole tree), and the clean link of all the
*un*-recompiled dependents is itself proof that the interface held.

## Abstract intermediate classes in a functionClass hierarchy
*Added 2026-09-17; verified 2026-09-17*

Moved to the developer guide: *Sharing code between implementations: abstract intermediate classes* in
`docs/manuals/developer-guide/creating-a-new-class.rst`. Read it before sharing code between implementations of a
class.

## Reference-counted objects, OpenMP, and gfortran-16 bugs
*Added 2026-08-03; verified 2026-11-19*

- **An `!$omp parallel` region is expensive; don't open one per loop.** (2026-11-09.) The build injects event-hook
  filter calls (`eventsHooksFilterCopyIn` etc., `Process/EventHooks.py`) into every region, and they serialize on
  libgfortran's unit lock (internal writes in `concatenate_varstr_integer`). A solver with ~40k `parallel do` calls
  took 75 s on 8 threads and 293 s on 32, against 21 s serial; one region with orphaned `!$omp do`/`single` took 2 s
  on 32. **Inside an orphaned `do`, every thread runs the whole procedure:** whole-array zeroing before the `do`
  races with the loop, so write results only inside it. **A contained procedure reaches the host's original
  variable, not a thread's `private` copy:** pass private variables (e.g. the loop index of the region) as
  arguments. Reading `omega(k)` through host association gave garbage `k` on the workers, and a floating point
  exception deep in the cosmology tables (2026-11-10).
- **OpenMP `private` clauses.** gfortran doesn't default-initialize the private copy of a derived type, so garbage
  pointers in reference-counted types (`integrator`, `interpolator`, anything with a `resourceManager`) deadlock or
  segfault on first assignment. Declare them `allocatable` and allocate inside the parallel region. Plain locals of
  a called procedure are safe. Root cause of #1321. Known unfixed instances of this and the next pattern are tracked
  as bugs; grep for `private` clauses naming such objects before assuming a file is clean.
- **Assigning onto used objects.** Assigning a constructor result onto an already-used object that has a
  `resourceManager` component double-finalizes it and corrupts the heap (GCC PR 110626; still present in gfortran
  16.0.1 20260416, checked 2026-11-19). Use `deallocate`/`allocate` around the assignment, hoist the constructor out
  of the loop, or give the type a type-bound `assignment(=)` (as `interpolator` and `nearestNeighbors` have).
- **Interpolators.** They carry mutable state (a GSL accelerator, and a cached bracket index), so never share one
  between threads. Share only read-only tables. A class that needs to own one must declare it
  `type(interpolator), allocatable` (an array is fine) and *assign* into it: the code generator then emits
  `deepCopyActions()` (which calls `GSLReallocate`, so each deep copy gets its own GSL objects) and
  `stateStore`/`stateRestore` for it automatically, and `resourceManager` frees the GSL objects when the parent is
  destroyed. A `pointer` to one, including a pointer inside a wrapper type, gets *none* of that: the generated deep
  copy copies it by intrinsic assignment, so every per-thread copy shares one accelerator, and nothing ever
  deallocates it. `Filter_Response_Function` allocates a new interpolator and hands ownership to the caller, so its
  result must be copied into an allocatable and the original deallocated, as `Filter_Response` does. Fixed this way
  in `luminosity_stellar_from_SED.F90` (PR #1689, 2026-10-26); `nodePropertyExtractorLuminosityEmissionLineAGN` is
  the model to copy.
- **Heap-allocated modifier lists.** `multiModifierList` nodes must be heap-allocated.
- **`ALLOCATE(SOURCE=)`** does an intrinsic copy, bypassing reference counting.
- **Never assign a constructor result to a *pointer* target of a type with finalizable components.**
  `allocate(p); p=rootFinder(...)` (also via a local pointer) compiles to a copy to address 0 (`movaps %xmm0,0x0`)
  and segfaults the first time it runs, at -O3 and -O0 (gfortran 16.0.1 20260322; 2026-10-05, fixed this way in
  PR #1569). A ~60-line standalone program reproduced it: a type with a pointer counter component that has a
  defined assignment and a final procedure, held twice inside another type. That reproducer was not kept, and it
  has not been reported upstream. An *allocatable* scalar is fine. Use `allocate(p,source=rootFinder(...))`, which
  is safe only while the new object holds no managed resources. A path that hasn't run yet can't show the bug, so
  scan the binary: `objdump -d` the function and grep `mov[a-z]* +[^,]+,0x[0-9a-f]+$` (hits mean miscompiled; check
  the grep against a known-bad build first).
- **A class constructor hands back an object with reference count *zero*.** `x => someClass(parameters,i)`
  therefore leaves the member unowned: the generated `assignment(=)` of the enclosing class increments it to 1, and
  finalizing the constructor's own function result drops it back to 0 and frees the member, which is then used
  dangling. Build members through `<objectBuilder class="..." name="list%member_" source="parameters" copy="i"/>`,
  which increments, and release them with `<objectDestructor>`, never a bare `deallocate`. `nodeOperatorMulti` is the
  model to copy. Root cause of the `outputTimes=union` and `mergerTreeBuildMasses=union` segfaults (2026-09-21);
  measured by printing `loc()` and `%referenceCount` at construction and at use.
- **Ordered comparisons with NaN crash `Galacticus.exe`, not just unit tests.** The main build also has
  `-ffpe-trap=invalid,zero,overflow` (`Makefile`), and `x < y`, `max(x, y)` or a `mask=a < b` with a quiet NaN
  operand signal "invalid", so a NaN used as a "not set" marker raises SIGFPE when compared. Use a `huge(0.0d0)`
  sentinel or a logical flag, or test with `ieee_is_nan` first (`x == y` and `x /= y` do not trap). Writing a NaN to
  output is fine. Seen 2026-10-28 in the then-new `branchingVariance` operator: single-threaded it was a clear SIGFPE,
  but multi-threaded it surfaced as an HDF5 "not a datatype" error in another operator.
- **HDF5 objects and `hdf5Access`.** Declare `hdf5File`/`hdf5Group` variables inside a `block` that ends *before*
  `!$ call hdf5Access%unset()`. Declared at procedure scope, they are finalized at procedure exit, after the lock is
  dropped, racing another thread that holds the lock while using HDF5 (from a review of the dust emission
  templates, 2026-09-15; pattern in `source/dust/emission/Draine_Li2007.F90`).

**Why:** these produce heap corruption or deadlocks far from their cause.

**How to apply:** check any change that copies, assigns or privatizes objects against this list.

## Smaller code traps
*Added 2026-08-07; verified 2026-09-23*

- **A `varying_string` expression in an argument list leaks.** gfortran does not finalize the temporary, so
  `File_Exists(path//'_x'//suffix)`, or any such concatenation or function result written directly as an actual
  argument, leaks on every call. Build it in an assignment and pass the variable. The repo documents this on
  `shellEscape` in `source/system/command.F90`. The *Memory-Leaks-MPI* CI check catches it, but reports the
  **innermost frame in `Galacticus.exe`** of the preprocessed file, and with `-O3 -flto` that can be a *neighboring
  line holding an identical expression* - on 2026-09-23 it named `CAMB.p.F90:542` for a leak whose real site was
  line 526. Reproduce locally in ~15 min rather than waiting ~1 hour for CI: run a small model under
  `valgrind --leak-check=full --xml=yes --xml-file=...`, then count `Leak_DefinitelyLost` records by their innermost
  `Galacticus.exe` frame (the logic is in `testSuite/test-memory-leaks.py`). A leak in `git2.c` line 82 is
  pre-existing and unrelated.
- **`varying_string` arrays are not allocated on assignment.** `ISO_Varying_String` gives them a defined
  (elemental) assignment, and allocate-on-assignment applies only to *intrinsic* assignment, so assigning a
  function's array result to an unallocated array segfaults inside `op_assign_vs_vs` ("Location address is not set"
  in gdb). Allocate it first. Scalars are fine.
- **`File_Unlock(lock)` creates files.** With the default `sync=.true.` it opens the *target* file with
  `status='unknown'`, creating it if absent. That breaks double-checked locking, so pass `sync=.false.` when
  releasing a shared lock to upgrade it.
- **Memoization keys.** Several classes memoize one previous result keyed on a few scalars (e.g. Diemer–Kravtsov on
  mass and time). Before passing a new argument (e.g. `node`) into such a routine, grep for `Previous` state and
  extend the key.
- **Type-opener regex.** It's duplicated across 8 files (`Fortran/Utils.py`, `DeepCopyActions.py`,
  `StateStorable.py`, `TypeScan.py`, `FunctionClass/Utils.py`, `Hierarchy.py`, `moduleDependencies.py`,
  `staticAnalyzer.py`), and most copies miss `bind(c)` (issue #1381). If you change type parsing, fix all eight
  together.
- **A dummy procedure in a functionClass method has no interface (gfortran-16).** Implementations compile as
  `module procedure` bodies in submodules, and in that form a dummy declared `procedure(<abstract interface>)` loses
  its explicit interface (confirmed 2026-09-21 with a 20-line reproducer; the full `module subroutine f(...)` form is
  fine). Keyword or optional arguments are rejected, and positional calls are *unchecked*: passing a real for a
  logical compiles silently. `computationalDomainVolumeIntegrator`'s `integrand` rides this path and happens to
  work. Put such a helper in a plain module (e.g. `Tasks_Build_Tool_Utilities`) instead, where the interface is
  enforced.
- **Never filter `^ *!` lines when reading or comparing Fortran.** OpenMP directives (`!$omp critical`,
  `!$omp threadprivate`) are comments syntactically, so a `grep -v '^ *!'` view or a comment-stripping diff hides
  them. On 2026-09-21 a template copied from such a view dropped an `!$omp critical` guarding a shared `save`
  variable, and the comment-stripping comparison meant to verify it passed; only the `01-claude-review` hook caught
  it. Keep `!$` lines in every normalization.

**Why:** each of these passed the test suite while being wrong.

**How to apply:** check this list when touching locking, memoized routines or the type parser.

## Reporting a numerical failure the run survives: `Warn` is invisible, and a GSL status is not just a flag
*Added 2026-10-26; verified 2026-10-26*

Moved to the developer guide: *Tolerating a numerical failure: which statuses leave a usable result* and *Warnings
that must be seen in a successful run*, in `docs/manuals/developer-guide/error-handling.rst`. Read them before adding a
user-facing warning, or a parameter that decides whether a numerical failure is fatal.
