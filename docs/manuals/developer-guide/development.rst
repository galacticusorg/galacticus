.. _manual-sec-development:

Developing Galacticus
=====================

The Galacticus Build System
---------------------------

Galacticus uses a GNU Make-based build system. Dependencies are automatically discovered, and extensive meta-programming and preprocessing of source files is undertaken as part of the build. As a result, the build system is complicated. This section describes it in detail: the Makefile and its configuration knobs, the incremental-build machinery, the special comment markers that the build scripts read from source files, the automatic-discovery and code-generation steps, and the files everything produces.

A build proceeds, conceptually, in four phases:

#. **Cataloging**: the source tree is scanned for directives, module declarations, ``use`` statements, ``include`` statements, and ``program`` units, producing the generated makefiles and XML catalogs under ``$(BUILDPATH)`` that drive the rest of the build (see :galacticus-ref:`buildDiscovery`).
#. **Preprocessing**: each Fortran source file is run through the source-tree preprocessor (see Section :galacticus-ref:`sourceTreePreprocessor`), which expands directives into generated code, producing ``*.p.F90`` files.
#. **Compilation**: preprocessed sources are compiled, with compiler diagnostics remapped back to original source line numbers.
#. **Linking**: for each executable, the set of required objects is read from its ``.d`` file, input-parameter and source-digest metadata objects are generated and compiled, and everything is linked with the required external libraries in the correct order.

The Makefile
~~~~~~~~~~~~

Files generated during the build (with the exception of final executables) are written to the directory specified by the ``BUILDPATH`` variable. This is normally ``work/build/`` (and is set by the ``Makefile``) unless a special build configuration is requested via ``GALACTICUS_BUILD_OPTION``:

``default``
   Standard build (``BUILDPATH=./work/build``).

``MPI``
   MPI-enabled build (``BUILDPATH=./work/buildMPI``); compilers default to ``mpif90``/``mpicc``/``mpic++`` (overridable via ``MPIFCCOMPILER`` etc.) and ``-DUSEMPI`` is defined.

``lib``
   Shared-library build of ``libgalacticus.so`` (``BUILDPATH=./work/buildLib``); adds ``-fPIC`` and heap trampolines, and activates the library-interface generators (see Section :galacticus-ref:`buildLibraryInterface`).

``gprof``, ``perf``, ``odeprof``
   Profiling variants (``-pg``, ``-fno-omit-frame-pointer``, and ``-DPROFILE`` respectively), each with its own build directory and executable ``SUFFIX``.

``compileprof``
   Build-profiling mode: every recipe is run through the ``profiler.sh`` wrapper SHELL to collect timing and memory data (see Section :galacticus-ref:`buildProfiling`).

.. _manual-sec-buildKnobs:

Configuration knobs
^^^^^^^^^^^^^^^^^^^

The following make variables (settable on the command line or from the environment) control the build:

``FCCOMPILER``, ``CCOMPILER``, ``CPPCOMPILER``
   The Fortran, C, and C++ compilers (defaults: ``gfortran``, ``gcc``, ``g++``; MPI builds default to the MPI wrappers, overridable via ``MPIFCCOMPILER`` etc.).

``GALACTICUS_FCFLAGS``, ``GALACTICUS_CFLAGS``, ``GALACTICUS_CPPFLAGS``, ``GALACTICUS_F77FLAGS``
   Extra per-language flags appended to the built-in flag sets — typically ``-I``/``-L`` paths for locally-installed dependencies. Note that the availability probes for ANN, qhull, and libmatheval compile against *only* ``GALACTICUS_CPPFLAGS`` (see :galacticus-ref:`buildConfigProbes`).

``LTO``
   ``enabled`` (default) or ``disabled``. Link-time optimization uses ``-flto=jobserver`` so that the parallelism of the link-time recompilation is governed by make's job server rather than the machine's CPU count; disable on platforms where LTO causes problems (e.g. ``dsymutil`` memory exhaustion on Apple Silicon).

``LOCKMD5``
   ``yes`` (default) or ``no``. Whether the source-digest computation takes ``flock`` locks on its per-file ``.md5`` sidecar caches. Locks are only actually used under parallel make (``-j``); they serialize concurrent builds sharing a build directory (e.g. several executables linking at once, each running ``sourceDigests.py``). The Makefile forces ``no`` when only a single target was named on the command line.

``OFDLOCKS``
   ``enabled`` (default) or ``disabled``: whether to use open-file-description locks (``-DOFDLOCKS`` vs. ``-DNOOFDLOCKS``) for run-time file locking.

``USEGIT2``
   Set to ``no`` to force the build to *not* use ``libgit2`` (falling back to the ``git`` command line) even when it is available. Static builds force this automatically.

``GALACTICUS_OBJECTS_DEBUG``
   Set to ``yes`` to define ``-DOBJECTDEBUG``, enabling object life-cycle debugging output.

``SUFFIX``
   Appended to executable names (set automatically by the profiling build options, e.g. ``Galacticus.exe_gprof``).

``PREPROCESSOR``
   Name of the C preprocessor, recorded in the build metadata (the ``.p.Inc``\ →\ ``.inc`` rule currently invokes ``cpp`` directly).

.. _manual-sec-buildConfigProbes:

Feature-availability probes
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Optional dependencies are detected by compile probes whose results are written to ``$(BUILDPATH)/Makefile_Config_<Name>`` fragments (included into the Makefile) that append ``-D<SYMBOL>AVAIL`` or ``-D<SYMBOL>UNAVAIL`` to the relevant flags variables. The identically-shaped probes (OFD file locking, FFTW3, ANN, qhull, libmatheval) share the ``CONFIG_PROBE`` canned recipe in the Makefile — adding a new probe is a ``-include`` line, a rule, and one ``$(call CONFIG_PROBE,...)``. The Proc probe (which must also *run* its compiled binary) and the libgit2 probe (which has static-build and opt-out variants) are hand-written.

Two caveats a maintainer should know:

* The C/C++ probes for ANN, qhull, and libmatheval compile against only ``GALACTICUS_CPPFLAGS``, *not* the full ``CPPFLAGS``. The latter adds ``-I$(BUILDPATH)/``, which can contain zero-byte header stubs left by an earlier unavailable-result build (created by the Makefile's generic ``%.h`` rule for headers that were preprocessed out). Such a stub satisfies the probe's ``#include`` and flips the probe to a spurious *available* result on a machine without the real headers — after which the real compile fails. The libgit2 probe hit this trap first and carries the original explanation.
* Probe results are cached in the ``Makefile_Config_*`` fragments and only regenerated when the probe *source* changes — not when the machine's installed libraries change. After installing or removing an optional dependency, delete the corresponding ``$(BUILDPATH)/Makefile_Config_*`` file (or the whole build directory) to re-probe.

Executable linking
^^^^^^^^^^^^^^^^^^

There is no generic ``%.exe`` pattern rule: ``findExecutables.py`` discovers every program in the source tree and writes an explicit rule for each into ``$(BUILDPATH)/Makefile_All_Execs``. Those generated rules contain only a dependency line and a ``$(call LINK_EXECUTABLE,...)``; the link procedure itself lives in two canned recipes in the main Makefile, so it has a single owner:

``LINK_METADATA``
   Enumerates every input parameter the target can accept (``parameterDependencies.py``, compiled into ``<stem>.parameters.o``) and computes per-type source digests (``sourceDigests.py``, compiled into ``<stem>.md5s.o``). Shared by executables and ``libgalacticus.so``.

``LINK_EXECUTABLE``
   ``LINK_METADATA`` plus the link step, which reads the object list from the target's ``.d`` file, appends the library flags computed by ``libraryDependencies.py``, and captures the linker's exit status so that a failed link fails the recipe (see :galacticus-ref:`buildStatusPropagation`).

.. _manual-sec-buildIncremental:

Incremental-Build Machinery
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because so much of the build is generated, naive use of file timestamps would cause huge rebuild cascades: regenerating a file whose content did not change would still advance its mtime and invalidate everything downstream. The build system uses three cooperating mechanisms to keep incremental builds minimal. All three follow the same principle — *a generated file's mtime advances only when its content actually changes* — so understanding them is essential before modifying any generator.

.. _manual-sec-buildUpdateFiles:

Update-sentinel (".up") files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a generated file ``X`` whose generator should re-run when its inputs change, but whose *consumers* should re-run only when ``X``'s content changes, the build uses a sentinel pair::

   $(BUILDPATH)/X.up : <inputs>
           <generator writes X only-if-changed, then touches X.up>
   $(BUILDPATH)/X : $(BUILDPATH)/X.up
           @true

Make compares the *inputs* against the sentinel ``X.up`` (touched on every generator run), while downstream targets compare against ``X`` itself (mtime preserved when content is unchanged). The only-if-changed write plus sentinel touch is implemented by ``Galacticus.Build.FileChanges.update(..., prove_update=True)``: the generator writes to a temporary file, which either replaces ``X`` (if different) or is discarded (if identical). This pattern is used for preprocessed sources (``*.p.F90``), generated include files (``*.Inc``), and the OpenMP critical-section catalog, among others.

.. _manual-sec-buildBlobCaches:

Per-file cache blobs
^^^^^^^^^^^^^^^^^^^^

The cataloging scripts each keep a pickle cache ("blob") under ``$(BUILDPATH)`` recording per-source-file scan results, so a re-run rescans only changed files. The shared protocol (implemented in ``Galacticus.Build.ScanCache`` and documented in its docstring) is:

* Each file's entry is keyed by a canonical identifier of its path, and records the list of files whose content contributed to the entry (the file itself plus any followed includes or directive-referenced files).
* The blob's own mtime is the freshness reference: any tracked file with mtime newer than the blob forces a rescan of that entry. Adding or removing source files forces a full rescan.
* A missing, corrupt, or wrong-format blob silently degrades to a full rescan — cache problems can never fail the build. Scripts whose recorded entry structure changes bump an embedded version stamp so old-format caches are discarded.
* Blobs are themselves rewritten only-if-changed (preserving the blob mtime when nothing changed), which keeps the freshness window identical to a serial run. ``codeDirectivesParse.py``'s worker results are canonicalized before pickling for exactly this reason — so a parallel scan pickles byte-identically to a serial one.

Module-file update avoidance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compiling a Fortran source regenerates the ``.mod`` files of every module it defines, which would invalidate every user of those modules even when the module *interfaces* did not change. The ``%.o`` rule therefore compiles with ``-J$(BUILDPATH)/moduleBuild/`` (a staging directory) and then, for each module listed in the source's ``.m`` file, compares the staged ``.mod`` against the current one — moving it into place only if it differs. The generated ``<module>.mod: <object>`` rules in ``Makefile_Module_Dependencies`` likewise rebuild the owning object only when the ``.mod`` file is actually missing. Together these stop recompilation cascades through the module graph.

The directive-catalog stamp
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The three directive-catalog scripts (``codeDirectivesParse.py``, ``stateStorables.py``, ``deepCopyActions.py``) run in a single recipe, in that order, because the latter two read the ``directiveLocations.xml`` written by the first. Their outputs (``Makefile_Directives``, ``directiveLocations.xml``, ``stateStorables.xml``, and ``deepCopyActions.xml``) are all products of that one recipe: a stamp file, ``$(BUILDPATH)/directiveCatalogs.stamp``, records when the scripts last ran, and each output that other rules name as a prerequisite has an empty-recipe rule on the stamp. Since the scripts write their outputs only-if-changed, downstream targets rebuild only when catalog *content* changes, while the stamp's mtime stops the scripts re-running when nothing is out of date.

.. _manual-sec-buildMarkers:

Source-Code Markers
~~~~~~~~~~~~~~~~~~~

The build scripts key on a number of special comment markers in the source. None of these affect compilation directly — they exist purely for the build system — and a maintainer grepping for "why does this file depend on that?" usually ends up at one of these:

``!![`` … ``!!]``
   Delimits an embedded XML directive block (a *code directive* — see :galacticus-ref:`buildDiscoveryDirectives`). Parsed by ``Galacticus.Build.Directives`` and by the source-tree preprocessor; skipped by the plain-text scanners.

``!!{`` … ``!!}``
   Delimits an embedded documentation block (reStructuredText/LaTeX). Skipped by every scanner (prose inside could otherwise look like code, e.g. a ``program`` statement).

``!/ exclude``
   In a file containing a ``program`` unit: build the executable, but exclude it from the ``all`` target (``findExecutables.py``).

``!; <library>``
   The containing object requires linking against ``<library>``. Recorded in the object's ``.fl`` file and consumed by ``libraryDependencies.py`` at link time (``useDependencies.py``).

``!: <path> …`` (Fortran) / ``//: <path> …`` (C/C++)
   Explicit object-file dependencies: the named objects (usually ``$(BUILDPATH)``-relative ``.o`` paths) are added to the containing object's dependencies — used where a dependency exists through ``bind(c)`` interfaces or other channels invisible to the ``use``-statement scan (``useDependencies.py``).

``! NO_USES``
   Appended to an ``include`` statement of a *generated* include file: opts that include out of the automatically-generated dependency chain that orders ``Makefile_Use_Dependencies`` after include generation (``includeDependencies.py``). Used for includes (such as build-metadata files) whose content must not retrigger dependency analysis.

``!$GLC …``
   Directives consumed by ``postprocess.py`` to suppress spurious compiler warnings: ``!$GLC function attributes unused ::`` (unused module procedures), ``!$GLC attributes unused ::`` (unused variables, per enclosing unit), ``!$GLC attributes initialized ::``, ``!$GLC attributes interoperable ::``, ``!$GLC ignore outlive ::``, and ``!$GLC ignore unused ::``.

``!-->``
   Instrumentation lines written by the source-introspection machinery (including the ``.lmap`` line maps); ignored by directive extraction.

Two ordinary Fortran constructs also carry hard-wired dependency implications in ``useDependencies.py``: a named ``!$omp critical(<name>)`` section implies a dependency on ``openmp_utilities_data.mod`` (the section enumeration — see :galacticus-ref:`buildOpenMPCritical`), a ``!$omp parallel`` construct implies ``events_filters.mod``, and a ``program`` unit implies ``iso_varying_string.mod``.

.. _manual-sec-buildDiscovery:

Automatic Discovery
~~~~~~~~~~~~~~~~~~~

Interdependencies between source files, together with requirements for auto-generated code, are automatically discovered during the build process by processing of source files. All such automatic discovery is described below. The scanning scripts share several conventions: they parallelize their per-file scans through ``Galacticus.Build.ParallelScan`` (a fork-based pool that preserves result order, so outputs are byte-identical to a serial run, and reads make's ``-j`` setting from ``MAKEFLAGS``); they cache per-file results in blobs (see :galacticus-ref:`buildBlobCaches`); and they emit their generated makefiles and catalogs only-if-changed.

.. _manual-sec-buildDiscoveryDirectives:

Code Directive Parsing
^^^^^^^^^^^^^^^^^^^^^^

The modular nature of Galacticus is enabled through the use of directives embedded within the source code (XML documents inside ``!![`` … ``!!]`` comment blocks) which instruct the build system how to link together the various functions. Parsing of these directives is handled by the ``codeDirectivesParse.py`` script.

The script also writes ``$(BUILDPATH)/Makefile_Directives``, which carries the extra build dependencies implied by certain directives: for each ``functionClass`` and for the ``componentBuilder`` directive it makes the preprocessed source depend on the other files that feed it (e.g. ``_class.p.F90`` depends on every ``<component>`` directive file and on the component generator's Python sources), and it orders ``Makefile_Use_Dependencies`` and ``Makefile_Component_Includes`` after the preprocessed ``_class`` file.

The script also outputs ``$(BUILDPATH)/directiveLocations.xml``, which lists the files containing each directive. This is used by the other cataloging scripts to permit rapid processing of files associated with each directive.

State storables and deep-copy actions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two further catalogs are generated in the same recipe (see the directive-catalog stamp above), and are read by the source-tree preprocessor when generating state-store and deep-copy code — which is why every preprocessing rule lists them as prerequisites:

* ``stateStorables.xml`` (from ``stateStorables.py``): every ``functionClass``, its concrete implementations, ``stateStorable`` types, and static event hooks.
* ``deepCopyActions.xml`` (from ``deepCopyActions.py``): every derived type participating in a ``deepCopyActions`` directive, expanded over each type's inheritance chain.

.. _manual-sec-buildOpenMPCritical:

OpenMP critical sections
^^^^^^^^^^^^^^^^^^^^^^^^

``enumerateOpenMPCriticalSections.py`` scans every Fortran source for named ``!$omp critical(<name>)`` sections and assigns each name an ID (in sorted-name order — a determinism contract the instrumented code relies on), writing ``$(BUILDPATH)/openMPCriticalSections.xml`` plus the ``openMPCriticalSections.count.inc`` and ``openMPCriticalSections.enumerate.inc`` include fragments used by the OpenMP wait-time instrumentation. The rule depends on all Fortran sources through an ``.up`` sentinel, so adding, removing, or renaming a critical section regenerates the enumeration, while a no-op regeneration does not cascade into re-preprocessing (the ``%.p.F90.up`` rules depend on the ``.xml``).

.. _manual-sec-buildExecutables:

Executable Files
^^^^^^^^^^^^^^^^

Files which produce executables (including the ``Galacticus.exe`` executable) are discovered by the ``scripts/build/findExecutables.py`` script: all Fortran files in the ``source`` directory (excluding vendored code under ``source/external/``) are parsed, and executable files are identified as those containing a ``program`` statement. For each, a rule invoking the Makefile's ``LINK_EXECUTABLE`` canned recipe is written to ``$(BUILDPATH)/Makefile_All_Execs``, and the executable is added to the ``all`` target unless the file carries a ``!/ exclude`` marker. Build artifacts mirror the hierarchical source tree (e.g. ``$(BUILDPATH)/tests/nodes.o``) while executables keep the historical flat, dot-separated names (``tests.nodes.exe``) so existing invocations keep working.

.. _manual-sec-buildModulesProvided:

Modules Provided
^^^^^^^^^^^^^^^^

Determination of which files provide which Fortran modules is carried out by the ``scripts/build/moduleDependencies.py`` script, which generates rules describing these dependencies (and how to build the module file by compiling the corresponding source file) and writes them to ``$(BUILDPATH)/Makefile_Module_Dependencies``.

Additionally, rules are generated which describe how to build the following classes of file:

``*.mod.d``
   Contains a list of all object files upon which the module depends. Used in constructing the final set of objects which must be linked to build an executable.

``*.mod.gv``
   Contains a list of all source files upon which the module depends. Used in constructing :term:`GraphViz` visualizations of dependencies.

``*.m``
   Contains a list of all modules provided by the corresponding object file. Used when building object files to test whether corresponding module files have been changed — allows avoidance of recompilation cascades if modules do not need to be updated (see :galacticus-ref:`buildIncremental`).

functionClass submodules
""""""""""""""""""""""""

Each concrete implementation of a ``functionClass`` is built *on the fly* as a Fortran submodule of the module hosting the class. ``moduleDependencies.py`` discovers these by walking each implementation file listed in ``directiveLocations.xml`` for the class's directive, synthesizing one submodule record per derived type (named ``<implementation>_``, with parent submodules following the type's ``extends`` chain). The generated rules encode two special behaviors: an implementation's preprocessed source (``<impl>_.p.F90``) is produced *as a side effect of preprocessing the functionClass base file*, so the rules gate each implementation behind its base (and force-regenerate the base if an implementation's file is missing); and ``<module>@<submodule>.smod`` targets tie submodule files to their owning objects.

.. _manual-sec-buildModulesUsed:

Modules Used
^^^^^^^^^^^^

Determination of which files use which Fortran modules is carried out by the ``scripts/build/useDependencies.py`` script, which generates rules describing these dependencies and writes them to ``$(BUILDPATH)/Makefile_Use_Dependencies``. Beyond literal ``use`` statements, the scan accounts for the implicit dependencies implied by directives (functionClass method modules, event hooks, global functions, enumerations, …) and by the source markers described in :galacticus-ref:`buildMarkers`.

Additionally, rules are generated which describe how to build the following classes of file:

``*.d``
   Contains a list of all object files upon which each object file depends. Used in constructing the final set of objects which must be linked to build an executable.

``*.gv``
   Contains a list of all source files upon which each file depends by virtue of ``use`` statements. Used in constructing :term:`GraphViz` visualizations of dependencies.

``*.fl``
   Contains a list of external libraries (one name per line) upon which each object file depends — populated from the module-to-library and include-to-library tables in ``Galacticus.Build.Libraries`` plus explicit ``!; <library>`` markers. Consumed by ``libraryDependencies.py`` to determine the final link line.

Conditional compilation
"""""""""""""""""""""""

The scanner honors preprocessor conditionals so that, for example, MPI-only ``use`` statements do not create dependencies in non-MPI builds. The set of defined macros is assembled from ``-D`` flags on active ``FCFLAGS`` lines of the Makefile and the ``Makefile_Config_*`` fragments (approximating ``ifdef``/``ifndef`` against the environment), from ``GALACTICUS_FCFLAGS``, and from the C compiler's predefined macros. Conditions that cannot be evaluated — ``ifeq``/``ifneq`` over make variables on the Makefile side, and ``#if <expression>`` or ``#elif`` chains on the source side — are treated as *active in every branch*: the scan deliberately over-approximates, since an unnecessary dependency merely costs a rebuild while a missed one silently produces stale builds. This behavior is pinned by the tests in ``scripts/build/tests/test_use_dependencies_conditionals.py``.

.. _manual-sec-buildIncludeDeps:

Included Files
^^^^^^^^^^^^^^

Dependencies on files due to the use of "``include``" statements (in both Fortran and C source files) are discovered by the ``scripts/build/includeDependencies.py`` script, which generates rules describing these dependencies and writes them to ``$(BUILDPATH)/Makefile_Include_Dependencies``.

All Fortran and C/C++ source files in the ``source`` directory are parsed. Files specified in any include statement in a source file are added as a dependency of that source file unless the source file is C/C++ and the named include file can be found in either ``source/``, a standard system include path, or in any include path specified in the ``GALACTICUS_CFLAGS`` or ``GALACTICUS_CPPFLAGS`` environment variables. This script also emits the ordering rule that makes ``Makefile_Use_Dependencies`` depend on every generated include file — the linchpin that guarantees generated includes exist before the use-scan reads them (opt out per include with ``! NO_USES``).

.. _manual-sec-parameterDependencies:

Parameter Dependencies
^^^^^^^^^^^^^^^^^^^^^^

All runtime parameters upon which a given executable may depend are discovered by the ``scripts/build/parameterDependencies.py`` script, which runs at link time (from ``LINK_METADATA``). It reads the target's ``.d`` file to find the participating objects, scans each corresponding source (preferring the preprocessed ``.p.F90``, and including the per-functionClass ``<stem>.p`` parameter listings dropped by the preprocessor) for ``inputParameter`` and ``objectBuilder`` directives, and writes ``$(BUILDPATH)/<stem>.parameters.F90`` — a ``knownParameterNames`` subroutine enumerating the de-duplicated, sorted parameter names, which the executable uses to detect unrecognized input parameters. Per-file results are cached in ``<stem>.blob``, and entries for files that leave the target are pruned.

.. _manual-sec-buildSourceDigests:

Source Digests
^^^^^^^^^^^^^^

To allow saved run states to be validated against the code that produced them, every ``functionClass``/``sourceDigest`` type participating in an executable carries an MD5 digest of its source. ``scripts/build/sourceDigests.py`` runs at link time (from ``LINK_METADATA``): for each type it hashes the defining source file plus every file it ``include``\ s (skipping non-deterministic inputs such as the build-environment and version includes), then composes each type's digest with those of its declared dependencies and inheritance chain, and emits ``$(BUILDPATH)/<stem>.md5s.c`` containing one ``char <type>MD5[]`` definition per type. Per-file digests are cached both in ``<stem>.md5.blob`` and in per-source ``.md5``/``.md5c`` sidecar files under ``$(BUILDPATH)``; sidecar access is serialized with ``flock`` locks when ``LOCKMD5`` is active (see :galacticus-ref:`buildKnobs`).

Code Generation
~~~~~~~~~~~~~~~

During build, code is generated automatically. Much of this is to link together functions within Galacticus (thereby permitting the modular nature of Galacticus), but also includes generation of various utility functions. All code generation is described below.

.. _manual-sec-codeBuildIncludeDirectives:

The node-component class hierarchy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``nodeComponent`` class hierarchy is generated from the ``<component>`` directives spread across ``source/objects/nodes/components/`` by the ``componentBuilder`` source-tree process hook (``python/Galacticus/Build/SourceTree/Process/ComponentBuilder.py``). A single ``<componentBuilder/>`` directive in ``source/objects/nodes/_class.F90`` triggers the hook while that file is preprocessed (see Section :galacticus-ref:`sourceTreePreprocessor`): it gathers, validates, and parses every ``<component>`` directive through the generator package ``python/Galacticus/Build/Components``, runs the generator's phased pipeline to synthesize the Fortran source, processes that fragment through the full preprocessor pipeline (expanding any directives the generator embedded in its own output), and grafts the result in place of the directive. The generated code therefore lands directly in ``$(BUILDPATH)/objects/nodes/_class.p.F90`` rather than in a separately-generated include file.

.. _manual-sec-codeBuildPrePostProcess:

Pre-processing and Post-processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before being compiled, all Fortran source files are passed through a tree-based preprocessor. The functionality of this preprocessor is described in Section :galacticus-ref:`sourceTreePreprocessor`. The preprocessor is invoked on each source file by the ``preprocess.py`` script, which takes a ``source/*.F90`` file as input and outputs a preprocessed file ``$(BUILDPATH)/*.p.F90`` (written only-if-changed — see :galacticus-ref:`buildUpdateFiles`), together with a ``*.p.F90.lmap`` line-map sidecar. Each ``.lmap`` record has the form::

   !--> <original-line> <preprocessed-line> "<original-source-file>"

marking the start of a run of lines in the preprocessed file that maps back to the given original location. The full set of source-tree process hooks is registered by importing ``Galacticus.Build.SourceTree.Process.all`` — every driver that runs the preprocessor imports this one module, so a new Process submodule is added in exactly one place.

Output from the compiler is captured to a diagnostics file and fed through the ``postprocess.py`` script, which:

* uses the ``.lmap`` to translate line numbers in errors and warnings from the preprocessed source back to the original source;
* filters out known-spurious warnings (e.g. GCC `PR58175 <https://gcc.gnu.org/bugzilla/show_bug.cgi?id=58175>`_ scalar-finalizer warnings) and warnings suppressed by ``!$GLC`` markers (see :galacticus-ref:`buildMarkers`);
* exits immediately when the compiler produced no diagnostics (the common case).

.. _manual-sec-buildStatusPropagation:

Diagnostics are written to a file and postprocessed *after* capturing the compiler's exit status — not piped — because in a pipeline the recipe's status would be the postprocessor's, and a hard compiler failure that emits no recognizable error line (an internal compiler error, a segfault, an out-of-memory kill) would be reported as success. The recipe fails if either the compiler or the postprocessor reports failure. Link steps are handled identically, with ``postprocessLinker.py`` filtering irrelevant static-link warnings and flagging ``undefined reference`` / ``ld returned N exit status`` diagnostics as a second line of defense. The ``profiler.sh`` wrapper SHELL used by ``compileprof`` builds likewise propagates the wrapped command's exit status.

.. _manual-sec-buildLibraryInterface:

The Shared-Library Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Library builds (``GALACTICUS_BUILD_OPTION=lib``) generate C-interoperable Fortran wrappers, and a Python ``ctypes`` interface, for every ``functionClass`` registered in ``source/libraryClasses.xml``:

``scripts/build/libraryInterfaces.py``
   The generator. Parses each implementation file once, classifies every constructor and method argument/return type against the rules in ``python/LibraryInterfaces/Classification.py`` (skipping, with a warning, anything the pipeline cannot translate), and emits: ``$(BUILDPATH)/libgalacticus.Inc`` (the library's initialization core), per-class wrapper units ``$(BUILDPATH)/libgalacticus/<class>.F90`` plus per-implementation constructor wrappers ``<class>__<impl>.F90`` (split so gfortran gets small compilation units), and the Python interface ``galacticus.py``. Wrapper units are written only-if-changed, so catalog changes recompile only genuinely affected wrappers.

``scripts/build/libraryInterfacesDependencies.py``
   Enumerates the generated wrapper units and writes ``$(BUILDPATH)/Makefile_Library_Dependencies`` with their preprocess/compile rules (mirroring the main ``%.o`` recipe) and the aggregate ``libgalacticus_classes.d``.

``scripts/build/libraryInterfacesAudit.py``
   A *manual* coverage-planning tool (invoked by nothing in the Makefile or CI): reports, for every ``functionClass`` in the source tree, whether it is ready to be registered in ``libraryClasses.xml``, blocked by a missing dependency, or blocked by an unsupported argument/return type. Because it imports the same ``Classification`` rules as the generator, its verdicts track the generator by construction.

Because the generators are slow, they are activated only for library builds.

Parameter Catalog and Schema
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Separate from the compile path, the build can generate machine-readable descriptions of every ``functionClass`` implementation's input parameters:

``make parameters-catalog``
   Runs ``scripts/build/parameterCatalog.py`` to produce ``$(BUILDPATH)/parameters.catalog.json`` — names, inferred types, defaults, and nesting for every parameter. Consumed by parameter-file validation tooling; CI checks its provenance coverage with ``parameterCatalogCheck.py``.

``make parameters-schema``
   Runs ``scripts/build/parameterSchema.py`` to regenerate the committed editor-assistance schema ``schema/parameters.xsd`` from the catalog; CI checks that it is up to date. ``parameterValidate.py`` validates parameter files against the catalog (used by CI and the migration tooling).

Build Files
~~~~~~~~~~~

All files involved in the build process are summarized below.

Scripts (all under ``scripts/build/``):

``buildProfiler.py``
   Generates an HTML report from a ``compileprof`` build log (see Section :galacticus-ref:`buildProfiling`);

``codeDirectivesParse.py``
   Discovers and parses directives embedded in source files. Generates ``make`` rules for the resulting dependencies (in ``Makefile_Directives``) and a file providing a mapping of which files contain each directive (``directiveLocations.xml``);

``deepCopyActions.py``
   Catalogs every derived type participating in a ``deepCopyActions`` directive (see :galacticus-ref:`buildDiscovery`);

``enumerateOpenMPCriticalSections.py``
   Enumerates named OpenMP critical sections (see Section :galacticus-ref:`buildOpenMPCritical`);

``findExecutables.py``
   Discovers (and generates rules for) files which generate executables (see Section :galacticus-ref:`buildExecutables`);

``includeDependencies.py``
   Discovers (and generates rules for) dependencies between files arising from the use of "``include``" statements (see Section :galacticus-ref:`buildIncludeDeps`);

``libraryDependencies.py``
   Determines the set of libraries that must be linked with each executable (from its ``.fl`` files and the dependency graph in ``Galacticus.Build.Libraries``), and ensures they are ordered correctly for static linking;

``libraryInterfaces.py``, ``libraryInterfacesDependencies.py``, ``libraryInterfacesAudit.py``
   Generate (and audit coverage of) the shared-library interface (see Section :galacticus-ref:`buildLibraryInterface`);

``moduleDependencies.py``
   Discovers (and generates rules for) dependencies of module files on their source file statements (see Section :galacticus-ref:`buildModulesProvided`);

``parameterCatalog.py``, ``parameterCatalogCheck.py``, ``parameterSchema.py``, ``parameterValidate.py``
   Generate and validate the typed parameter catalog and editor schema;

``parameterDependencies.py``
   Determines all parameters upon which a given executable depends (see Section :galacticus-ref:`parameterDependencies`);

``postprocess.py``, ``postprocessLinker.py``
   Postprocess compiler and linker diagnostics (see Section :galacticus-ref:`codeBuildPrePostProcess`);

``preprocess.py``
   Preprocesses Fortran source code to provide various extended functionality (see Section :galacticus-ref:`codeBuildPrePostProcess`);

``profiler.sh``
   SHELL wrapper used by ``compileprof`` builds to time every recipe (see Section :galacticus-ref:`buildProfiling`);

``sourceDigests.py``
   Computes per-type source MD5 digests for state-file validation (see Section :galacticus-ref:`buildSourceDigests`);

``stateStorables.py``
   Catalogs functionClasses, their implementations, and stateStorable types (see :galacticus-ref:`buildDiscovery`);

``staticRelinker.py``
   CI-only helper (macOS static builds): re-runs a logged link command with system dynamic libraries hidden so static archives are picked up;

``useDependencies.py``
   Discovers (and generates rules for) dependencies between files originating from module ``use`` statements (see Section :galacticus-ref:`buildModulesUsed`).

Generated files (all under ``$(BUILDPATH)`` unless noted):

``Makefile_All_Execs``
   Rules for building each executable — generated by ``findExecutables.py`` (see Section :galacticus-ref:`buildExecutables`);

``Makefile_Component_Includes``
   Dependencies of the ``nodeComponent`` class hierarchy module (``objects/nodes/_class``) on the per-implementation ``bound_functions`` include files; written by the component generator as ``_class.F90`` is preprocessed (see Section :galacticus-ref:`codeBuildIncludeDirectives`);

``Makefile_Config_*``
   Feature-availability probe results (see :galacticus-ref:`buildConfigProbes`);

``Makefile_Directives``
   Extra build dependencies implied by ``functionClass`` and ``componentBuilder`` directives (see Section :galacticus-ref:`buildDiscoveryDirectives`);

``Makefile_Include_Dependencies``
   Dependencies of source files on files included via "``include``" statements — generated by ``includeDependencies.py`` (see Section :galacticus-ref:`buildIncludeDeps`);

``Makefile_Library_Dependencies``
   Preprocess/compile rules for the generated library wrapper units (library builds only; see Section :galacticus-ref:`buildLibraryInterface`);

``Makefile_Module_Dependencies``
   Rules which describe how to build module files from their source files (see Section :galacticus-ref:`buildModulesProvided`);

``Makefile_Use_Dependencies``
   Dependencies between source files originating from module ``use`` statements (see Section :galacticus-ref:`buildModulesUsed`);

``directiveCatalogs.stamp``
   Timestamp of the last directive-catalog run (see :galacticus-ref:`buildIncremental`);

``directiveLocations.xml``
   A map of which files contain each code directive (see Section :galacticus-ref:`buildDiscoveryDirectives`);

``stateStorables.xml``, ``deepCopyActions.xml``
   Catalogs read by the preprocessor when generating state-store and deep-copy code (see :galacticus-ref:`buildDiscovery`);

``openMPCriticalSections.xml``, ``openMPCriticalSections.count.inc``, ``openMPCriticalSections.enumerate.inc``
   The named-critical-section enumeration (see Section :galacticus-ref:`buildOpenMPCritical`);

``*.d``
   Lists of dependencies on object files, accumulated to find the final set of files which must be linked to build each executable;

``*.m``
   Lists of modules provided by each object file, used when testing whether modules have been updated (see :galacticus-ref:`buildIncremental`);

``*.fl``
   Lists of libraries upon which each object file depends (one name per line), used to determine the final set of libraries to link;

``*.gv``
   Lists of source files upon which each file depends, used in the construction of :term:`GraphViz` representations of file dependencies;

``*.p.Inc`` / ``*.inc``
   Preprocessor-intermediate / fully-processed include files generated from ``source/*.Inc`` files (e.g. the per-implementation ``bound_functions.Inc`` sources). The intermediate is deliberately named ``.p.Inc`` rather than ``.Inc``: the latter would differ from the processed ``.inc`` only by case and the two collide on case-insensitive filesystems (macOS APFS);

``*.p.F90``
   Preprocessed Fortran source files;

``*.p.F90.lmap``
   Line maps from preprocessed back to original sources (see Section :galacticus-ref:`codeBuildPrePostProcess`);

``*.up``
   Update sentinels recording when a generator last ran (see :galacticus-ref:`buildUpdateFiles`);

``*.blob``, ``*.md5.blob``
   Per-file pickle caches used by the cataloging and metadata scripts (see :galacticus-ref:`buildBlobCaches`);

``*.md5``, ``*.md5c``
   Per-source digest sidecar caches (see Section :galacticus-ref:`buildSourceDigests`);

``<stem>.parameters.F90``
   Generated ``knownParameterNames`` subroutine for an executable (see Section :galacticus-ref:`parameterDependencies`);

``<stem>.md5s.c``
   Generated per-type source-digest definitions for an executable (see Section :galacticus-ref:`buildSourceDigests`);

``libgalacticus.Inc``, ``libgalacticus/``, ``galacticus.py`` (repository root)
   The generated shared-library interface (see Section :galacticus-ref:`buildLibraryInterface`).

.. _manual-sec-buildProfiling:

Profiling the Build
~~~~~~~~~~~~~~~~~~~

Tooling is provided to help profile the build process. This is useful to identify bottlenecks during build. Profiling can be switched on using the ``compileprof`` build option. For example:

.. code-block:: none

   make -j16 GALACTICUS_BUILD_OPTION=compileprof Galacticus.exe >& build.log

Every command run by ``Make`` will be timed, its peak memory usage measured, and the results output. In the above these outputs are collected to the ``build.log`` file. Result lines have the form::

   ++Task: {<start>|<stop>|<maxRSS>} '<command>'

with RFC 3339 start/stop times and the peak resident set size in kilobytes (:math:`-1` when memory usage could not be measured, for example because GNU ``time`` is unavailable) — this format is the contract between ``profiler.sh`` (producer) and ``buildProfiler.py`` (consumer). While you can look through these manually a script is provided that generates a report from these timing data as a web page. For example:

.. code-block:: none

   ./scripts/build/buildProfiler.py build.log profile.html --durationMinimum 10

will parse timing data from ``build.log`` and generate a web page ``profile.html`` showing the results. In this case the ``--durationMinimum 10`` option specifies that only commands which took 10 or more seconds to run should be included in the report. (Reducing this limit will lead to a very long report.)

The profile report begins by listing the total time taken for the build, together with an estimate of the peak memory used by the build---found by summing, for each second of the build, the peak memory of all included tasks running at that time. Following that is a list of all commands performed ordered by *completion* time. Next to each command is a bar which extends from the start to the end time of the command. Each second of the bar is colored to indicate the degree of build parallelism at that time---green shows maximum parallelism, while red shows no parallelism (i.e. a single command was running at that time).

Finally, each command is assigned a "cost", :math:`\chi`. This is defined as:

.. math::

   \chi = \sum_{i_\mathrm{start}}^{i_\mathrm{end}} N_i^{-1},

where :math:`i_\mathrm{start}` and :math:`i_\mathrm{end}` are the start and end times of the command, and :math:`N_i` is the number of commands running in parallel at time :math:`i`. Cost is therefore higher for commands which run longer and for commands which are executed with less parallelism. A ranked list of commands, from most to least costly is included in the report. The report also shows the peak memory used by each command, together with a ranking of commands from most to least memory-hungry, to help identify the largest memory consumers in the build.

A historical note
~~~~~~~~~~~~~~~~~

The build scripts were ported file-by-file from Perl originals in 2026. Docstrings citing ``<name>.pl:<line>`` refer to those removed Perl sources — they document which Perl idiom a piece of Python deliberately mirrors, and remain useful when investigating subtle behaviors inherited from the port. The Perl files are retrievable from git history (they were removed in commits ``5278aadf8`` and ``33dfc7939``).
