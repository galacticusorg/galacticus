.. _manual-sec-development:

Developing Galacticus
=====================

The Galacticus Build System
---------------------------

Galacticus use a GNU Make-based build system. Dependencies are automatically discovered, and extensive meta-programming and preprocessing of source files is undertaken as part of the build. As a result, the build system is complicated. In this section the build system is described and detailed.

The Makefile
~~~~~~~~~~~~

Files generated during the build (with the exception of final executables) are written to the directory specified by the ``BUILDPATH`` environment variable. This is normally ``work/build/`` (and is set by the ``Makefile``) unless a special build configuration is requested.

Automatic Discovery
~~~~~~~~~~~~~~~~~~~

Interdependencies between source files, together with requirements for auto-generated code are automatically discovered during the build process by processing of source files. All such automatic discovery is described below.

.. _manual-sec-buildDiscoveryDirectives:

Code Directive Parsing
^^^^^^^^^^^^^^^^^^^^^^

The modular nature of Galacticus is enabled through the use of directives embedded within the source code (XML documents embedded as comment lines starting with a special tag) which instruct the build system how to link together the various functions. Parsing of these directives is handled by the ``codeDirectivesParse.py`` script.

For ``include`` directives (which generate files that get included into the source), the script generates ``$(BUILDPATH)/*.xml`` files which describe the directive, and ``$(BUILDPATH)/Makefile_Directives`` which contains rules for building those include files.

The script also outputs a file ``$(BUILDPATH)/directiveLocations.xml`` which provides list of files that contain each directive. This is used by other scripts to permit rapid processing of files associated with each directive.

.. _manual-sec-buildExecutables:

Executable Files
^^^^^^^^^^^^^^^^

Files which produce executables (including the ``Galacticus.exe`` executable) are discovered by the ``scripts/build/findExecutables.py`` script, which generates rules describing how these files should be built, along with all of their dependencies, and writes them to ``$(BUILDPATH)/Makefile_All_Execs``. All Fortran files in the ``source`` directory are parsed, and executable files are identified as those containing a ``program`` statement. All executables so identified are also added to the list of objects that will be built by ``make all``.

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
   Contains a list of all modules provided by the corresponding object file. Used when building object files to test whether corresponding module files have been changed---allows avoidance of recompilation cascades if modules do not need to be updated (i.e. if the modules do not differ).

.. _manual-sec-buildModulesUsed:

Modules Used
^^^^^^^^^^^^

Determination of which files use which Fortran modules is carried out by the ``scripts/build/useDependencies.py`` script, which generates rules describing these dependencies and writes them to ``$(BUILDPATH)/Makefile_Use_Dependencies``.

Additionally, rules are generated which describe how to build the following classes of file:

``*.d``
   Contains a list of all object files upon which each object file depends. Used in constructing the final set of objects which must be linked to build an executable.

``*.gv``
   Contains a list of all source files upon which each file depends by virtue of ``use`` statements. Used in constructing :term:`GraphViz` visualizations of dependencies.

``*.fl``
   Contains a list of libraries upon which each object file depends. Used in constructing the final set of libraries which must be linked with the executable.

.. _manual-sec-buildIncludeDeps:

Included Files
^^^^^^^^^^^^^^

Dependencies on files due to the use of "``include``" statements (in both Fortran and C source files) are discovered by the ``scripts/build/includeDependencies.py`` script, which generates rules describing these dependencies and writes them to ``$(BUILDPATH)/Makefile_Include_Dependencies``.

All Fortran and C/C++ source files in the ``source`` directory are parsed. Files specified in any include statement in a source file are added as a dependency of that source file if unless the source file is C/C++ and the named include file can be found in either ``source/``, a standard system include path, or in any include path specified in the ``GALACTICUS_CFLAGS`` or ``GALACTICUS_CPPFLAGS`` environment variables.

.. _manual-sec-parameterDependencies:

Parameter Dependencies
^^^^^^^^^^^^^^^^^^^^^^

All runtime parameters upon which a given executable may depend are discovered by the ``scripts/build/parameterDependencies.py`` script, which generates XML files listing these parameters and writes them to ``$(BUILDPATH)/<executableName>.parameters.xml``. All Fortran and C++ files in the ``source`` directory are parsed, and embedded parameter definitions (i.e. embedded XML blocks in lines beginning "!\@" for Fortran or "//\@" for C++ with root element ``inputParameter``) are identified. Any files included into a source file are also searched. In the case of Fortran source files, ``source/*.F90``, the corresponding preprocessed file, ``$(BUILDPATH)/*.p.F90``, will be searched for parameter dependencies.

Code Generation
~~~~~~~~~~~~~~~

During build, code is generated automatically. Much of this is to link together functions within Galacticus (thereby permitting the modular nature of Galacticus), but also includes generation of various utility functions. All code generation is described below.

.. _manual-sec-codeBuildIncludeDirectives:

Directives
^^^^^^^^^^

Generation of code to implement the functionality of ``include`` directives is carried out by the ``buildCode.py`` script. This script simply reads an XML file generated by the ``codeDirectivesParse.py`` (see Section :galacticus-ref:`buildDiscoveryDirectives`), calls the appropriate function to generate the necessary code, passes this through the preprocessor (see Section :galacticus-ref:`sourceTreePreprocessor`), and outputs the result to the appropriate include file.

.. _manual-sec-codeBuildPrePostProcess:

Pre-processing and Post-processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before being compiled, all Fortran source files are passed through a tree-based preprocessor. The functionality of this preprocessor is described in Section :galacticus-ref:`sourceTreePreprocessor`. The preprocessor is invoked on each source file by the ``preprocess.py`` script, which takes a ``source/*.F90`` file as input and outputs a preprocessed file ``$(BUILDPATH)/*.p.F90``.

Output from the ``gfortran`` compiler is directed through the ``postprocess.py`` script. This script currently performs the following functions:

* Builds a map from line numbers in the preprocessed source file back to line numbers in the unpreprocessed source, and translates line numbers in any error or warnings messages from the compiler into line numbers in the unpreprocessed source;
* Filters out incorrect warning messages about the lack of scalar finalizers (see GCC `PR58175 <https://gcc.gnu.org/bugzilla/show_bug.cgi?id=58175>`_);
* Checks for unused function attributes in the source and filters out warnings about these unused functions.

Build Files
~~~~~~~~~~~

All files involved in the build process are summarized below.

``scripts/build/buildCode.py``:
   Acts upon ``include`` directives embedded in the source code and generates the corresponding include file (see Section :galacticus-ref:`codeBuildIncludeDirectives`);

``scripts/build/codeDirectivesParse.py``:
   Discovers and parses directives embedded in source files. Generates ``make`` rules for the resulting dependencies, a file providing a mapping of which files contain each directive, and rule files describing how to build each include file resulting from a ``include`` directive;

``scripts/build/includeDependencies.py``:
   Discovers (and generates rules for) dependencies between files arising from the use of "``include``" statements (see Section :galacticus-ref:`buildIncludeDeps`);

``scripts/build/findExecutables.py``:
   Discovers (and generates rules for) files which generate executables (see Section :galacticus-ref:`buildExecutables`);

``scripts/build/libraryDependencies.py``:
   Determines the set of libraries that must be linked with each executable, and ensures they are ordered correctly for static linking;

``scripts/build/moduleDependencies.py``:
   Discovers (and generates rules for) dependencies of module files on their source file statements (see Section :galacticus-ref:`buildModulesProvided`);

``scripts/build/parameterDependencies.py``:
   Determines all parameters upon which a given executable depends (see Section :galacticus-ref:`parameterDependencies`);

``scripts/build/postprocess.py``:
   Postprocesses the output of the ``gfortran`` compiler to remove spurious warnings and provide line numbers in warning/error reports that point into the unpreprocessed source files (see Section :galacticus-ref:`codeBuildPrePostProcess`);

``scripts/build/preprocess.py``:
   Preprocesses Fortran source code to provide various extended functionality (see Section :galacticus-ref:`codeBuildPrePostProcess`);

``scripts/build/useDependencies.py``:
   Discovers (and generates rules for) dependencies between files originating from module ``use`` statements (see Section :galacticus-ref:`buildModulesUsed`);

``$(BUILDPATH)/Makefile_All_Execs``:
   Contains rules describing how to build the executable file for each distinct program---generated by ``scripts/build/findExecutables.py`` (see Section :galacticus-ref:`buildExecutables`);

``$(BUILDPATH)/Makefile_Component_Includes``:
   Describes dependencies on the ``nodeComponent`` class hierarchy module on include files which implement specific functionality for individual node component implementations;

``$(BUILDPATH)/Makefile_Directives``:
   Describes rules for building include files for ``include`` directives;

``$(BUILDPATH)/Makefile_Include_Dependencies``:
   Contains rules describing the dependencies of source file on files included via "``include``" statements---generated by ``scripts/build/includeDependencies.py`` (see Section :galacticus-ref:`buildIncludeDeps`);

``$(BUILDPATH)/Makefile_Module_Dependencies``:
   Contains rules which describe how to build module files from their source files (see Section :galacticus-ref:`buildModulesProvided`);

``$(BUILDPATH)/Makefile_Use_Dependencies``:
   Contains rules which describe dependencies between source file originating from module ``use`` statements (see Section :galacticus-ref:`buildModulesUsed`);

``$(BUILDPATH)/directiveLocations.xml``:
   Contains a map of which files contain each code directive (see Section :galacticus-ref:`buildDiscoveryDirectives`);

``$(BUILDPATH)/*.xml``:
   Contain rules for building files for ``include`` files (see Section :galacticus-ref:`buildDiscoveryDirectives`);

``$(BUILDPATH)/*.d``:
   Contain lists of dependencies on object files, and are accumulated to find the final set of files which must be linked to build each executable (see Section :galacticus-ref:`buildModulesProvided` & Section :galacticus-ref:`buildModulesUsed`);

``$(BUILDPATH)/*.m``:
   Contain lists of modules provided by each object file, and are used during testing whether modules have been updated (see Section :galacticus-ref:`buildModulesProvided`);

``$(BUILDPATH)/*.fl``:
   Contain lists of libraries upon which each object file depends, and are used to determine the final set of libraries which must be linked with each executable (see Section :galacticus-ref:`buildModulesUsed`);

``$(BUILDPATH)/*.gv``:
   Contain lists of source files upon which each file depends and are used in the construction of :term:`GraphViz` representations of file dependencies (see Section :galacticus-ref:`buildModulesProvided` & Section :galacticus-ref:`buildModulesUsed`);

``$(BUILDPATH)/*.Inc``:
   Unpreprocessed include files generated from ``include`` dependencies (see Section :galacticus-ref:`buildDiscoveryDirectives`);

``$(BUILDPATH)/*.inc``:
   Preprocessed include files generated from ``include`` dependencies (see Section :galacticus-ref:`buildDiscoveryDirectives`);

``$(BUILDPATH)/*.p.F90``:
   Preprocessed Fortran source files (see Section :galacticus-ref:`codeBuildPrePostProcess`);

``$(BUILDPATH)/*.parameters.xml``:
   Contain lists of parameters upon which an executable depends (see Section :galacticus-ref:`parameterDependencies`).

Profiling the Build
~~~~~~~~~~~~~~~~~~~

Tooling is provided to help profile the build process. This is useful to identify bottlenecks during build. Profiling can be switched on using the ``compileprof`` build option. For example:

.. code-block:: none

   make -j16 GALACTICUS_BUILD_OPTION=compileprof Galacticus.exe >& build.log

Every command run by ``Make`` will be timed, its peak memory usage measured, and the results output. In the above these outputs are collected to the ``build.log`` file. Result lines begin ``++Task:`` followed by the start and end times of the task, the peak resident set size of the task (in kilobytes; a value of :math:`-1` indicates that memory usage could not be measured, for example because GNU ``time`` is unavailable), and then the command run. While you can look through these manually a script is provided that generates a report from these timing data as a web page. For example:

.. code-block:: none

   ./scripts/build/buildProfiler.py build.log profile.html --durationMinimum 10

will parse timing data from ``build.log`` and generate a web page ``profile.html`` showing the results. In this case the ``--durationMinimum 10`` option specifies that only commands which took 10 or more seconds to run should be included in the report. (Reducing this limit will lead to a very long report.)

The profile report begins by listing the total time taken for the build, together with an estimate of the peak memory used by the build---found by summing, for each second of the build, the peak memory of all included tasks running at that time. Following that is a list of all commands performed ordered by *completion* time. Next to each command is a bar which extends from the start to the end time of the command. Each second of the bar is colored to indicate the degree of build parallelism at that time---green shows maximum parallelism, while red shows no parallelism (i.e. a single command was running at that time).

Finally, each command is assigned a "cost", :math:`\chi`. This is defined as:

.. math::

   \chi = \sum_{i_\mathrm{start}}^{i_\mathrm{end}} N_i^{-1},

where :math:`i_\mathrm{start}` and :math:`i_\mathrm{end}` are the start and end times of the command, and :math:`N_i` is the number of commands running in parallel at time :math:`i`. Cost is therefore higher for commands which run longer and for commands which are executed with less parallelism. A ranked list of commands, from most to least costly is included in the report. The report also shows the peak memory used by each command, together with a ranking of commands from most to least memory-hungry, to help identify the largest memory consumers in the build.
