"""Build-system support modules used while compiling Galacticus.

The build entry points (``scripts/build/preprocess.py``,
``scripts/build/libraryInterfaces.py``, and
``scripts/build/codeDependencies.py``) drive the modules in this package to
generate Fortran source from the component hierarchy and to process the
``!![ ... !!]`` XML directives embedded in Galacticus source files.

Submodules:

* :mod:`Galacticus.Build.Components` -- parsing of component-class
  definitions and helpers for emitting their Fortran representation.
* :mod:`Galacticus.Build.Dependencies` -- module-level dependency extraction
  used to build per-source dependency lists.
* :mod:`Galacticus.Build.Directives` -- on-disk extraction of ``!![ ... !!]``
  XML directive blocks from a Fortran source file.
* :mod:`Galacticus.Build.FileChanges` -- utilities for writing files only when
  their content has changed (so dependent build steps stay cached).
* :mod:`Galacticus.Build.FortranUtils` -- small Fortran-aware string helpers
  (line continuation, comment stripping, etc.).
* :mod:`Galacticus.Build.SourceTree` -- tree-based parser for Fortran source,
  used when a directive's context within the file is needed.
* :mod:`Galacticus.Build.StateStorables` -- support for the state-storable
  serialization metadata.
"""
