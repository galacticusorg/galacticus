"""Galacticus Python support package.

Hosts the build-system and constraint-pipeline modules that were ported from
the legacy ``perl/Galacticus/`` tree:

* :mod:`Galacticus.Build` -- code generation and source-tree processing used
  while compiling Galacticus (component hierarchies, directive extraction,
  Fortran source-tree walks, dependency tracking).
* :mod:`Galacticus.Constraints` -- helpers for the dark-matter constraint
  pipeline (MCMC parameter chains, simulation metadata).
* :mod:`Galacticus._logging` -- shared logging configuration used by the
  build pipeline and constraint scripts.

This package is distinct from the generated ``galacticus`` module produced by
``make libgalacticus.so``, which provides the ``ctypes`` interface to the
shared library (see ``doc/Python_Interface.tex``).
"""
