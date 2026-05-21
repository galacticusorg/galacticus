"""Parsers that walk an already-parsed Fortran source tree.

Each submodule provides a ``parse_*`` (or equivalent) function that takes a
:class:`Galacticus.Build.SourceTree` AST and annotates or extracts a specific
syntactic element:

* :mod:`~Galacticus.Build.SourceTree.Parse.Declarations` -- variable and type
  declarations.
* :mod:`~Galacticus.Build.SourceTree.Parse.Directives` -- ``!![ ... !!]``
  XML directive blocks (the in-tree counterpart of
  :mod:`Galacticus.Build.Directives`, which works directly on the file).
* :mod:`~Galacticus.Build.SourceTree.Parse.ModuleUses` -- ``use`` statements
  and the symbols imported from each module.
* :mod:`~Galacticus.Build.SourceTree.Parse.OpenMP` -- OpenMP directives.
* :mod:`~Galacticus.Build.SourceTree.Parse.Visibilities` -- ``public`` /
  ``private`` access specifiers.

These parsers mirror the corresponding ``perl/Galacticus/Build/SourceTree/Parse/``
modules.
"""
