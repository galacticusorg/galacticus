"""LibraryInterfaces — Python package for Galacticus C/Fortran/Python interface generation.

Andrew Benson (ported to Python with assistance from Claude 2026)

Sub-modules:
  ArgSpec   — ArgSpec dataclass (intermediate representation for a single argument)
  Pipeline  — four pipeline stages that progressively enrich ArgSpec objects
  Emitters  — ten emitter functions that render ArgSpec lists to code strings
"""
