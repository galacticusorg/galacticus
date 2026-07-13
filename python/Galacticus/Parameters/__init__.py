"""Typed parameter catalog for Galacticus.

Andrew Benson (2026)

This package extracts a static, machine-readable catalog of the input
parameters accepted by every ``functionClass`` implementation, by harvesting
the ``<inputParameter>`` and ``<objectBuilder>`` directives embedded in the
Fortran source and inferring a type for each parameter.  The catalog is the
foundation for typed validation of parameter files (and, longer term, a typed
Python configuration layer).

Modules
-------
inference : Fortran-literal and declaration-based type inference.
catalog   : Source-tree harvesting that assembles the catalog dict.
"""

from Galacticus.Parameters.catalog import build_catalog
from Galacticus.Parameters.inference import infer_parameter_type

__all__ = ['build_catalog', 'infer_parameter_type']
