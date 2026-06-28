"""End-user launcher for Galacticus.

Galacticus is a Fortran/C++ code, but its pre-built static binaries plus the
``datasets`` repository and run-time ``tools`` are all published as GitHub
release assets.  This package turns ``pip install galacticus`` into a working
model run: it downloads the right artefacts for the host platform into a
per-user location, sets the environment variables the binary reads
(``GALACTICUS_EXEC_PATH``, ``GALACTICUS_DATA_PATH``, ``GALACTICUS_TOOLS_PATH``,
``GALACTICUS_DYNAMIC_DATA_PATH``), and dispatches parameter files to the
executable.

The ``galacticus`` console script (see :mod:`galacticus_launcher.cli`) is the
entry point.  The launcher also works as a thin dispatch front-end on top of a
non-pip install (a git clone built from source): when a usable build is already
present in the environment, downloads are skipped entirely (see
:mod:`galacticus_launcher.paths`).
"""

from importlib import metadata as _metadata

try:
    __version__ = _metadata.version("galacticus")
except _metadata.PackageNotFoundError:  # pragma: no cover - source checkouts
    __version__ = "0.0.0"

__all__ = ["__version__"]
