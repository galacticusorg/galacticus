"""Parameter-file validation hook, invoked before dispatch.

This is intentionally a thin, swappable layer.  Today it reuses what already
exists in the tree:

* the catalog-aware Python checker
  :func:`Galacticus.Parameters.validate.validate_file`, when a
  ``parameters.catalog.json`` can be found; and
* otherwise the executable's own ``--dry-run`` mode, which parses and
  structurally validates the parameter file without running a model.

The public surface is a single :func:`validate` call returning a
:class:`Result`, so a richer middle-layer validator can replace the
implementation without any change to :mod:`galacticus_launcher.cli`.
"""

import os
import subprocess
import sys
from collections import namedtuple
from pathlib import Path

# A finding mirrors Galacticus.Parameters.validate.Finding so callers can treat
# both code paths uniformly.
Finding = namedtuple("Finding", ["level", "kind", "path", "message"])

Result = namedtuple("Result", ["ok", "findings", "method", "detail"])


def find_catalog(install):
    """Locate a ``parameters.catalog.json``, or None.

    Honors ``GALACTICUS_PARAMETER_CATALOG`` first, then looks beside the
    executable's exec path.
    """
    override = os.environ.get("GALACTICUS_PARAMETER_CATALOG")
    if override and Path(override).is_file():
        return Path(override)
    candidate = Path(install.exec_path) / "parameters.catalog.json"
    return candidate if candidate.is_file() else None


def validate(param_file, install, *, structural=False):
    """Validate `param_file` for `install`; return a :class:`Result`.

    Prefers the catalog-aware Python checker; falls back to the binary's
    ``--dry-run``.  ``ok`` is False only on an error-level finding or a parse /
    dry-run failure (type-level findings are warnings and do not fail).
    """
    catalog_path = find_catalog(install)
    if catalog_path is not None:
        return _validate_with_catalog(param_file, install, catalog_path, structural)
    return _validate_with_dry_run(param_file, install)


def _validate_with_catalog(param_file, install, catalog_path, structural):
    import json

    # The checker lives under the exec path's python/ tree; make it importable.
    python_dir = Path(install.exec_path) / "python"
    if python_dir.is_dir() and str(python_dir) not in sys.path:
        sys.path.insert(0, str(python_dir))
    try:
        from Galacticus.Parameters.validate import validate_file
    except ImportError as error:
        return Result(True, [], "catalog-unavailable", str(error))

    with open(catalog_path) as handle:
        catalog = json.load(handle)
    findings_raw, parse_error = validate_file(
        str(param_file), catalog, structural=structural
    )
    if parse_error:
        return Result(False, [Finding("error", "parse", str(param_file), parse_error)],
                      "catalog", parse_error)
    findings = [Finding(f.level, f.kind, f.path, f.message) for f in findings_raw]
    ok = not any(f.level == "error" for f in findings)
    return Result(ok, findings, "catalog", str(catalog_path))


def _validate_with_dry_run(param_file, install):
    env = install.environ()
    try:
        completed = subprocess.run(
            [str(install.binary), str(param_file), "--dry-run"],
            env=env, capture_output=True, text=True,
        )
    except OSError as error:
        return Result(False, [Finding("error", "exec", str(param_file), str(error))],
                      "dry-run", str(error))
    if completed.returncode == 0:
        return Result(True, [], "dry-run", "ok")
    detail = (completed.stderr or completed.stdout or "").strip()
    return Result(False, [Finding("error", "dry-run", str(param_file), detail)],
                  "dry-run", detail)
