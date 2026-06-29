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

    Honors ``GALACTICUS_PARAMETER_CATALOG`` first; then the exec path root (where
    a managed install generates it at provision time); then the default build
    output dir ``work/build/`` (where ``make parameters-catalog`` writes it for a
    build-from-source install). A non-default ``BUILDPATH`` is supported via the
    env override.
    """
    override = os.environ.get("GALACTICUS_PARAMETER_CATALOG")
    if override and Path(override).is_file():
        return Path(override)
    exec_path = Path(install.exec_path)
    for candidate in (exec_path / "parameters.catalog.json",
                      exec_path / "work" / "build" / "parameters.catalog.json"):
        if candidate.is_file():
            return candidate
    return None


def _galacticus_python_dir(install):
    return Path(install.exec_path) / "python"


def _ensure_importable(install):
    """Put the exec-path ``python/`` tree (matching the binary) on ``sys.path``."""
    python_dir = _galacticus_python_dir(install)
    if python_dir.is_dir() and str(python_dir) not in sys.path:
        sys.path.insert(0, str(python_dir))


def _import_resolver(install):
    _ensure_importable(install)
    from Galacticus.Parameters import resolve
    return resolve


def resolve_to_file(param_file, install, change_files=(), *, output,
                    conditionals=True):
    """Resolve ``param_file`` (+ change files) and write ``output``.

    Raises ``RuntimeError`` (which the CLI reports) on a resolution error or if
    the resolver cannot be imported.
    """
    try:
        resolve = _import_resolver(install)
    except ImportError as error:
        raise RuntimeError(f"parameter resolver unavailable: {error}")
    try:
        resolve.resolve_file(param_file, list(change_files), output=output,
                             conditionals=conditionals)
    except resolve.ResolveError as error:
        raise RuntimeError(f"cannot resolve {param_file}: {error}")


def validate(param_file, install, *, structural=False, change_files=()):
    """Validate `param_file` for `install`; return a :class:`Result`.

    Prefers the catalog-aware Python checker (run on the fully RESOLVED tree --
    XInclude, change files, and conditionals applied -- so it checks the
    structure Galacticus will actually build); falls back to the binary's
    ``--dry-run``.  ``ok`` is False only on an error-level finding or a parse /
    resolve / dry-run failure (type-level findings are warnings and do not fail).
    """
    catalog_path = find_catalog(install)
    if catalog_path is not None:
        return _validate_with_catalog(param_file, install, catalog_path, structural,
                                      change_files)
    return _validate_with_dry_run(param_file, install, change_files)


def _validate_with_catalog(param_file, install, catalog_path, structural, change_files):
    import json
    import tempfile

    _ensure_importable(install)
    try:
        from Galacticus.Parameters.validate import validate_file
    except ImportError as error:
        return Result(True, [], "catalog-unavailable", str(error))

    # Resolve first so validation sees the structure Galacticus will build. If the
    # resolver is unavailable (very old install), fall back to the raw file --
    # validate_file still expands XInclude itself.
    target = str(param_file)
    resolved_tmp = None
    try:
        resolve = _import_resolver(install)
    except ImportError:
        resolve = None
    if resolve is not None:
        handle = tempfile.NamedTemporaryFile(suffix=".xml", delete=False)
        handle.close()
        try:
            resolve.resolve_file(param_file, list(change_files), output=handle.name)
        except resolve.ResolveError as error:
            os.unlink(handle.name)
            return Result(False,
                          [Finding("error", "resolve", str(param_file), str(error))],
                          "catalog", str(error))
        target, resolved_tmp = handle.name, handle.name

    try:
        with open(catalog_path) as handle:
            catalog = json.load(handle)
        findings_raw, parse_error = validate_file(
            target, catalog, structural=structural)
    finally:
        if resolved_tmp is not None:
            os.unlink(resolved_tmp)

    if parse_error:
        return Result(False, [Finding("error", "parse", str(param_file), parse_error)],
                      "catalog", parse_error)
    findings = [Finding(f.level, f.kind, f.path, f.message) for f in findings_raw]
    ok = not any(f.level == "error" for f in findings)
    return Result(ok, findings, "catalog", str(catalog_path))


def _validate_with_dry_run(param_file, install, change_files=()):
    env = install.environ()
    try:
        completed = subprocess.run(
            [str(install.binary), str(param_file), *change_files, "--dry-run"],
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
