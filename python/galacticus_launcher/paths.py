"""Resolve where Galacticus and its data live, and the env it runs under.

Two storage roots are used for a managed (pip) install, chosen so that
refreshing the install and purging the cache are independent operations:

* **data root** (``platformdirs.user_data_dir``) -- durable, managed by
  ``galacticus install``/``update``: the executable, the source tree (for
  ``parameters/`` + ``aux/`` + ``scripts/``), the ``datasets`` repository, and
  the pre-built run-time ``tools``.  Laid out per release tag so multiple
  versions can coexist: ``<data>/<tag>/{exec,datasets,tools}``.
* **cache root** (``platformdirs.user_cache_dir``) -- regenerable compute cache
  (transfer functions, SSP spectra, cooling tables, ...): ``<cache>/<tag>``.
  Safe to delete; ``galacticus clean`` only ever touches this.

Tools deliberately live under the *data* root (not the cache): a binary-only
pip user has no compilers, so if the pre-built tools were lost to a cache purge
Galacticus would try -- and fail -- to build CAMB/Cloudy/etc. from source.
``GALACTICUS_TOOLS_PATH`` keeps them out of the purgeable cache.

Resolution order (the first that matches wins), which lets the same launcher
drive a non-pip, build-from-source install with downloads skipped:

1. ``GALACTICUS_HOME`` -- a Galacticus build/clone tree containing
   ``Galacticus.exe``; data taken from ``GALACTICUS_DATA_PATH`` if set, else
   from a sibling ``datasets`` directory, else the managed data root.
2. ``GALACTICUS_EXEC_PATH`` + ``GALACTICUS_DATA_PATH`` already set in the
   environment with the executable present -- use them as-is.
3. The managed download cache (default for pip users).
"""

import os
from collections import namedtuple
from pathlib import Path

import platformdirs

from . import platforms

APP_NAME = "galacticus"

# Source of the resolved install, for reporting by `galacticus info`.
SOURCE_ENVIRONMENT = "environment"   # GALACTICUS_EXEC_PATH/DATA_PATH preset
SOURCE_HOME = "home"                 # GALACTICUS_HOME build/clone tree
SOURCE_MANAGED = "managed"           # downloaded into platformdirs cache


class Install(namedtuple(
    "Install",
    ["source", "tag", "exec_path", "data_path", "tools_path",
     "dynamic_path", "binary", "assets"],
)):
    """A resolved Galacticus install and the environment it runs under.

    Paths are :class:`pathlib.Path`.  ``assets`` is the
    :class:`~galacticus_launcher.platforms.PlatformAssets` for the host, or
    ``None`` for a local (non-managed) install where assets are never fetched.
    ``managed`` is True only when the launcher owns the install and may
    download into it.
    """

    @property
    def managed(self):
        return self.source == SOURCE_MANAGED

    @property
    def tools_separated(self):
        """True if the pre-built tools live on a path distinct from the dynamic
        (regenerable) data path, so the dynamic path can be purged safely.

        Always true for a managed install (tools under the data root). For an
        environment/home install it is true only when ``GALACTICUS_TOOLS_PATH``
        was set to something other than the dynamic path; otherwise tools are
        colocated under the dynamic path (the binary's default) and purging the
        dynamic path would also delete them.
        """
        if self.managed:
            return True
        return self.tools_path is not None and self.tools_path != self.dynamic_path

    def environ(self, base=None):
        """Return a copy of `base` (default ``os.environ``) with the Galacticus
        path variables set to this install's locations.

        A managed install imposes the full layout (exec/data/tools/dynamic). An
        environment/home install only fixes exec/data and otherwise passes the
        user's environment through unchanged, so ``GALACTICUS_TOOLS_PATH`` /
        ``GALACTICUS_DYNAMIC_DATA_PATH`` stay exactly as the user set them (or
        unset -> the binary's own defaults) and ``run`` behaves identically to
        invoking the binary directly.
        """
        env = dict(os.environ if base is None else base)
        env["GALACTICUS_EXEC_PATH"] = str(self.exec_path)
        if self.data_path is not None:
            env["GALACTICUS_DATA_PATH"] = str(self.data_path)
        if self.managed:
            env["GALACTICUS_TOOLS_PATH"] = str(self.tools_path)
            env["GALACTICUS_DYNAMIC_DATA_PATH"] = str(self.dynamic_path)
        return env


def release_tag(version=None):
    """Resolve the GitHub release tag to fetch assets from.

    ``GALACTICUS_RELEASE_TAG`` overrides everything.  Otherwise a clean
    ``X.Y.Z`` package version (other than the ``0.0.0`` development placeholder)
    maps to the ``vX.Y.Z`` release; anything else -- dev/local checkouts -- maps
    to the rolling ``bleeding-edge`` release.
    """
    override = os.environ.get("GALACTICUS_RELEASE_TAG")
    if override:
        return override
    if version is None:
        from . import __version__ as version
    parts = version.split(".")
    is_final = (
        version != "0.0.0"
        and len(parts) == 3
        and all(part.isdigit() for part in parts)
    )
    return f"v{version}" if is_final else "bleeding-edge"


def data_root(tag):
    return Path(platformdirs.user_data_dir(APP_NAME)) / tag


def cache_root(tag):
    return Path(platformdirs.user_cache_dir(APP_NAME)) / tag


def _has_binary(directory):
    binary = Path(directory) / platforms.LOCAL_BINARY_NAME
    return binary if binary.is_file() else None


def resolve(version=None):
    """Resolve the active install per the documented order above.

    Never performs I/O beyond ``stat`` checks; for a managed install it returns
    the *intended* layout whether or not the artefacts have been downloaded yet
    (callers provision via :mod:`galacticus_launcher.download`).
    """
    # 1. GALACTICUS_HOME -- an explicit build/clone tree.
    home = os.environ.get("GALACTICUS_HOME")
    if home:
        home = Path(home)
        binary = _has_binary(home) or (home / platforms.LOCAL_BINARY_NAME)
        data = os.environ.get("GALACTICUS_DATA_PATH")
        if data:
            data = Path(data)
        elif (home.parent / "datasets" / "static").is_dir():
            data = home.parent / "datasets"
        else:
            data = None
        return Install(
            source=SOURCE_HOME, tag=None, exec_path=home, data_path=data,
            tools_path=_env_path("GALACTICUS_TOOLS_PATH"),
            dynamic_path=_effective_dynamic(data),
            binary=binary, assets=None,
        )

    # 2. Environment already configured for an existing (source) install.
    exec_env = os.environ.get("GALACTICUS_EXEC_PATH")
    data_env = os.environ.get("GALACTICUS_DATA_PATH")
    if exec_env and data_env:
        binary = _has_binary(exec_env)
        if binary is not None:
            return Install(
                source=SOURCE_ENVIRONMENT, tag=None, exec_path=Path(exec_env),
                data_path=Path(data_env),
                tools_path=_env_path("GALACTICUS_TOOLS_PATH"),
                dynamic_path=_effective_dynamic(Path(data_env)),
                binary=binary, assets=None,
            )

    # 3. Managed download install.
    assets = platforms.detect()
    tag = release_tag(version)
    data = data_root(tag)
    cache = cache_root(tag)
    return Install(
        source=SOURCE_MANAGED, tag=tag,
        exec_path=data / "exec",
        data_path=data / "datasets",
        tools_path=data / "tools",
        dynamic_path=cache / "dynamic",
        binary=data / "exec" / platforms.LOCAL_BINARY_NAME,
        assets=assets,
    )


def _env_path(name):
    value = os.environ.get(name)
    return Path(value) if value else None


def _effective_dynamic(data_path):
    """The dynamic (regenerable) data path the binary will actually use:
    ``GALACTICUS_DYNAMIC_DATA_PATH`` if set, else ``<data_path>/dynamic`` (the
    binary's default), else None. Used so ``info``/``clean`` reflect where the
    regenerable data really lives for an environment/home install."""
    explicit = _env_path("GALACTICUS_DYNAMIC_DATA_PATH")
    if explicit is not None:
        return explicit
    return (data_path / "dynamic") if data_path is not None else None
