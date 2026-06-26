"""Fetch and unpack Galacticus release artefacts into a managed install.

Provisioning is idempotent: each component drops a sentinel file once complete,
and a present sentinel means the component is skipped.  Only managed installs
are ever provisioned -- a local build/clone is used in place, untouched.

Four components make up a runnable install:

* **exec**     -- the repository source archive (whole tree), giving
  ``GALACTICUS_EXEC_PATH`` its ``parameters/``, ``aux/`` and ``scripts/``,
  plus the platform executable copied in as ``Galacticus.exe``.
* **datasets** -- the ``galacticusorg/datasets`` archive (static data root).
* **tools**    -- the pre-built run-time tools archive.  Its entries are
  prefixed ``dynamic/`` (the location the un-relocated binary expects); we strip
  that prefix so the contents land directly under ``GALACTICUS_TOOLS_PATH``.
"""

import os
import shutil
import tarfile
import tempfile
import time
import zipfile
from pathlib import Path

import requests

REPO = "galacticusorg/galacticus"
DATASETS_REPO = "galacticusorg/datasets"

# Tool executables whose execute bit must be restored after unpacking a zip
# (macOS tools ship as .zip; tar archives already preserve permissions).
_TOOL_EXECUTABLE_NAMES = frozenset({
    "camb", "class", "recfast.exe", "cloudy.exe", "autosps.exe",
    "harmonize", "ransack", "emu.exe",
})

_CHUNK = 1 << 20  # 1 MiB


def asset_url(tag, asset):
    return f"https://github.com/{REPO}/releases/download/{tag}/{asset}"


def source_url(tag):
    """URL of the Galacticus source archive matching `tag`."""
    if tag == "bleeding-edge":
        return f"https://github.com/{REPO}/archive/refs/heads/master.zip"
    return f"https://github.com/{REPO}/archive/refs/tags/{tag}.zip"


def datasets_url():
    """URL of the datasets archive (``GALACTICUS_DATASETS_REF`` selects a ref)."""
    ref = os.environ.get("GALACTICUS_DATASETS_REF", "master")
    return f"https://github.com/{DATASETS_REPO}/archive/refs/heads/{ref}.zip"


def provision(install, *, force=False, log=print):
    """Download and unpack any missing components of a managed `install`.

    Returns the list of component names that were (re)provisioned.  Raises if
    `install` is not managed.
    """
    if not install.managed:
        raise ValueError("only managed installs can be provisioned")
    if install.assets is None:  # pragma: no cover - managed always has assets
        raise ValueError("managed install is missing platform assets")

    done = []
    for path in (install.exec_path, install.data_path, install.tools_path,
                 install.dynamic_path):
        path.mkdir(parents=True, exist_ok=True)

    if _provision_exec(install, force=force, log=log):
        done.append("exec")
    if _provision_datasets(install, force=force, log=log):
        done.append("datasets")
    if _provision_tools(install, force=force, log=log):
        done.append("tools")
    return done


def _sentinel(directory, name):
    return Path(directory) / f".galacticus-{name}"


def _provision_exec(install, *, force, log):
    sentinel = _sentinel(install.exec_path, "exec")
    if sentinel.exists() and not force:
        return False
    log(f"Fetching Galacticus source ({install.tag}) ...")
    with tempfile.TemporaryDirectory() as work:
        archive = Path(work) / "source.zip"
        _download(source_url(install.tag), archive, log=log)
        _extract_strip_top(archive, install.exec_path)
    log(f"Fetching executable {install.assets.binary} ...")
    binary = install.exec_path / "Galacticus.exe"
    _download(asset_url(install.tag, install.assets.binary), binary, log=log)
    binary.chmod(0o755)
    sentinel.write_text(install.tag)
    return True


def _provision_datasets(install, *, force, log):
    sentinel = _sentinel(install.data_path, "datasets")
    if sentinel.exists() and not force:
        return False
    log("Fetching datasets ...")
    with tempfile.TemporaryDirectory() as work:
        archive = Path(work) / "datasets.zip"
        _download(datasets_url(), archive, log=log)
        _extract_strip_top(archive, install.data_path)
    sentinel.write_text(os.environ.get("GALACTICUS_DATASETS_REF", "master"))
    return True


def _provision_tools(install, *, force, log):
    sentinel = _sentinel(install.tools_path, "tools")
    if sentinel.exists() and not force:
        return False
    log(f"Fetching tools {install.assets.tools} ...")
    with tempfile.TemporaryDirectory() as work:
        archive = Path(work) / install.assets.tools
        _download(asset_url(install.tag, install.assets.tools), archive, log=log)
        staging = Path(work) / "unpacked"
        _extract(archive, staging, install.assets.tools_format)
        # Tools archives are rooted at "dynamic/"; lift that subtree up so the
        # contents sit directly under GALACTICUS_TOOLS_PATH.
        root = staging / "dynamic"
        source = root if root.is_dir() else staging
        for child in source.iterdir():
            destination = install.tools_path / child.name
            if destination.exists():
                shutil.rmtree(destination) if destination.is_dir() else destination.unlink()
            shutil.move(str(child), str(destination))
    _restore_executable_bits(install.tools_path)
    sentinel.write_text(install.tag)
    return True


def _download(url, dest, *, log=print, retries=4):
    """Stream `url` to `dest` with exponential-backoff retries."""
    last_error = None
    for attempt in range(retries):
        try:
            with requests.get(url, stream=True, timeout=60) as response:
                response.raise_for_status()
                tmp = Path(str(dest) + ".part")
                with open(tmp, "wb") as handle:
                    for chunk in response.iter_content(chunk_size=_CHUNK):
                        if chunk:
                            handle.write(chunk)
                tmp.replace(dest)
            return
        except (requests.RequestException, OSError) as error:  # pragma: no cover - network
            last_error = error
            if attempt == retries - 1:
                break
            wait = 2 ** (attempt + 1)
            log(f"  download failed ({error}); retrying in {wait}s ...")
            time.sleep(wait)
    raise RuntimeError(f"failed to download {url}: {last_error}")


def _extract(archive, dest, fmt):
    """Unpack `archive` (`fmt` in {'zip','tar.bz2','tar.gz'}) into `dest`."""
    dest.mkdir(parents=True, exist_ok=True)
    if fmt == "zip":
        with zipfile.ZipFile(archive) as zf:
            zf.extractall(dest)
    else:
        mode = "r:bz2" if fmt == "tar.bz2" else "r:gz"
        with tarfile.open(archive, mode) as tf:
            tf.extractall(dest)


def _extract_strip_top(archive, dest):
    """Extract a GitHub source/datasets zip, stripping its single top-level
    directory (``galacticus-master/`` / ``datasets-master/``) so contents land
    directly in `dest`."""
    dest.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as work:
        with zipfile.ZipFile(archive) as zf:
            zf.extractall(work)
        entries = [child for child in Path(work).iterdir()]
        top = entries[0] if len(entries) == 1 and entries[0].is_dir() else Path(work)
        for child in top.iterdir():
            destination = dest / child.name
            if destination.exists():
                shutil.rmtree(destination) if destination.is_dir() else destination.unlink()
            shutil.move(str(child), str(destination))


def _restore_executable_bits(tools_path):
    for dirpath, _dirnames, filenames in os.walk(tools_path):
        for name in filenames:
            if name in _TOOL_EXECUTABLE_NAMES or name.endswith(".exe"):
                file_path = Path(dirpath) / name
                try:
                    file_path.chmod(file_path.stat().st_mode | 0o111)
                except OSError:  # pragma: no cover - best effort
                    pass
