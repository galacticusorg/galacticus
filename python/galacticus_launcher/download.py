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

import hashlib
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import zipfile
from pathlib import Path

import requests

REPO = "galacticusorg/galacticus"
DATASETS_REPO = "galacticusorg/datasets"

# Release asset listing the SHA-256 of each of the other assets, used to verify what we download.
CHECKSUM_ASSET = "SHA256SUMS"

# Tool executables in the pre-built tools archive whose execute bit must be
# restored after unpacking a zip (macOS tools ship as .zip; tar archives already
# preserve permissions). This list matches the binaries the CI `Build-Tools`
# jobs pack into tools.tar.bz2 / toolsMacOS*.zip; tools built at run time and
# NOT shipped (CosmicEmu's emu.exe, AGN_Spectrum) are intentionally absent.
_TOOL_EXECUTABLE_NAMES = frozenset({
    "camb", "class", "recfast.exe", "cloudy.exe", "autosps.exe", "harmonize",
})

_CHUNK = 1 << 20  # 1 MiB


def asset_url(tag, asset):
    return f"https://github.com/{REPO}/releases/download/{tag}/{asset}"


def source_url(tag):
    """URL of the Galacticus source archive matching `tag`."""
    if tag == "bleeding-edge":
        return f"https://github.com/{REPO}/archive/refs/heads/master.zip"
    return f"https://github.com/{REPO}/archive/refs/tags/{tag}.zip"


def datasets_url(ref="master"):
    """URL of the datasets archive at `ref` (a branch, tag, or commit SHA).

    GitHub's ``archive/<ref>.zip`` endpoint resolves all three, so a pinned
    commit SHA yields exactly that snapshot.
    """
    return f"https://github.com/{DATASETS_REPO}/archive/{ref}.zip"


def resolve_datasets_ref(tag, *, log=print):
    """Resolve which datasets ref to fetch for a release `tag`.

    ``GALACTICUS_DATASETS_REF`` overrides everything. Otherwise a versioned
    release pins datasets to the commit recorded in its ``datasets.ref`` asset
    (written by the publish workflow), so a given package version always fetches
    the same data. Dev/bleeding-edge installs, or a release with no pin, track
    datasets ``master``.
    """
    override = os.environ.get("GALACTICUS_DATASETS_REF")
    if override:
        return override
    if tag and tag != "bleeding-edge":
        pinned = _read_remote_text(asset_url(tag, "datasets.ref"))
        if pinned:
            return pinned.strip().split()[0]
        log(f"  no datasets pin on release {tag}; using datasets master.")
    return "master"


def _read_remote_text(url):
    """Return the text body of `url`, or None if it is absent / unreachable."""
    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            return response.text
    except requests.RequestException:  # pragma: no cover - network
        pass
    return None


def load_checksums(tag, *, log=print):
    """Return {asset name: SHA-256} from the ``SHA256SUMS`` asset of release `tag`, or None if it has none.

    The release workflow writes this file alongside the binary and tools archives (see the `Deploy` job in
    ``.github/workflows/cicd.yml``). It is fetched over HTTPS from the same release as the assets it covers, so it does not
    protect against a compromise of the account or workflow which publishes them -- what it does establish is that the artefacts
    actually received are the ones that release published, so a corrupted, truncated, or substituted asset is caught rather than
    being executed. Releases published before this file existed simply have no checksums; that is reported and provisioning
    continues, since refusing would break every existing install.
    """
    body = _read_remote_text(asset_url(tag, CHECKSUM_ASSET))
    if not body:
        return None
    checksums = {}
    for line in body.splitlines():
        fields = line.split()
        if len(fields) == 2:
            digest, name = fields
            checksums[name.lstrip("*")] = digest.lower()
    return checksums or None


def _sha256(path):
    """Return the SHA-256 of the file at `path`, as lowercase hexadecimal."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def _verify(path, name, checksums, *, log=print):
    """Check the file at `path` against its expected SHA-256, if one is published.

    A mismatch is always fatal, and the offending file is removed: these artefacts are executed, so a file which is not what the
    release published must not be left on disk where a later run could find it and skip the download.
    """
    if checksums is None:
        return
    expected = checksums.get(name)
    if expected is None:
        log(f"  warning: release publishes no checksum for {name}; skipping verification.")
        return
    actual = _sha256(path)
    if actual != expected:
        Path(path).unlink(missing_ok=True)
        raise RuntimeError(
            f"checksum mismatch for {name}:\n"
            f"  expected SHA-256: {expected}\n"
            f"  actual   SHA-256: {actual}\n"
            "The downloaded file is not the one published by this release and has been removed."
        )
    log(f"  verified {name} (SHA-256).")


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

    # Fetch the published checksums once, and verify each downloaded artefact against them below.
    checksums = load_checksums(install.tag, log=log)
    if checksums is None:
        log(f"  note: release {install.tag} publishes no {CHECKSUM_ASSET}; "
            "downloaded artefacts can not be verified against it.")

    if _provision_exec(install, force=force, log=log, checksums=checksums):
        done.append("exec")
    if _provision_datasets(install, force=force, log=log):
        done.append("datasets")
    if _provision_tools(install, force=force, log=log, checksums=checksums):
        done.append("tools")
    if _provision_catalog(install, force=force, log=log):
        done.append("parameter catalog")
    return done


def _provision_catalog(install, *, force, log):
    """Generate the parameter catalog from the just-provisioned source tree.

    The catalog (used by ``galacticus validate``) is NOT a committed/shipped
    artifact; it is generated on demand here so it always matches the installed
    source. Best-effort: a failure is logged but does not fail the install
    (validation then falls back to the executable's ``--dry-run``).
    """
    catalog = Path(install.exec_path) / "parameters.catalog.json"
    generator = Path(install.exec_path) / "scripts" / "build" / "parameterCatalog.py"
    if not generator.is_file():
        return False                       # source tree not present
    if catalog.is_file() and not force:
        return False
    log("Generating parameter catalog ...")
    try:
        subprocess.run(
            [sys.executable, str(generator), str(install.exec_path), str(catalog)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        log(f"warning: could not generate the parameter catalog ({error}); "
            "validation will fall back to the executable's --dry-run.")
        return False
    return catalog.is_file()


def _sentinel(directory, name):
    return Path(directory) / f".galacticus-{name}"


def _provision_exec(install, *, force, log, checksums=None):
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
    # Verify before making the file executable, so that a file which fails the check is never left runnable.
    _verify(binary, install.assets.binary, checksums, log=log)
    binary.chmod(0o755)
    sentinel.write_text(install.tag)
    return True


def _provision_datasets(install, *, force, log):
    sentinel = _sentinel(install.data_path, "datasets")
    if sentinel.exists() and not force:
        return False
    ref = resolve_datasets_ref(install.tag, log=log)
    log(f"Fetching datasets ({ref}) ...")
    with tempfile.TemporaryDirectory() as work:
        archive = Path(work) / "datasets.zip"
        _download(datasets_url(ref), archive, log=log)
        _extract_strip_top(archive, install.data_path)
    sentinel.write_text(ref)
    return True


def _provision_tools(install, *, force, log, checksums=None):
    sentinel = _sentinel(install.tools_path, "tools")
    if sentinel.exists() and not force:
        return False
    log(f"Fetching tools {install.assets.tools} ...")
    with tempfile.TemporaryDirectory() as work:
        archive = Path(work) / install.assets.tools
        _download(asset_url(install.tag, install.assets.tools), archive, log=log)
        # Verify before unpacking: this archive contains executables which are later run.
        _verify(archive, install.assets.tools, checksums, log=log)
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
    """Stream `url` to `dest` with exponential-backoff retries and progress."""
    last_error = None
    for attempt in range(retries):
        try:
            with requests.get(url, stream=True, timeout=60) as response:
                response.raise_for_status()
                total = _content_length(response)
                tmp = Path(str(dest) + ".part")
                progress = _Progress(total, log=log)
                with open(tmp, "wb") as handle:
                    for chunk in response.iter_content(chunk_size=_CHUNK):
                        if chunk:
                            handle.write(chunk)
                            progress.update(len(chunk))
                progress.finish()
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


def _content_length(response):
    """Total size in bytes from the response headers, or None if not advertised."""
    value = response.headers.get("Content-Length")
    try:
        return int(value) if value is not None else None
    except (TypeError, ValueError):  # pragma: no cover - malformed header
        return None


class _Progress:
    """Render download progress.

    On an interactive terminal a single carriage-return-updated bar is drawn to
    stderr.  Otherwise (logs, pipes, non-default ``log``) progress is emitted as
    occasional milestone lines, so it stays readable without a TTY.  Dependency
    free -- no ``tqdm`` -- to keep the launcher's footprint minimal.
    """

    _BAR_WIDTH = 30
    _MILESTONE = 0.10  # log every 10% when not a TTY

    def __init__(self, total, *, log=print):
        self._total = total if total and total > 0 else None
        self._log = log
        self._done = 0
        self._last_render = 0.0
        self._next_milestone = self._MILESTONE
        self._tty = (
            log is print
            and hasattr(sys.stderr, "isatty")
            and sys.stderr.isatty()
        )
        self._active = False

    def update(self, count):
        self._done += count
        if self._tty:
            self._render_bar()
        elif self._total is not None:
            fraction = self._done / self._total
            if fraction >= self._next_milestone:
                while self._next_milestone <= fraction:
                    self._next_milestone += self._MILESTONE
                self._log(f"  ... {int(fraction * 100)}% "
                          f"({_human(self._done)} / {_human(self._total)})")

    def finish(self):
        if self._tty and self._active:
            self._render_bar(force=True)
            sys.stderr.write("\n")
            sys.stderr.flush()

    def _render_bar(self, *, force=False):
        now = time.monotonic()
        if not force and now - self._last_render < 0.1:
            return
        self._last_render = now
        self._active = True
        if self._total is not None:
            fraction = min(1.0, self._done / self._total)
            filled = int(self._BAR_WIDTH * fraction)
            bar = "#" * filled + "-" * (self._BAR_WIDTH - filled)
            text = (f"\r  [{bar}] {int(fraction * 100):3d}% "
                    f"{_human(self._done)} / {_human(self._total)}")
        else:
            text = f"\r  {_human(self._done)} downloaded"
        sys.stderr.write(text)
        sys.stderr.flush()


def _human(num_bytes):
    value = float(num_bytes)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024.0 or unit == "TiB":
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024.0


def _extract(archive, dest, fmt):
    """Unpack `archive` (`fmt` in {'zip','tar.bz2','tar.gz'}) into `dest`.

    Tar archives are unpacked with the ``data`` filter, which refuses members whose paths would escape `dest` (via an absolute
    path, a ``..`` component, or a symlink pointing outside the destination). Without it, unpacking a tar archive lets the archive
    choose where its contents land -- see https://docs.python.org/3/library/tarfile.html#tarfile-extraction-filter . Python 3.14
    makes this the default; setting it explicitly means the behaviour does not depend on the interpreter version. ``zipfile``
    already sanitizes member paths itself, so needs no equivalent.
    """
    dest.mkdir(parents=True, exist_ok=True)
    if fmt == "zip":
        with zipfile.ZipFile(archive) as zf:
            zf.extractall(dest)
    else:
        mode = "r:bz2" if fmt == "tar.bz2" else "r:gz"
        with tarfile.open(archive, mode) as tf:
            tf.extractall(dest, filter="data")


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
