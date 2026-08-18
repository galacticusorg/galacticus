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
  This is by far the largest artefact, and models which use none of the tools it
  carries (CAMB, CLASS, Cloudy, ...) do not need it, so it is the one component
  which can be skipped -- see :func:`provision`.
* **parameter catalog** -- the machine-readable catalog of accepted input
  parameters used by ``galacticus validate``, downloaded from the release when it
  publishes one and generated from the installed source otherwise.

Archives are staged alongside their destination rather than in the system
temporary directory, so that unpacking is followed by a rename rather than a
copy across filesystems, and so that a multi-gigabyte download does not have to
fit in a (often small, often ``tmpfs``) ``/tmp``.
"""

import hashlib
import json
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

# Release asset carrying the pre-built parameter catalog (see `_provision_catalog`).
CATALOG_ASSET = "parameters.catalog.json"

# Sentinel recording that the tools archive is deliberately not wanted.
TOOLS_SKIPPED = "tools-skipped"

# Tool executables in the pre-built tools archive whose execute bit must be
# restored after unpacking a zip (macOS tools ship as .zip; tar archives already
# preserve permissions). This list matches the binaries the CI `Build-Tools`
# jobs pack into tools.tar.bz2 / toolsMacOS*.zip; tools built at run time and
# NOT shipped (CosmicEmu's emu.exe, AGN_Spectrum) are intentionally absent.
_TOOL_EXECUTABLE_NAMES = frozenset({
    "camb", "class", "recfast.exe", "cloudy.exe", "autosps.exe", "harmonize",
})

_CHUNK = 1 << 20  # 1 MiB

# Prefix for the staging directories archives are downloaded into and unpacked
# in. They sit beside their destination rather than in the system temporary
# directory (see the module docstring), so name them recognizably: an install
# killed part way through leaves one behind, and it should be obvious what it is.
_STAGING = ".galacticus-staging-"


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


def provision(install, *, force=False, log=print, tools=None):
    """Download and unpack any missing components of a managed `install`.

    `tools` selects whether the (large) pre-built tools archive is fetched:
    ``True`` fetches it, ``False`` skips it and records that choice, and ``None``
    (the default) fetches it unless an earlier ``--no-tools`` install recorded
    that it is not wanted.

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
    if _tools_wanted(install, tools, log=log) and \
            _provision_tools(install, force=force, log=log, checksums=checksums):
        done.append("tools")
    if _provision_catalog(install, force=force, log=log, checksums=checksums):
        done.append("parameter catalog")
    return done


def tools_installed(install):
    """True if the pre-built tools archive has been unpacked for `install`."""
    return install.tools_path is not None and _sentinel(install.tools_path, "tools").exists()


def tools_skipped(install):
    """True if the tools archive was deliberately skipped for `install`."""
    return install.tools_path is not None and _sentinel(install.tools_path, TOOLS_SKIPPED).exists()


def _tools_wanted(install, tools, *, log):
    """Decide whether to fetch the tools archive, recording an explicit choice.

    The choice has to persist, because `run` and `validate` provision on demand:
    without a marker on disk the very next command would download the archive the
    user just asked to skip.  Asking for tools clears the marker; asking to skip
    them when they are already installed leaves them alone (nothing is deleted).
    """
    marker = _sentinel(install.tools_path, TOOLS_SKIPPED)
    if tools is None:
        return not marker.exists()          # silent: this runs on every `run`
    if tools:
        marker.unlink(missing_ok=True)
        return True
    if tools_installed(install):
        log("  tools are already installed; leaving them in place.")
        return False
    marker.write_text(install.tag)
    log(f"Skipping {install.assets.tools}; models needing CAMB, CLASS, Cloudy, "
        "or the other pre-built tools will not run.\n"
        "  Run `galacticus install --tools` to add them later.")
    return False


def _provision_catalog(install, *, force, log, checksums=None):
    """Install the parameter catalog used by ``galacticus validate``.

    The catalog is not a committed artifact -- it is derived from the source
    tree -- but each release publishes one built from the very commit it was cut
    from, so it is downloaded when available.  Generating it instead costs tens
    of seconds of CPU on a modest machine, which would otherwise be paid on the
    first ``galacticus run``.  A release which publishes no catalog (or any
    failure to fetch one) falls back to generating it locally, and that in turn
    is best effort: without a catalog, validation falls back to the executable's
    ``--dry-run``.
    """
    catalog = Path(install.exec_path) / CATALOG_ASSET
    if catalog.is_file() and not force:
        return False
    if _fetch_catalog(install, catalog, log=log, checksums=checksums):
        return True
    return _generate_catalog(install, catalog, log=log)


def _fetch_catalog(install, catalog, *, log, checksums):
    """Download the pre-built parameter catalog published by the release.

    Returns True if `catalog` was written.  Every failure mode -- the release
    publishing no catalog, a download error, a checksum mismatch, or a file which
    is not valid JSON -- returns False so the caller can generate one instead.
    Unlike the executable and the tools, the catalog is data and is never run, so
    a bad copy is a reason to fall back rather than to abort the install; the
    generator reproduces the same catalog from the source tree just installed.
    """
    log("Fetching parameter catalog ...")
    staged = Path(str(catalog) + ".part-download")
    try:
        if not _download(asset_url(install.tag, CATALOG_ASSET), staged, log=log,
                         missing_ok=True):
            log(f"  release {install.tag} publishes no {CATALOG_ASSET}.")
            return False
        _verify(staged, CATALOG_ASSET, checksums, log=log)
        with open(staged) as handle:
            json.load(handle)
    except (RuntimeError, OSError, ValueError) as error:
        log(f"  could not use the published parameter catalog ({error}).")
        staged.unlink(missing_ok=True)
        return False
    staged.replace(catalog)
    return True


def _generate_catalog(install, catalog, *, log):
    """Generate the parameter catalog from the installed source tree."""
    generator = Path(install.exec_path) / "scripts" / "build" / "parameterCatalog.py"
    if not generator.is_file():
        return False                       # source tree not present
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
    with tempfile.TemporaryDirectory(dir=install.exec_path.parent, prefix=_STAGING) as work:
        archive = Path(work) / "source.zip"
        _download(source_url(install.tag), archive, log=log)
        _extract_strip_top(archive, install.exec_path, log=log)
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
    with tempfile.TemporaryDirectory(dir=install.data_path.parent, prefix=_STAGING) as work:
        archive = Path(work) / "datasets.zip"
        _download(datasets_url(ref), archive, log=log)
        _extract_strip_top(archive, install.data_path, log=log)
    sentinel.write_text(ref)
    return True


def _provision_tools(install, *, force, log, checksums=None):
    sentinel = _sentinel(install.tools_path, "tools")
    if sentinel.exists() and not force:
        return False
    log(f"Fetching tools {install.assets.tools} ...")
    with tempfile.TemporaryDirectory(dir=install.tools_path.parent, prefix=_STAGING) as work:
        archive = Path(work) / install.assets.tools
        _download(asset_url(install.tag, install.assets.tools), archive, log=log)
        # Verify before unpacking: this archive contains executables which are later run.
        _verify(archive, install.assets.tools, checksums, log=log)
        staging = Path(work) / "unpacked"
        log(f"Unpacking {install.assets.tools} ...")
        _extract(archive, staging, install.assets.tools_format, log=log)
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


def _download(url, dest, *, log=print, retries=4, missing_ok=False):
    """Stream `url` to `dest` with exponential-backoff retries and progress.

    Returns True once `dest` is written.  With `missing_ok`, a 404 -- the release
    simply does not publish this asset -- returns False immediately instead of
    retrying an absence through the full backoff schedule.
    """
    last_error = None
    for attempt in range(retries):
        try:
            with requests.get(url, stream=True, timeout=60) as response:
                if missing_ok and response.status_code == 404:
                    return False
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
            return True
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

    def set_position(self, position):
        """Advance to an absolute `position`, in bytes; never moves backwards."""
        if position > self._done:
            self.update(position - self._done)

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


def _extract(archive, dest, fmt, *, log=print):
    """Unpack `archive` (`fmt` in {'zip','tar.bz2','tar.gz'}) into `dest`.

    Tar archives are unpacked with the ``data`` filter, which refuses members whose paths would escape `dest` (via an absolute
    path, a ``..`` component, or a symlink pointing outside the destination). Without it, unpacking a tar archive lets the archive
    choose where its contents land -- see https://docs.python.org/3/library/tarfile.html#tarfile-extraction-filter . Python 3.14
    makes this the default; setting it explicitly means the behaviour does not depend on the interpreter version. ``zipfile``
    already sanitizes member paths itself, so needs no equivalent.

    Unpacking the release archives takes long enough to look like a hang, so both
    formats report progress through the same bar the downloads use.
    """
    dest.mkdir(parents=True, exist_ok=True)
    if fmt == "zip":
        _extract_zip(archive, dest, log=log)
    else:
        _extract_tar(archive, dest, "r:bz2" if fmt == "tar.bz2" else "r:gz", log=log)


def _extract_zip(archive, dest, *, log=print):
    """Unpack a zip member by member, reporting progress by unpacked bytes.

    Equivalent to ``ZipFile.extractall`` -- which is itself this loop -- other
    than the progress reporting; member paths are sanitized by ``zipfile`` in
    either case.  The member list comes from the central directory, so the total
    is known up front without reading the archive.
    """
    with zipfile.ZipFile(archive) as zf:
        members = zf.infolist()
        progress = _Progress(sum(member.file_size for member in members), log=log)
        for member in members:
            zf.extract(member, dest)
            progress.update(member.file_size)
        progress.finish()


def _extract_tar(archive, dest, mode, *, log=print):
    """Unpack a compressed tar with the ``data`` filter, reporting progress.

    Progress is measured by how far into the *compressed* file the reader has
    reached, because the member list of a compressed tar is only knowable by
    decompressing the whole archive -- counting unpacked bytes instead would mean
    paying for the expensive part twice.  Members are pulled through a generator
    so the offset can be sampled as each one is reached, while ``extractall``
    still applies the extraction filter to every member.
    """
    total = Path(archive).stat().st_size
    progress = _Progress(total, log=log)
    with open(archive, "rb") as handle:
        with tarfile.open(fileobj=handle, mode=mode) as tf:
            def members():
                for member in tf:
                    progress.set_position(handle.tell())
                    yield member
            tf.extractall(dest, members=members(), filter="data")
    progress.set_position(total)
    progress.finish()


def _extract_strip_top(archive, dest, *, log=print):
    """Extract a GitHub source/datasets zip, stripping its single top-level
    directory (``galacticus-master/`` / ``datasets-master/``) so contents land
    directly in `dest`.

    The staging directory sits beside `dest` (not in the system temporary
    directory) so that lifting the contents out of it is a rename rather than a
    copy of the whole tree onto another filesystem.
    """
    dest.mkdir(parents=True, exist_ok=True)
    log(f"Unpacking {Path(archive).name} ...")
    with tempfile.TemporaryDirectory(dir=dest.parent, prefix=_STAGING) as work:
        _extract_zip(archive, Path(work), log=log)
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
