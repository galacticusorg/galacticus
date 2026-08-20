"""Fetch and unpack Galacticus release artefacts into a managed install.

Provisioning is idempotent: each component drops a sentinel file once complete,
and a present sentinel means the component is skipped.  Only managed installs
are ever provisioned -- a local build/clone is used in place, untouched.

Four components make up a runnable install:

* **exec**     -- the repository source archive (whole tree), giving
  ``GALACTICUS_EXEC_PATH`` its ``parameters/``, ``aux/`` and ``scripts/``,
  plus the platform executable copied in as ``Galacticus.exe``.
* **datasets** -- the run-time static data root.  Releases publish a snapshot of
  ``galacticusorg/datasets`` as an asset; older ones are fetched from the
  repository archive instead (see :func:`_plan_datasets`).
* **tools**    -- the pre-built run-time tools archive.  Its entries are
  prefixed ``dynamic/`` (the location the un-relocated binary expects); we strip
  that prefix so the contents land directly under ``GALACTICUS_TOOLS_PATH``.
  This is by far the largest artefact, and models which use none of the tools it
  carries (CAMB, CLASS, Cloudy, ...) do not need it, so it is the one component
  which can be skipped -- see :func:`provision`.
* **parameter catalog** -- the machine-readable catalog of accepted input
  parameters used by ``galacticus validate``, downloaded from the release when it
  publishes one and generated from the installed source otherwise.

A first install moves gigabytes, so the transfer is arranged to keep the link
busy rather than to be simple:

* Every component's archive is downloaded **concurrently** with the others.  The
  repository archive endpoint throttles hard per connection, so overlapping the
  components is worth far more than it costs -- separate connections each get
  their own share.
* A single large asset is split across several **byte-range** requests where the
  server supports them (release assets do; the on-the-fly repository archives do
  not), which lifts a single-stream transfer to roughly the link rate.
* Archives are staged alongside their destination rather than in the system
  temporary directory, so that unpacking is followed by a rename rather than a
  copy across filesystems, and so that a multi-gigabyte download does not have to
  fit in a (often small, often ``tmpfs``) ``/tmp``.
"""

import contextlib
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import threading
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

REPO = "galacticusorg/galacticus"
DATASETS_REPO = "galacticusorg/datasets"

# Release asset listing the SHA-256 of each of the other assets, used to verify what we download.
CHECKSUM_ASSET = "SHA256SUMS"

# Release asset carrying the pre-built parameter catalog (see `_plan_catalog`).
CATALOG_ASSET = "parameters.catalog.json"

# Release assets carrying the datasets snapshot and the commit it was taken from.
DATASETS_ASSET = "datasets.tar.zst"
DATASETS_REF_ASSET = "datasets.ref"

# Sentinel recording that the tools archive is deliberately not wanted.
TOOLS_SKIPPED = "tools-skipped"

# Tool executables in the pre-built tools archive whose execute bit must be
# restored after unpacking a zip.  Tar archives -- the current format on every
# platform -- record the mode themselves, so this only matters when falling back
# to the macOS `.zip` tools archive published by older releases.  This list
# matches the binaries the CI `Build-Tools` jobs pack; tools built at run time
# and NOT shipped (CosmicEmu's emu.exe, AGN_Spectrum) are intentionally absent.
_TOOL_EXECUTABLE_NAMES = frozenset({
    "camb", "class", "recfast.exe", "cloudy.exe", "autosps.exe", "harmonize",
})

_CHUNK = 1 << 20  # 1 MiB

# Prefix for the staging directories archives are downloaded into and unpacked
# in. They sit beside their destination rather than in the system temporary
# directory (see the module docstring), so name them recognizably: an install
# killed part way through leaves one behind, and it should be obvious what it is.
_STAGING = ".galacticus-staging-"

# Byte-range requests to split one asset across.  Four is where the measured
# gain flattens out; `GALACTICUS_DOWNLOAD_CONNECTIONS` overrides it, and 1
# disables splitting altogether (for a link or proxy which dislikes it).
_CONNECTIONS = 4

# Ceiling on connections in flight at once across every concurrent download, so
# that splitting several components at the same time stays polite.
_MAX_CONNECTIONS = 8

# Assets smaller than this are fetched in one request: splitting them costs more
# in round trips than it saves in transfer time.
_SPLIT_THRESHOLD = 64 << 20  # 64 MiB

_CONNECTION_SLOTS = threading.Semaphore(_MAX_CONNECTIONS)


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

    ``GALACTICUS_DATASETS_REF`` overrides everything. Otherwise a release pins
    datasets to the commit recorded in its ``datasets.ref`` asset (written by the
    workflow which publishes the snapshot), so a given release always installs
    the same data. A release with no pin tracks datasets ``master``.
    """
    override = os.environ.get("GALACTICUS_DATASETS_REF")
    if override:
        return override
    pinned = _read_remote_text(asset_url(tag, DATASETS_REF_ASSET))
    if pinned:
        return pinned.strip().split()[0]
    if tag and tag != "bleeding-edge":
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

    It doubles as the manifest of what a release carries: an asset listed here is
    one the release has, which is how provisioning chooses between the current
    and the legacy form of an artefact without probing for a 404.
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


def _publishes(checksums, asset):
    """True if release `checksums` list `asset` -- i.e. the release carries it."""
    return checksums is not None and asset in checksums


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


class _Plan:
    """What one component needs: files to download, then how to install them.

    Provisioning is split in two so that the slow part can be overlapped.  Each
    component's planner registers its downloads with the shared :class:`_Fetcher`
    and returns the callable which unpacks them; every download then runs
    concurrently, and the unpack callables run afterwards, in component order.
    """

    def __init__(self, name, unpack):
        self.name = name
        self.unpack = unpack


class _Context:
    """Everything a component planner needs, gathered once per `provision` call."""

    def __init__(self, *, install, checksums, force, log, fetcher, stack):
        self.install = install
        self.checksums = checksums
        self.force = force
        self.log = log
        self.fetcher = fetcher
        self.stack = stack

    def workspace(self, beside):
        """A staging directory next to `beside`, removed when provisioning ends."""
        beside.mkdir(parents=True, exist_ok=True)
        return Path(self.stack.enter_context(
            tempfile.TemporaryDirectory(dir=beside.parent, prefix=_STAGING)))


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

    for path in (install.exec_path, install.data_path, install.tools_path,
                 install.dynamic_path):
        path.mkdir(parents=True, exist_ok=True)

    # Fetch the published checksums once, and verify each downloaded artefact against them below.
    checksums = load_checksums(install.tag, log=log)
    if checksums is None:
        log(f"  note: release {install.tag} publishes no {CHECKSUM_ASSET}; "
            "downloaded artefacts can not be verified against it.")

    done = []
    fetcher = _Fetcher(log=log)
    with contextlib.ExitStack() as stack:
        context = _Context(install=install, checksums=checksums, force=force,
                           log=log, fetcher=fetcher, stack=stack)
        plans = [plan for plan in (_plan_exec(context),
                                   _plan_datasets(context),
                                   _plan_tools(context, tools),
                                   _plan_catalog(context))
                 if plan is not None]
        if not plans:
            return done
        fetcher.run()
        for plan in plans:
            if plan.unpack() is not False:
                done.append(plan.name)
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


def _plan_catalog(context):
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
    install, log = context.install, context.log
    catalog = Path(install.exec_path) / CATALOG_ASSET
    if catalog.is_file() and not context.force:
        return None
    staged = Path(str(catalog) + ".part-download")
    job = None
    if context.checksums is None or CATALOG_ASSET in context.checksums:
        log("Fetching parameter catalog ...")
        # Optional: the catalog can be rebuilt from the source tree, so neither
        # a release which does not publish one nor a failure to fetch it may
        # abort an install whose other components succeeded.
        job = context.fetcher.enqueue(asset_url(install.tag, CATALOG_ASSET), staged,
                                      label=CATALOG_ASSET, missing_ok=True,
                                      optional=True)
    else:
        log(f"  release {install.tag} publishes no {CATALOG_ASSET}.")

    def unpack():
        if job is not None and job.present and _accept_catalog(
                staged, context.checksums, log=log):
            staged.replace(catalog)
            return True
        staged.unlink(missing_ok=True)
        return _generate_catalog(install, catalog, log=log)

    return _Plan("parameter catalog", unpack)


def _accept_catalog(staged, checksums, *, log):
    """True if the downloaded catalog at `staged` is fit to install.

    Every failure mode -- a checksum mismatch, a file which is not valid JSON --
    returns False so the caller can generate one instead.  Unlike the executable
    and the tools, the catalog is data and is never run, so a bad copy is a reason
    to fall back rather than to abort the install; the generator reproduces the
    same catalog from the source tree just installed.
    """
    try:
        _verify(staged, CATALOG_ASSET, checksums, log=log)
        with open(staged) as handle:
            json.load(handle)
    except (RuntimeError, OSError, ValueError) as error:
        log(f"  could not use the published parameter catalog ({error}).")
        return False
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


def _plan_exec(context):
    install, log = context.install, context.log
    sentinel = _sentinel(install.exec_path, "exec")
    if sentinel.exists() and not context.force:
        return None
    log(f"Fetching Galacticus source ({install.tag}) ...")
    work = context.workspace(install.exec_path)
    archive = work / "source.zip"
    # The executable is staged beside the archive rather than downloaded straight
    # into `exec_path`, because unpacking the source tree into that directory
    # happens afterwards and replaces whatever it finds under each of its own
    # top-level names.
    staged_binary = work / "Galacticus.exe"
    context.fetcher.enqueue(source_url(install.tag), archive, label="source.zip")
    context.fetcher.enqueue(asset_url(install.tag, install.assets.binary),
                            staged_binary, label=install.assets.binary)

    def unpack():
        _extract_strip_top(archive, install.exec_path, log=log)
        # Verify before making the file executable, so that a file which fails the check is never left runnable.
        _verify(staged_binary, install.assets.binary, context.checksums, log=log)
        binary = install.exec_path / "Galacticus.exe"
        binary.unlink(missing_ok=True)
        shutil.move(str(staged_binary), str(binary))
        binary.chmod(0o755)
        sentinel.write_text(install.tag)
        return True

    return _Plan("exec", unpack)


def _plan_datasets(context):
    """Install the static data root, preferring the snapshot the release publishes.

    The datasets repository archive is generated on the fly by GitHub and served
    from a throttled endpoint which supports neither byte ranges nor a content
    length, so it is the slowest artefact of an install by a wide margin.  A
    release which publishes a `datasets.tar.zst` snapshot is used instead: it is a
    release asset like any other, so it can be split across connections, verified
    against the release checksums, and unpacked far faster than a zip.  It also
    pins the data to the commit the release was built against rather than to
    whatever ``master`` happens to be, which is why the sentinel records the
    snapshot's own ref.  Setting ``GALACTICUS_DATASETS_REF`` asks for a specific
    ref and so always takes the repository-archive path.
    """
    install, log = context.install, context.log
    sentinel = _sentinel(install.data_path, "datasets")
    if sentinel.exists() and not context.force:
        return None
    override = os.environ.get("GALACTICUS_DATASETS_REF")
    if not override and _publishes(context.checksums, DATASETS_ASSET):
        log(f"Fetching datasets snapshot ({install.tag}) ...")
        work = context.workspace(install.data_path)
        archive = work / DATASETS_ASSET
        context.fetcher.enqueue(asset_url(install.tag, DATASETS_ASSET), archive,
                                label=DATASETS_ASSET)

        def unpack():
            _verify(archive, DATASETS_ASSET, context.checksums, log=log)
            log(f"Unpacking {DATASETS_ASSET} ...")
            # The snapshot is packed from the root of the datasets tree, so there
            # is no top-level directory to strip -- but it is still unpacked to
            # one side and moved in, so that re-provisioning replaces each
            # top-level entry outright instead of merging into a stale tree.
            staging = work / "unpacked"
            _extract(archive, staging, "tar.zst", log=log)
            _replace_children(staging, install.data_path)
            sentinel.write_text(_snapshot_ref(install.tag) or install.tag)
            return True

        return _Plan("datasets", unpack)

    ref = resolve_datasets_ref(install.tag, log=log)
    log(f"Fetching datasets ({ref}) ...")
    work = context.workspace(install.data_path)
    archive = work / "datasets.zip"
    context.fetcher.enqueue(datasets_url(ref), archive, label="datasets.zip")

    def unpack():
        _extract_strip_top(archive, install.data_path, log=log)
        sentinel.write_text(ref)
        return True

    return _Plan("datasets", unpack)


def _snapshot_ref(tag):
    """The datasets commit a release's snapshot was taken from, or None."""
    pinned = _read_remote_text(asset_url(tag, DATASETS_REF_ASSET))
    return pinned.strip().split()[0] if pinned else None


def _tools_archive(install, checksums):
    """Return the (asset name, format) of the tools archive to fetch.

    Releases cut before the switch to zstd carry only the older archive, and it
    can not be regenerated for them, so the choice is made from what the release
    actually publishes rather than from the platform alone.
    """
    assets = install.assets
    legacy = (assets.tools_legacy, assets.tools_legacy_format)
    current = (assets.tools, assets.tools_format)
    if checksums is None:
        # Old enough to publish no checksums at all is old enough to predate the
        # current archive; try the legacy name, which is all such a release has.
        return legacy if assets.tools_legacy else current
    if assets.tools in checksums:
        return current
    if assets.tools_legacy in checksums:
        return legacy
    return current      # published neither: ask for the current one and fail loudly


def _plan_tools(context, tools):
    install, log = context.install, context.log
    if not _tools_wanted(install, tools, log=log):
        return None
    sentinel = _sentinel(install.tools_path, "tools")
    if sentinel.exists() and not context.force:
        return None
    name, fmt = _tools_archive(install, context.checksums)
    log(f"Fetching tools {name} ...")
    work = context.workspace(install.tools_path)
    archive = work / name
    context.fetcher.enqueue(asset_url(install.tag, name), archive, label=name)

    def unpack():
        # Verify before unpacking: this archive contains executables which are later run.
        _verify(archive, name, context.checksums, log=log)
        staging = work / "unpacked"
        log(f"Unpacking {name} ...")
        _extract(archive, staging, fmt, log=log)
        # Tools archives are rooted at "dynamic/"; lift that subtree up so the
        # contents sit directly under GALACTICUS_TOOLS_PATH.
        root = staging / "dynamic"
        _replace_children(root if root.is_dir() else staging, install.tools_path)
        if fmt == "zip":
            _restore_executable_bits(install.tools_path)
        sentinel.write_text(install.tag)
        return True

    return _Plan("tools", unpack)


class _Job:
    """One queued download, and whether it turned out to be there.

    `optional` marks an artefact the install can do without.  Such a job records
    a failure rather than raising it, because the components are fetched together
    and a failure escaping the fetch phase discards every staging directory --
    losing gigabytes which did arrive over one small artefact which did not.
    """

    def __init__(self, url, dest, label, missing_ok, optional):
        self.url = url
        self.dest = dest
        self.label = label
        self.missing_ok = missing_ok
        self.optional = optional
        self.present = False
        self.error = None


class _Fetcher:
    """Collect the downloads a provisioning run needs, then run them together.

    The components are independent, and the repository archive endpoint throttles
    each connection rather than the account, so fetching them at the same time is
    close to free: the slowest component sets the wall clock instead of the sum.
    """

    def __init__(self, *, log=print):
        self._jobs = []
        self._log = log

    def enqueue(self, url, dest, *, label, missing_ok=False, optional=False):
        job = _Job(url, dest, label, missing_ok, optional)
        self._jobs.append(job)
        return job

    def run(self):
        if not self._jobs:
            return
        with _display(log=self._log) as display:
            if len(self._jobs) == 1:
                self._fetch(self._jobs[0], display)
                return
            with ThreadPoolExecutor(len(self._jobs)) as pool:
                # `map` is consumed eagerly so that the first failure is raised
                # here, once every job has been given its chance to finish.
                list(pool.map(lambda job: self._fetch(job, display), self._jobs))

    def _fetch(self, job, display):
        progress = display.add(job.label)
        try:
            job.present = _download(job.url, job.dest, log=self._log,
                                    missing_ok=job.missing_ok, progress=progress)
        except (RuntimeError, OSError) as error:
            if not job.optional:
                raise
            job.error = error
            self._log(f"  could not fetch {job.label} ({error}).")
        finally:
            progress.finish()


def _connections():
    """How many byte-range requests to split one large asset across."""
    try:
        value = int(os.environ.get("GALACTICUS_DOWNLOAD_CONNECTIONS", _CONNECTIONS))
    except ValueError:
        return _CONNECTIONS
    return max(1, min(value, _MAX_CONNECTIONS))


def _download(url, dest, *, log=print, retries=4, missing_ok=False, progress=None):
    """Fetch `url` to `dest`, splitting it across connections where that helps.

    Returns True once `dest` is written.  With `missing_ok`, a 404 -- the release
    simply does not publish this asset -- returns False immediately instead of
    retrying an absence through the full backoff schedule.

    The first request asks for ``bytes=0-``: a server which honours byte ranges
    answers 206 and states the total size, which is the signal that the transfer
    can be split (and, for a large enough asset, it then is).  One which does not
    -- GitHub's on-the-fly repository archives -- answers 200, and its response is
    the transfer, so nothing is wasted probing.
    """
    own_display = None
    if progress is None:
        own_display = _Display(log=log)
        progress = own_display.add(Path(dest).name)
    try:
        return _fetch(url, dest, progress=progress, log=log, retries=retries,
                      missing_ok=missing_ok)
    finally:
        if own_display is not None:
            progress.finish()
            own_display.close()


def _fetch(url, dest, *, progress, log, retries, missing_ok):
    last_error = None
    for attempt in range(retries):
        try:
            with _CONNECTION_SLOTS:
                response = requests.get(url, stream=True, timeout=60,
                                        headers={"Range": "bytes=0-"})
                with response:
                    if missing_ok and response.status_code == 404:
                        return False
                    response.raise_for_status()
                    total = _total_length(response)
                    progress.set_total(total)
                    splittable = (response.status_code == 206
                                  and total is not None
                                  and total >= _SPLIT_THRESHOLD
                                  and _connections() > 1)
                    if not splittable:
                        _stream_to(response, dest, progress)
                        return True
            # Splittable: the probe response is closed (and its connection slot
            # released) before the range workers ask for slots of their own.
            _fetch_split(url, dest, total, progress=progress, retries=retries)
            return True
        except (requests.RequestException, OSError) as error:  # pragma: no cover - network
            last_error = error
            if attempt == retries - 1:
                break
            wait = 2 ** (attempt + 1)
            log(f"  download failed ({error}); retrying in {wait}s ...")
            progress.reset()
            time.sleep(wait)
    raise RuntimeError(f"failed to download {url}: {last_error}")


def _stream_to(response, dest, progress):
    """Write the body of `response` to `dest` through a `.part` file."""
    tmp = Path(str(dest) + ".part")
    with open(tmp, "wb") as handle:
        for chunk in response.iter_content(chunk_size=_CHUNK):
            if chunk:
                handle.write(chunk)
                progress.update(len(chunk))
    tmp.replace(dest)


def _fetch_split(url, dest, total, *, progress, retries):
    """Fetch `url` as several byte ranges written into one preallocated file.

    Each worker owns a disjoint span and writes it at its own offset, so no
    reassembly step is needed.  A worker which is cut off resumes from where it
    stopped rather than restarting its span, which keeps both the transfer and
    the progress total honest across a retry.
    """
    connections = _connections()
    tmp = Path(str(dest) + ".part")
    with open(tmp, "wb") as handle:
        handle.truncate(total)
    span = -(-total // connections)         # ceiling division
    spans = [(start, min(start + span, total) - 1)
             for start in range(0, total, span)]
    with ThreadPoolExecutor(len(spans)) as pool:
        list(pool.map(lambda bounds: _fetch_span(url, tmp, *bounds,
                                                 progress=progress, retries=retries),
                      spans))
    tmp.replace(dest)


def _fetch_span(url, path, start, end, *, progress, retries):
    position = start
    last_error = None
    for attempt in range(retries):
        try:
            with _CONNECTION_SLOTS:
                headers = {"Range": f"bytes={position}-{end}"}
                with requests.get(url, stream=True, timeout=60, headers=headers) as response:
                    response.raise_for_status()
                    with open(path, "r+b") as handle:
                        handle.seek(position)
                        for chunk in response.iter_content(chunk_size=_CHUNK):
                            if chunk:
                                handle.write(chunk)
                                position += len(chunk)
                                progress.update(len(chunk))
            if position > end:
                return
            last_error = OSError(f"stream ended at byte {position}, short of {end}")
        except (requests.RequestException, OSError) as error:  # pragma: no cover - network
            last_error = error
        if attempt < retries - 1:
            time.sleep(2 ** (attempt + 1))
    raise RuntimeError(f"failed to download bytes {start}-{end} of {url}: {last_error}")


def _total_length(response):
    """Total size of the resource in bytes, or None if it is not advertised.

    A 206 states the total after the ``/`` of ``Content-Range``; a 200 states it
    in ``Content-Length`` (which on a 206 would be the length of the range only).
    """
    content_range = response.headers.get("Content-Range")
    if content_range and "/" in content_range:
        total = content_range.rsplit("/", 1)[1].strip()
        if total.isdigit():
            return int(total)
    value = response.headers.get("Content-Length")
    try:
        return int(value) if value is not None else None
    except (TypeError, ValueError):  # pragma: no cover - malformed header
        return None


class _Progress:
    """One labelled bar's worth of state; rendered by its :class:`_Display`."""

    def __init__(self, label, total, *, display):
        self.label = label
        self._total = total if total and total > 0 else None
        self._display = display
        self._done = 0
        self._next_milestone = _Display.MILESTONE

    def set_total(self, total):
        """Record the size once the response headers state it."""
        with self._display.lock:
            self._total = total if total and total > 0 else None

    def reset(self):
        """Start counting again, after a failed attempt is retried."""
        with self._display.lock:
            self._done = 0
            self._next_milestone = _Display.MILESTONE

    def update(self, count):
        with self._display.lock:
            self._done += count
            milestone = self._milestone()
        self._display.touch(self, milestone)

    def set_position(self, position):
        """Advance to an absolute `position`, in bytes; never moves backwards."""
        with self._display.lock:
            if position <= self._done:
                return
            self._done = position
            milestone = self._milestone()
        self._display.touch(self, milestone)

    def finish(self):
        self._display.touch(self, None, force=True)

    def _milestone(self):
        """The percentage to report as a line, or None; caller holds the lock."""
        if self._total is None:
            return None
        fraction = self._done / self._total
        if fraction < self._next_milestone:
            return None
        while self._next_milestone <= fraction:
            self._next_milestone += _Display.MILESTONE
        return f"{int(fraction * 100)}% ({_human(self._done)} / {_human(self._total)})"

    def render(self):
        """The bar as a single line of text; caller holds the lock."""
        if self._total is None:
            return f"  {self.label}: {_human(self._done)}"
        fraction = min(1.0, self._done / self._total)
        filled = int(_Display.BAR_WIDTH * fraction)
        bar = "#" * filled + "-" * (_Display.BAR_WIDTH - filled)
        return (f"  [{bar}] {int(fraction * 100):3d}% "
                f"{_human(self._done)} / {_human(self._total)}  {self.label}")


class _Display:
    """Render progress for one or more transfers at once.

    On an interactive terminal each item gets a line, and the whole block is
    redrawn in place; otherwise progress is emitted as occasional milestone lines
    tagged with the item's label, so a log or a pipe stays readable.  Dependency
    free -- no ``tqdm`` -- to keep the launcher's footprint minimal.
    """

    BAR_WIDTH = 30
    MILESTONE = 0.10  # log every 10% when not a TTY
    INTERVAL = 0.1    # seconds between redraws

    def __init__(self, *, log=print):
        self._log = log
        self._items = []
        self._drawn = 0
        self._last_render = 0.0
        self.lock = threading.RLock()
        self._tty = (
            log is print
            and hasattr(sys.stderr, "isatty")
            and sys.stderr.isatty()
        )

    def add(self, label, total=None):
        item = _Progress(label, total, display=self)
        with self.lock:
            self._items.append(item)
        return item

    def touch(self, item, milestone, *, force=False):
        """Note that `item` advanced; redraw, or log its `milestone` if given."""
        if self._tty:
            self._render(force=force)
        elif milestone is not None:
            self._log(f"  ... {item.label} {milestone}")

    def _render(self, *, force=False):
        with self.lock:
            now = time.monotonic()
            if not force and now - self._last_render < self.INTERVAL:
                return
            self._last_render = now
            text = ""
            if self._drawn:
                text += f"\033[{self._drawn}A"    # back to the top of the block
            for item in self._items:
                text += f"\r\033[K{item.render()}\n"
            self._drawn = len(self._items)
            sys.stderr.write(text)
            sys.stderr.flush()

    def close(self):
        """Leave the final state on screen and stop drawing over it."""
        with self.lock:
            if self._tty and self._items:
                self._render(force=True)
            self._items = []
            self._drawn = 0


@contextlib.contextmanager
def _display(*, log=print):
    display = _Display(log=log)
    try:
        yield display
    finally:
        display.close()


@contextlib.contextmanager
def _bar(label, total, *, log=print):
    """A single progress bar, closed when the block ends."""
    with _display(log=log) as display:
        progress = display.add(label, total)
        yield progress
        progress.finish()


def _human(num_bytes):
    value = float(num_bytes)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024.0 or unit == "TiB":
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024.0


def _extract(archive, dest, fmt, *, log=print):
    """Unpack `archive` (`fmt` in {'zip','tar.zst','tar.bz2','tar.gz'}) into `dest`.

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
        _extract_tar(archive, dest, fmt, log=log)


def _extract_zip(archive, dest, *, log=print):
    """Unpack a zip member by member, reporting progress by unpacked bytes.

    Equivalent to ``ZipFile.extractall`` -- which is itself this loop -- other
    than the progress reporting; member paths are sanitized by ``zipfile`` in
    either case.  The member list comes from the central directory, so the total
    is known up front without reading the archive.
    """
    with zipfile.ZipFile(archive) as zf:
        members = zf.infolist()
        total = sum(member.file_size for member in members)
        with _bar(Path(archive).name, total, log=log) as progress:
            for member in members:
                zf.extract(member, dest)
                progress.update(member.file_size)


def _zstd_stream(handle):
    """A file-like object decompressing zstd from the open file `handle`."""
    try:
        from compression.zstd import ZstdFile        # Python 3.14 and later
    except ImportError:
        pass
    else:
        return ZstdFile(handle, "rb")
    try:
        import zstandard
    except ImportError as error:                     # pragma: no cover - packaging
        raise RuntimeError(
            "unpacking this archive needs zstd support, which this interpreter "
            "does not have: install it with `pip install zstandard`, or use "
            "Python 3.14 or later."
        ) from error
    return zstandard.ZstdDecompressor().stream_reader(handle)


def _extract_tar(archive, dest, fmt, *, log=print):
    """Unpack a compressed tar with the ``data`` filter, reporting progress.

    Progress is measured by how far into the *compressed* file the reader has
    reached, because the member list of a compressed tar is only knowable by
    decompressing the whole archive -- counting unpacked bytes instead would mean
    paying for the expensive part twice.  Members are pulled through a generator
    so the offset can be sampled as each one is reached, while ``extractall``
    still applies the extraction filter to every member.

    zstd archives are read as a stream (``r|``) layered on a decompressor, since
    the standard library's tar reader has no zstd mode of its own; bzip2 and gzip
    it opens directly.
    """
    total = Path(archive).stat().st_size
    with _bar(Path(archive).name, total, log=log) as progress:
        with open(archive, "rb") as handle:
            with contextlib.ExitStack() as stack:
                if fmt == "tar.zst":
                    stream = stack.enter_context(contextlib.closing(_zstd_stream(handle)))
                    tf = stack.enter_context(tarfile.open(fileobj=stream, mode="r|"))
                else:
                    mode = "r:bz2" if fmt == "tar.bz2" else "r:gz"
                    tf = stack.enter_context(tarfile.open(fileobj=handle, mode=mode))

                def members():
                    for member in tf:
                        progress.set_position(handle.tell())
                        yield member

                tf.extractall(dest, members=members(), filter="data")
        progress.set_position(total)


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
        _replace_children(top, dest)


def _replace_children(source, dest):
    """Move each entry of `source` into `dest`, replacing any of the same name.

    The staging directory always sits on the same filesystem as `dest`, so each
    move is a rename rather than a copy of the subtree.
    """
    for child in Path(source).iterdir():
        destination = Path(dest) / child.name
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
