"""Unit tests for provisioning a managed install.

Network-free: `_download` is replaced by a copy out of a locally built fixture
"release", so the component selection logic, the ``--no-tools`` choice and its
persistence, the choice between an artefact's current and legacy form, and the
parameter-catalog download/generate fallback are all exercised against real
archives on disk.
"""

import io
import json
import tarfile
import zipfile
from pathlib import Path

import pytest

pytest.importorskip("platformdirs")
pytest.importorskip("requests")

from galacticus_launcher import download, paths, platforms


# The stub the fixture release ships as `scripts/build/parameterCatalog.py`, so
# that the fallback path really does run the generator as a subprocess.
_GENERATOR = """\
import json, sys
json.dump({"source": "generated"}, open(sys.argv[2], "w"))
"""


@pytest.fixture
def release(tmp_path):
    """A directory of fixture assets, keyed by the file name each URL ends in."""
    assets = tmp_path / "release"
    assets.mkdir()

    tree = tmp_path / "tree" / "galacticus-master"
    (tree / "scripts" / "build").mkdir(parents=True)
    (tree / "parameters").mkdir()
    (tree / "scripts" / "build" / "parameterCatalog.py").write_text(_GENERATOR)
    (tree / "parameters" / "quickTest.xml").write_text("<parameters/>")
    _zip(assets / "master.zip", tree.parent)

    data = tmp_path / "data" / "datasets-master" / "static"
    data.mkdir(parents=True)
    (data / "table.txt").write_text("data")
    _zip(assets / "datasets.zip", data.parent.parent)

    tools = tmp_path / "tools" / "dynamic"
    tools.mkdir(parents=True)
    (tools / "camb").write_text("#!/bin/sh\n")
    with tarfile.open(assets / "tools.tar.bz2", "w:bz2") as archive:
        archive.add(tools, arcname="dynamic")
    _tar_zst(assets / "tools.tar.zst",
             lambda archive: archive.add(tools, arcname="dynamic"))

    # The datasets snapshot is packed from the root of the datasets tree, so it
    # has no top-level directory to strip (unlike the repository archive above).
    _tar_zst(assets / "datasets.tar.zst",
             lambda archive: archive.add(data, arcname="static"))

    (assets / "Galacticus.exe").write_text("#!/bin/true\n")
    (assets / "parameters.catalog.json").write_text(
        json.dumps({"source": "published"}))
    return assets


def _tar_zst(path, add):
    """Write a zstd-compressed tar built by `add`, as the release assets are."""
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as archive:
        add(archive)
    path.write_bytes(_zstd(buffer.getvalue()))


def _zstd(data):
    try:
        from compression.zstd import compress        # Python 3.14 and later
    except ImportError:
        import zstandard
        return zstandard.ZstdCompressor().compress(data)
    return compress(data)


def _zip(path, root):
    with zipfile.ZipFile(path, "w") as archive:
        for item in sorted(Path(root).rglob("*")):
            archive.write(item, item.relative_to(root))


@pytest.fixture
def install(tmp_path):
    root = tmp_path / "install"
    return paths.Install(
        source=paths.SOURCE_MANAGED, tag="v1.2.3",
        exec_path=root / "exec", data_path=root / "datasets",
        tools_path=root / "tools", dynamic_path=tmp_path / "cache" / "dynamic",
        binary=root / "exec" / "Galacticus.exe",
        assets=platforms.PlatformAssets("Galacticus.exe",
                                        "tools.tar.zst", "tar.zst",
                                        "tools.tar.bz2", "tar.bz2", "test"),
    )


@pytest.fixture
def fetched(monkeypatch, release):
    """Serve `_download` from the fixture release; record every URL requested."""
    requested = []

    def _fake_download(url, dest, *, log=print, retries=4, missing_ok=False,
                       progress=None):
        requested.append(url)
        # Both repository archives are `.../archive/...zip`, so key on the repo;
        # every other URL is named for the release asset it carries.
        if f"/{download.DATASETS_REPO}/archive/" in url:
            name = "datasets.zip"
        elif "/archive/" in url:
            name = "master.zip"
        else:
            name = Path(url).name
        source = release / name
        if not source.is_file():
            if missing_ok:
                return False
            raise RuntimeError(f"no fixture asset {name}")
        Path(dest).write_bytes(source.read_bytes())
        return True

    monkeypatch.setattr(download, "_download", _fake_download)
    monkeypatch.setattr(download, "load_checksums", lambda tag, log=print: None)
    monkeypatch.setattr(download, "resolve_datasets_ref",
                        lambda tag, log=print: "master")
    # Nothing may reach the network: the small text assets (the datasets pin)
    # are read through this, and a release which publishes none is the default.
    monkeypatch.setattr(download, "_read_remote_text", lambda url: None)
    return requested


def _catalog(install):
    return json.loads((install.exec_path / "parameters.catalog.json").read_text())


def _quiet(message):
    pass


# --- the full set of components -------------------------------------------

def test_provision_fetches_every_component(install, fetched):
    done = download.provision(install, log=_quiet)
    assert done == ["exec", "datasets", "tools", "parameter catalog"]
    assert (install.exec_path / "parameters" / "quickTest.xml").is_file()
    assert (install.data_path / "static" / "table.txt").is_file()
    # Tools are lifted out of the archive's `dynamic/` prefix.
    assert (install.tools_path / "camb").is_file()
    # The published catalog is preferred over generating one locally.
    assert _catalog(install) == {"source": "published"}


def test_provision_is_idempotent(install, fetched):
    download.provision(install, log=_quiet)
    assert download.provision(install, log=_quiet) == []


# --- skipping the tools archive -------------------------------------------

def test_no_tools_skips_the_archive_and_records_the_choice(install, fetched):
    done = download.provision(install, tools=False, log=_quiet)
    assert "tools" not in done
    assert not any("tools.tar.bz2" in url for url in fetched)
    assert not (install.tools_path / "camb").exists()
    assert download.tools_skipped(install) and not download.tools_installed(install)


def test_recorded_skip_survives_a_later_provision(install, fetched):
    download.provision(install, tools=False, log=_quiet)
    fetched.clear()
    # A bare `run`/`validate` provisions on demand and must honor the choice --
    # otherwise the next command downloads the archive just declined.
    assert download.provision(install, log=_quiet) == []
    assert not any("tools.tar.bz2" in url for url in fetched)


def test_asking_for_tools_clears_the_skip(install, fetched):
    download.provision(install, tools=False, log=_quiet)
    done = download.provision(install, tools=True, log=_quiet)
    assert done == ["tools"]
    assert (install.tools_path / "camb").is_file()
    assert download.tools_installed(install) and not download.tools_skipped(install)


def test_no_tools_never_removes_installed_tools(install, fetched):
    download.provision(install, log=_quiet)
    messages = []
    assert download.provision(install, tools=False, log=messages.append) == []
    assert (install.tools_path / "camb").is_file()
    assert not download.tools_skipped(install)
    assert any("already installed" in message for message in messages)


# --- the parameter catalog ------------------------------------------------

def test_catalog_is_generated_when_the_release_publishes_none(install, fetched,
                                                              release):
    (release / "parameters.catalog.json").unlink()
    done = download.provision(install, log=_quiet)
    assert "parameter catalog" in done
    assert _catalog(install) == {"source": "generated"}


def test_bad_catalog_checksum_falls_back_to_generating_one(install, fetched,
                                                           monkeypatch):
    monkeypatch.setattr(download, "load_checksums",
                        lambda tag, log=print: {"parameters.catalog.json": "0" * 64})
    # A mismatch is fatal for artefacts which get executed; the catalog is data,
    # and a local build of it is equivalent, so the install continues.
    done = download.provision(install, log=_quiet)
    assert "parameter catalog" in done
    assert _catalog(install) == {"source": "generated"}


def test_malformed_catalog_falls_back_to_generating_one(install, fetched, release):
    (release / "parameters.catalog.json").write_text("{ not json")
    download.provision(install, log=_quiet)
    assert _catalog(install) == {"source": "generated"}


# --- choosing between an artefact's current and legacy form ----------------

def _checksums(*assets):
    """A release asset listing, with digests no test verifies against."""
    return {asset: "" for asset in assets}


def test_tools_come_from_the_zstd_archive_when_the_release_publishes_one(
        install, fetched, monkeypatch):
    monkeypatch.setattr(download, "load_checksums",
                        lambda tag, log=print: _checksums("tools.tar.zst"))
    monkeypatch.setattr(download, "_verify", lambda *a, **k: None)
    assert "tools" in download.provision(install, log=_quiet)
    assert any(url.endswith("tools.tar.zst") for url in fetched)
    assert not any(url.endswith("tools.tar.bz2") for url in fetched)
    assert (install.tools_path / "camb").is_file()


def test_tools_fall_back_to_the_archive_an_older_release_published(
        install, fetched, monkeypatch):
    """A release cut before the switch to zstd carries only the older archive,
    and it can not be regenerated, so its name has to remain reachable."""
    monkeypatch.setattr(download, "load_checksums",
                        lambda tag, log=print: _checksums("tools.tar.bz2"))
    monkeypatch.setattr(download, "_verify", lambda *a, **k: None)
    assert "tools" in download.provision(install, log=_quiet)
    assert any(url.endswith("tools.tar.bz2") for url in fetched)
    assert not any(url.endswith("tools.tar.zst") for url in fetched)
    assert (install.tools_path / "camb").is_file()


# --- the datasets snapshot -------------------------------------------------

def test_datasets_come_from_the_release_snapshot_when_published(
        install, fetched, monkeypatch):
    monkeypatch.setattr(download, "load_checksums",
                        lambda tag, log=print: _checksums("datasets.tar.zst"))
    monkeypatch.setattr(download, "_verify", lambda *a, **k: None)
    monkeypatch.setattr(download, "_snapshot_ref", lambda tag: "abc123")
    assert "datasets" in download.provision(install, log=_quiet)
    assert any(url.endswith("datasets.tar.zst") for url in fetched)
    assert not any(f"/{download.DATASETS_REPO}/archive/" in url for url in fetched)
    assert (install.data_path / "static" / "table.txt").is_file()
    # The sentinel records the commit the snapshot was taken from, so that
    # `galacticus info` and a later `update` can tell which data is installed.
    assert (install.data_path / ".galacticus-datasets").read_text() == "abc123"


def test_datasets_fall_back_to_the_repository_archive(install, fetched,
                                                      monkeypatch):
    monkeypatch.setattr(download, "load_checksums",
                        lambda tag, log=print: _checksums("Galacticus.exe"))
    monkeypatch.setattr(download, "_verify", lambda *a, **k: None)
    assert "datasets" in download.provision(install, log=_quiet)
    assert any(f"/{download.DATASETS_REPO}/archive/" in url for url in fetched)
    assert (install.data_path / "static" / "table.txt").is_file()


def test_an_explicit_datasets_ref_bypasses_the_snapshot(install, fetched,
                                                        monkeypatch):
    """`GALACTICUS_DATASETS_REF` asks for a specific commit, which the release's
    own snapshot can not satisfy -- so that request must reach the repository."""
    monkeypatch.setenv("GALACTICUS_DATASETS_REF", "somebranch")
    monkeypatch.setattr(download, "load_checksums",
                        lambda tag, log=print: _checksums("datasets.tar.zst"))
    monkeypatch.setattr(download, "_verify", lambda *a, **k: None)
    download.provision(install, log=_quiet)
    assert any(f"/{download.DATASETS_REPO}/archive/" in url for url in fetched)
    assert not any(url.endswith("datasets.tar.zst") for url in fetched)


def test_a_failed_catalog_download_does_not_abort_the_install(install, fetched,
                                                              monkeypatch):
    """The catalog is the one artefact provisioning can rebuild itself, so a
    transfer failure has to fall back to generating it -- not discard the
    components which already downloaded."""
    real = download._download

    def failing(url, dest, **kwargs):
        if url.endswith(download.CATALOG_ASSET):
            raise RuntimeError("connection reset by peer")
        return real(url, dest, **kwargs)

    monkeypatch.setattr(download, "_download", failing)
    done = download.provision(install, log=_quiet)
    assert done == ["exec", "datasets", "tools", "parameter catalog"]
    assert (install.exec_path / "parameters" / "quickTest.xml").is_file()
    assert (install.tools_path / "camb").is_file()
    assert _catalog(install) == {"source": "generated"}
