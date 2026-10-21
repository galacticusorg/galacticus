"""Unit tests for provisioning a managed install.

Network-free: `_download` is replaced by a copy out of a locally built fixture
"release", so the component selection logic, the ``--no-tools`` choice and its
persistence, the choice between an artifact's current and legacy form, and the
parameter-catalog download/generate fallback are all exercised against real
archives on disk.
"""

import hashlib
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
    # A mismatch is fatal for artifacts which get executed; the catalog is data,
    # and a local build of it is equivalent, so the install continues.
    done = download.provision(install, log=_quiet)
    assert "parameter catalog" in done
    assert _catalog(install) == {"source": "generated"}


def test_malformed_catalog_falls_back_to_generating_one(install, fetched, release):
    (release / "parameters.catalog.json").write_text("{ not json")
    download.provision(install, log=_quiet)
    assert _catalog(install) == {"source": "generated"}


# --- choosing between an artifact's current and legacy form ----------------

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
    """The catalog is the one artifact provisioning can rebuild itself, so a
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


# --- a platform whose builds have been retired -----------------------------

@pytest.fixture
def retired_install(install):
    """The same install, for a platform Galacticus no longer builds."""
    return install._replace(
        assets=install.assets._replace(retired="No more of those."))


def test_a_retired_platform_installs_from_a_release_which_carries_its_binary(
        retired_install, fetched, monkeypatch):
    """Retirement must not break the releases cut while the platform was still
    built: they publish its binary, and have to stay installable."""
    monkeypatch.setattr(download, "load_checksums",
                        lambda tag, log=print: _checksums("Galacticus.exe"))
    monkeypatch.setattr(download, "_verify", lambda *a, **k: None)
    assert "exec" in download.provision(retired_install, log=_quiet)


def test_a_retired_platform_is_refused_by_a_release_which_carries_no_binary(
        retired_install, fetched, monkeypatch):
    monkeypatch.setattr(download, "load_checksums",
                        lambda tag, log=print: _checksums("tools.tar.zst"))
    with pytest.raises(platforms.UnsupportedPlatform) as raised:
        download.provision(retired_install, log=_quiet)
    # The reason, not a bare download failure.
    assert "No more of those." in str(raised.value)
    assert not fetched


def test_a_release_publishing_no_checksums_predates_every_retirement(
        retired_install, fetched, monkeypatch):
    """Such a release has no manifest to consult, and is old enough that the
    binary is there; the default `load_checksums` stub already returns None."""
    monkeypatch.setattr(download, "_verify", lambda *a, **k: None)
    assert "exec" in download.provision(retired_install, log=_quiet)


def test_a_retired_platform_already_installed_keeps_working(
        retired_install, fetched, monkeypatch):
    """The check gates the download of the executable, so an install provisioned
    while the binary was still published is never refused afterwards."""
    monkeypatch.setattr(download, "_verify", lambda *a, **k: None)
    download.provision(retired_install, log=_quiet)
    monkeypatch.setattr(download, "load_checksums",
                        lambda tag, log=print: _checksums("tools.tar.zst"))
    assert download.provision(retired_install, log=_quiet) == []


# --- the profile-guided optimization (PGO) build of the executable --------

_COMMIT = "0123456789abcdef0123456789abcdef01234567"
_STANDARD = "#!/bin/true\n"
_PGO = "#!/bin/true\n# the PGO build\n"


@pytest.fixture
def pgo_install(install):
    """A managed install on a platform which has a PGO build."""
    return install._replace(assets=install.assets._replace(binary_pgo="Galacticus.exe_PGO"))


@pytest.fixture
def pgo_release(release, fetched, monkeypatch):
    """Add to the fixture release a PGO executable built from the release's own commit, as `pgo.yml` publishes it."""
    binary = release / "Galacticus.exe_PGO"
    binary.write_text(_PGO)
    digest = hashlib.sha256(binary.read_bytes()).hexdigest()
    (release / "Galacticus.exe_PGO.sha256").write_text(f"{digest}  Galacticus.exe_PGO\n")
    (release / "Galacticus.exe_PGO.commit").write_text(_COMMIT + "\n")
    read = []

    def _text(url):
        read.append(url)
        source = release / Path(url).name
        return source.read_text() if source.is_file() else None

    monkeypatch.setattr(download, "_read_remote_text", _text)
    monkeypatch.setattr(download, "_release_commit", lambda tag: _COMMIT)
    return read


def _executable(install):
    return (install.exec_path / "Galacticus.exe").read_text()


def test_the_pgo_build_is_installed_where_the_release_carries_one(pgo_install, fetched, pgo_release):
    done = download.provision(pgo_install, log=_quiet)
    assert "exec" in done
    assert _executable(pgo_install) == _PGO
    assert download.pgo_installed(pgo_install)
    assert any(url.endswith("/Galacticus.exe_PGO") for url in fetched)
    assert not any(url.endswith("/Galacticus.exe") for url in fetched)


def test_a_pgo_build_failing_its_checksum_is_fatal(pgo_install, fetched, pgo_release, release):
    (release / "Galacticus.exe_PGO.sha256").write_text("0" * 64 + "  Galacticus.exe_PGO\n")
    with pytest.raises(RuntimeError, match="checksum mismatch"):
        download.provision(pgo_install, log=_quiet)
    assert not (pgo_install.exec_path / "Galacticus.exe").exists()


def test_a_pgo_build_from_another_commit_is_not_installed(pgo_install, fetched, pgo_release, monkeypatch):
    """`bleeding-edge` keeps a PGO executable from an older commit once it moves on: it must not be paired with the newer
    release's other assets."""
    monkeypatch.setattr(download, "_release_commit", lambda tag: "f" * 40)
    messages = []
    download.provision(pgo_install, log=messages.append)
    assert _executable(pgo_install) == _STANDARD
    assert not download.pgo_installed(pgo_install)
    assert any("built from commit 0123456789" in message for message in messages)


def test_the_standard_build_is_used_if_the_release_commit_is_unknown(pgo_install, fetched, pgo_release, monkeypatch):
    monkeypatch.setattr(download, "_release_commit", lambda tag: None)
    download.provision(pgo_install, log=_quiet)
    assert _executable(pgo_install) == _STANDARD


def test_the_standard_build_is_used_without_a_valid_pgo_checksum(pgo_install, fetched, pgo_release, release):
    (release / "Galacticus.exe_PGO.sha256").unlink()
    download.provision(pgo_install, log=_quiet)
    assert _executable(pgo_install) == _STANDARD


def test_a_release_without_a_pgo_build_installs_the_standard_one(pgo_install, fetched):
    download.provision(pgo_install, log=_quiet)
    assert _executable(pgo_install) == _STANDARD
    assert not download.pgo_installed(pgo_install)


def test_a_platform_without_a_pgo_build_ignores_a_published_one(install, fetched, pgo_release):
    download.provision(install, log=_quiet)
    assert _executable(install) == _STANDARD


def test_no_pgo_installs_the_standard_build_and_records_the_choice(pgo_install, fetched, pgo_release):
    download.provision(pgo_install, pgo=False, log=_quiet)
    assert _executable(pgo_install) == _STANDARD
    assert download.pgo_skipped(pgo_install)
    # An update (re-provisioning everything) keeps to the recorded choice.
    download.provision(pgo_install, force=True, log=_quiet)
    assert _executable(pgo_install) == _STANDARD


def test_declining_pgo_replaces_an_installed_pgo_build(pgo_install, fetched, pgo_release):
    download.provision(pgo_install, log=_quiet)
    assert download.pgo_installed(pgo_install)
    assert "exec" in download.provision(pgo_install, pgo=False, log=_quiet)
    assert _executable(pgo_install) == _STANDARD
    assert not download.pgo_installed(pgo_install)


def test_asking_for_pgo_replaces_an_installed_standard_build(pgo_install, fetched, pgo_release):
    download.provision(pgo_install, pgo=False, log=_quiet)
    assert "exec" in download.provision(pgo_install, pgo=True, log=_quiet)
    assert _executable(pgo_install) == _PGO
    assert not download.pgo_skipped(pgo_install)


def test_asking_for_pgo_keeps_the_standard_build_if_none_is_available(pgo_install, fetched, pgo_release,
                                                                     monkeypatch):
    download.provision(pgo_install, pgo=False, log=_quiet)
    monkeypatch.setattr(download, "_release_commit", lambda tag: "f" * 40)
    assert "exec" not in download.provision(pgo_install, pgo=True, log=_quiet)
    assert _executable(pgo_install) == _STANDARD


def test_provisioning_on_demand_does_not_look_for_a_pgo_build(pgo_install, fetched, pgo_release, release):
    """`galacticus run` provisions on demand (pgo=None): once installed, it must make no further requests."""
    (release / "Galacticus.exe_PGO.commit").unlink()
    download.provision(pgo_install, log=_quiet)       # No PGO build available: the standard one is installed.
    (release / "Galacticus.exe_PGO.commit").write_text(_COMMIT + "\n")
    pgo_release.clear()
    assert download.provision(pgo_install, log=_quiet) == []
    assert not pgo_release
    assert _executable(pgo_install) == _STANDARD


@pytest.mark.parametrize("responses, expected", [
    # A lightweight tag points straight at the commit.
    ({"ref": {"object": {"type": "commit", "sha": "a" * 40}}}, "a" * 40),
    # An annotated tag points at a tag object, which points at the commit.
    ({"ref": {"object": {"type": "tag", "sha": "b" * 40, "url": "tagObject"}},
      "tagObject": {"object": {"type": "commit", "sha": "c" * 40}}}, "c" * 40),
    # No such tag, or the API unreachable.
    ({}, None),
])
def test_release_commit_resolves_lightweight_and_annotated_tags(monkeypatch, responses, expected):
    def _api(url):
        return responses.get("ref") if "/git/ref/tags/" in url else responses.get(url)

    monkeypatch.setattr(download, "_github_api_json", _api)
    assert download._release_commit("bleeding-edge") == expected


def test_only_linux_x86_64_has_a_pgo_build():
    assert platforms.detect("Linux", "x86_64").binary_pgo == "Galacticus.exe_PGO"
    assert platforms.detect("Darwin", "arm64").binary_pgo is None

