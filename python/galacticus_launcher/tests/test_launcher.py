"""Unit tests for the Galacticus launcher.

Network-free: every test exercises the pure logic (platform mapping, install
resolution, version->tag resolution, and the cache-clean selection logic).
"""

import os

import pytest

from galacticus_launcher import platforms, paths, cli, download


# --- platform mapping ------------------------------------------------------

@pytest.mark.parametrize("system,machine,binary,fmt", [
    ("Linux", "x86_64", "Galacticus.exe", "tar.bz2"),
    ("Linux", "AMD64", "Galacticus.exe", "tar.bz2"),
    ("Darwin", "x86_64", "Galacticus_MacOS.exe", "zip"),
    ("Darwin", "arm64", "Galacticus_MacOS-M1.exe", "zip"),
    ("Darwin", "aarch64", "Galacticus_MacOS-M1.exe", "zip"),
])
def test_detect_supported(system, machine, binary, fmt):
    assets = platforms.detect(system, machine)
    assert assets.binary == binary
    assert assets.tools_format == fmt


@pytest.mark.parametrize("system,machine", [
    ("Linux", "ppc64le"),
    ("Darwin", "i386"),
    ("Windows", "AMD64"),
    ("SunOS", "sparc"),
])
def test_detect_unsupported(system, machine):
    with pytest.raises(platforms.UnsupportedPlatform):
        platforms.detect(system, machine)


# --- version -> release tag ------------------------------------------------

@pytest.mark.parametrize("version,expected", [
    ("1.2.3", "v1.2.3"),
    ("0.0.0", "bleeding-edge"),          # development placeholder
    ("1.2.3.dev4+g999", "bleeding-edge"),
    ("0.1.dev575+g58c2.d20260626", "bleeding-edge"),
    ("1.2", "bleeding-edge"),            # not X.Y.Z
])
def test_release_tag(monkeypatch, version, expected):
    monkeypatch.delenv("GALACTICUS_RELEASE_TAG", raising=False)
    assert paths.release_tag(version) == expected


def test_release_tag_override(monkeypatch):
    monkeypatch.setenv("GALACTICUS_RELEASE_TAG", "v9.9.9")
    assert paths.release_tag("1.2.3") == "v9.9.9"


# --- install resolution ----------------------------------------------------

def _clear_galacticus_env(monkeypatch):
    for name in ("GALACTICUS_HOME", "GALACTICUS_EXEC_PATH", "GALACTICUS_DATA_PATH",
                 "GALACTICUS_TOOLS_PATH", "GALACTICUS_DYNAMIC_DATA_PATH",
                 "GALACTICUS_RELEASE_TAG"):
        monkeypatch.delenv(name, raising=False)


def test_resolve_managed(monkeypatch):
    _clear_galacticus_env(monkeypatch)
    monkeypatch.setattr(platforms, "detect",
                        lambda: platforms.PlatformAssets("Galacticus.exe", "tools.tar.bz2", "tar.bz2", "test"))
    install = paths.resolve("1.2.3")
    assert install.source == paths.SOURCE_MANAGED
    assert install.managed
    assert install.tag == "v1.2.3"
    # Tools live under the durable data root; dynamic under the cache root.
    assert install.tools_path.parent == install.exec_path.parent
    assert install.dynamic_path != install.tools_path
    assert install.tools_separated  # managed always separates tools from cache
    env = install.environ()
    assert env["GALACTICUS_TOOLS_PATH"] == str(install.tools_path)
    assert env["GALACTICUS_DYNAMIC_DATA_PATH"] == str(install.dynamic_path)


def test_resolve_environment(monkeypatch, tmp_path):
    _clear_galacticus_env(monkeypatch)
    exec_dir = tmp_path / "build"
    exec_dir.mkdir()
    (exec_dir / "Galacticus.exe").write_text("#!/bin/true\n")
    data_dir = tmp_path / "datasets"
    data_dir.mkdir()
    monkeypatch.setenv("GALACTICUS_EXEC_PATH", str(exec_dir))
    monkeypatch.setenv("GALACTICUS_DATA_PATH", str(data_dir))
    install = paths.resolve()
    assert install.source == paths.SOURCE_ENVIRONMENT
    assert not install.managed
    assert install.binary == exec_dir / "Galacticus.exe"
    # Dynamic data path defaults to <data>/dynamic when unset, so info/clean
    # reflect where the binary actually writes regenerable data.
    assert install.dynamic_path == data_dir / "dynamic"
    # Tools are not separated (no GALACTICUS_TOOLS_PATH), so they share dynamic.
    assert install.tools_path is None
    assert not install.tools_separated
    # environ() must NOT impose tools/dynamic in environment mode.
    env = install.environ({"PATH": "/usr/bin"})
    assert "GALACTICUS_TOOLS_PATH" not in env
    assert "GALACTICUS_DYNAMIC_DATA_PATH" not in env


def test_resolve_environment_separated_tools(monkeypatch, tmp_path):
    _clear_galacticus_env(monkeypatch)
    exec_dir = tmp_path / "build"
    exec_dir.mkdir()
    (exec_dir / "Galacticus.exe").write_text("#!/bin/true\n")
    (tmp_path / "datasets").mkdir()
    tools_dir = tmp_path / "tools"
    monkeypatch.setenv("GALACTICUS_EXEC_PATH", str(exec_dir))
    monkeypatch.setenv("GALACTICUS_DATA_PATH", str(tmp_path / "datasets"))
    monkeypatch.setenv("GALACTICUS_TOOLS_PATH", str(tools_dir))
    install = paths.resolve()
    assert install.tools_separated  # explicit tools path differs from dynamic


def test_resolve_environment_without_binary_falls_through(monkeypatch, tmp_path):
    _clear_galacticus_env(monkeypatch)
    monkeypatch.setattr(platforms, "detect",
                        lambda: platforms.PlatformAssets("Galacticus.exe", "tools.tar.bz2", "tar.bz2", "test"))
    monkeypatch.setenv("GALACTICUS_EXEC_PATH", str(tmp_path / "nope"))
    monkeypatch.setenv("GALACTICUS_DATA_PATH", str(tmp_path / "datasets"))
    install = paths.resolve("1.2.3")
    assert install.source == paths.SOURCE_MANAGED  # no binary -> managed


def test_resolve_home(monkeypatch, tmp_path):
    _clear_galacticus_env(monkeypatch)
    home = tmp_path / "galacticus"
    home.mkdir()
    (home / "Galacticus.exe").write_text("#!/bin/true\n")
    (tmp_path / "datasets" / "static").mkdir(parents=True)
    monkeypatch.setenv("GALACTICUS_HOME", str(home))
    install = paths.resolve()
    assert install.source == paths.SOURCE_HOME
    assert install.exec_path == home
    # Datasets discovered from the sibling directory.
    assert install.data_path == tmp_path / "datasets"


# --- url builders ----------------------------------------------------------

def test_source_url():
    assert download.source_url("bleeding-edge").endswith("/heads/master.zip")
    assert download.source_url("v1.2.3").endswith("/tags/v1.2.3.zip")


def test_asset_url():
    url = download.asset_url("v1.2.3", "tools.tar.bz2")
    assert url == ("https://github.com/galacticusorg/galacticus/releases/"
                   "download/v1.2.3/tools.tar.bz2")


def test_datasets_url():
    assert download.datasets_url().endswith("/datasets/archive/master.zip")
    assert download.datasets_url("deadbeef").endswith("/datasets/archive/deadbeef.zip")


def test_resolve_datasets_ref_env_override(monkeypatch):
    monkeypatch.setenv("GALACTICUS_DATASETS_REF", "my-branch")
    assert download.resolve_datasets_ref("v1.2.3") == "my-branch"


def test_resolve_datasets_ref_bleeding(monkeypatch):
    monkeypatch.delenv("GALACTICUS_DATASETS_REF", raising=False)
    # bleeding-edge tracks master and never hits the network for a pin.
    assert download.resolve_datasets_ref("bleeding-edge") == "master"


def test_resolve_datasets_ref_pinned(monkeypatch):
    monkeypatch.delenv("GALACTICUS_DATASETS_REF", raising=False)
    monkeypatch.setattr(download, "_read_remote_text", lambda url: "abc123def\n")
    assert download.resolve_datasets_ref("v1.2.3") == "abc123def"


def test_resolve_datasets_ref_unpinned_falls_back(monkeypatch):
    monkeypatch.delenv("GALACTICUS_DATASETS_REF", raising=False)
    monkeypatch.setattr(download, "_read_remote_text", lambda url: None)
    assert download.resolve_datasets_ref("v1.2.3") == "master"


# --- cache cleaning --------------------------------------------------------

def test_purge_tree_older_than(tmp_path):
    old = tmp_path / "old.dat"
    new = tmp_path / "new.dat"
    old.write_bytes(b"x" * 100)
    new.write_bytes(b"y" * 50)
    # Age `old` to 10 days; keep `new` fresh.
    ten_days = 10 * 86400
    now = os.path.getmtime(new)
    os.utime(old, (now - ten_days, now - ten_days))
    # Cutoff = 5 days ago: only `old` qualifies.
    cutoff = now - 5 * 86400
    freed, removed = cli._purge_tree(tmp_path, cutoff=cutoff, dry_run=True)
    assert removed == 1 and freed == 100
    assert old.exists()  # dry-run deletes nothing
    freed, removed = cli._purge_tree(tmp_path, cutoff=cutoff, dry_run=False)
    assert removed == 1 and freed == 100
    assert not old.exists() and new.exists()


def test_purge_tree_all(tmp_path):
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "a.dat").write_bytes(b"a" * 10)
    (tmp_path / "b.dat").write_bytes(b"b" * 20)
    freed, removed = cli._purge_tree(tmp_path, cutoff=None, dry_run=False)
    assert freed == 30
    assert removed == 2  # one dir + one file at top level
    assert list(tmp_path.iterdir()) == []


def test_human_readable():
    assert cli._human(0) == "0 B"
    assert cli._human(1023) == "1023 B"
    assert cli._human(1024).endswith("KiB")
    assert cli._human(1024 ** 3).endswith("GiB")
