"""Unit tests for the integrity checks applied to downloaded artefacts.

Network-free: the published checksums and the archives are constructed locally.

These cover the two properties which stand between a bad download and code being executed on the user's machine:

* an artefact whose content does not match the checksum published by its release is rejected, and removed rather than left on
  disk where a later run might use it; and
* unpacking an archive can not write outside the directory it is being unpacked into.
"""

import hashlib
import io
import tarfile

import pytest

pytest.importorskip("platformdirs")
pytest.importorskip("requests")

from galacticus_launcher import download


# --- checksum verification -------------------------------------------------

def _write(path, content=b"pretend artefact"):
    path.write_bytes(content)
    return hashlib.sha256(content).hexdigest()


def test_matching_checksum_accepted(tmp_path):
    target = tmp_path / "Galacticus.exe"
    digest = _write(target)
    download._verify(target, "Galacticus.exe", {"Galacticus.exe": digest}, log=lambda message: None)
    assert target.exists()


def test_mismatched_checksum_is_fatal_and_removes_the_file(tmp_path):
    target = tmp_path / "Galacticus.exe"
    _write(target)
    with pytest.raises(RuntimeError, match="checksum mismatch"):
        download._verify(target, "Galacticus.exe", {"Galacticus.exe": "0" * 64}, log=lambda message: None)
    assert not target.exists(), "an artefact which failed verification must not be left on disk"


def test_release_without_checksums_still_provisions(tmp_path):
    """Releases published before `SHA256SUMS` existed must remain installable."""
    target = tmp_path / "Galacticus.exe"
    _write(target)
    download._verify(target, "Galacticus.exe", None, log=lambda message: None)
    assert target.exists()


def test_missing_entry_for_this_asset_is_reported(tmp_path):
    target = tmp_path / "Galacticus.exe"
    digest = _write(target)
    messages = []
    download._verify(target, "Galacticus.exe", {"tools.tar.bz2": digest}, log=messages.append)
    assert any("no checksum" in message for message in messages)


@pytest.mark.parametrize("line,expected", [
    ("{d}  Galacticus.exe", "Galacticus.exe"),   # text mode, as written by `sha256sum`
    ("{d} *tools.tar.bz2", "tools.tar.bz2"),     # binary mode marker
])
def test_checksum_file_parsing(monkeypatch, line, expected):
    digest = "a" * 64
    monkeypatch.setattr(download, "_read_remote_text", lambda url: line.format(d=digest) + "\n")
    assert download.load_checksums("v1.0.0") == {expected: digest}


def test_absent_checksum_asset_yields_none(monkeypatch):
    monkeypatch.setattr(download, "_read_remote_text", lambda url: None)
    assert download.load_checksums("v1.0.0") is None


# --- archive extraction ----------------------------------------------------

def _tar_with_member(path, name, content=b"payload"):
    with tarfile.open(path, "w:gz") as archive:
        info = tarfile.TarInfo(name)
        info.size = len(content)
        archive.addfile(info, io.BytesIO(content))


def test_extract_refuses_to_escape_the_destination(tmp_path):
    archive = tmp_path / "evil.tar.gz"
    _tar_with_member(archive, "../ESCAPED")
    with pytest.raises(Exception):
        download._extract(archive, tmp_path / "dest", "tar.gz")
    assert not (tmp_path / "ESCAPED").exists()


def test_extract_unpacks_a_benign_archive(tmp_path):
    archive = tmp_path / "tools.tar.gz"
    _tar_with_member(archive, "dynamic/camb", b"fine")
    destination = tmp_path / "dest"
    download._extract(archive, destination, "tar.gz")
    assert (destination / "dynamic" / "camb").read_bytes() == b"fine"
