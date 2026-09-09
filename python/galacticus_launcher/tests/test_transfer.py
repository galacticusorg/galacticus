"""Unit tests for how the launcher moves bytes: byte-range splitting, the
choice not to split, and running several downloads at once.

Network-free: ``requests.get`` is replaced by a fake origin server which serves
one byte string and honours ``Range`` the way GitHub's release assets do (or
refuses to, the way its on-the-fly repository archives do).
"""

import threading

import pytest

pytest.importorskip("requests")

import requests

from galacticus_launcher import download


class _FakeResponse:
    """Just enough of `requests.Response` for the download path."""

    def __init__(self, status_code, body, headers):
        self.status_code = status_code
        self.headers = headers
        self._body = body

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"status {self.status_code}")

    def iter_content(self, chunk_size=1):
        for start in range(0, len(self._body), chunk_size):
            yield self._body[start:start + chunk_size]

    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def _origin(body, *, ranges=True, missing=False):
    """A fake ``requests.get`` serving `body`, recording the ranges asked for."""
    asked = []

    def get(url, stream=False, timeout=None, headers=None):
        requested = (headers or {}).get("Range")
        asked.append(requested)
        if missing:
            return _FakeResponse(404, b"", {})
        if ranges and requested:
            span = requested.split("=", 1)[1]
            start_text, _, end_text = span.partition("-")
            start = int(start_text)
            end = int(end_text) if end_text else len(body) - 1
            end = min(end, len(body) - 1)
            chunk = body[start:end + 1]
            return _FakeResponse(206, chunk, {
                "Content-Range": f"bytes {start}-{end}/{len(body)}",
                "Content-Length": str(len(chunk)),
            })
        # A server which does not do ranges answers the whole body, whatever was
        # asked for -- which is exactly what makes the "ask for bytes=0-" probe
        # free: this response *is* the transfer.
        return _FakeResponse(200, body, {"Content-Length": str(len(body))})

    return get, asked


@pytest.fixture
def quiet_progress():
    return download._Display(log=lambda message: None).add("asset")


# --- reading the size off the response --------------------------------------

def test_total_length_prefers_the_content_range_total():
    """On a 206, Content-Length is the length of the *range*; only Content-Range
    states the size of the whole resource."""
    response = _FakeResponse(206, b"", {"Content-Range": "bytes 0-9/1000",
                                        "Content-Length": "10"})
    assert download._total_length(response) == 1000


def test_total_length_falls_back_to_content_length():
    assert download._total_length(_FakeResponse(200, b"", {"Content-Length": "42"})) == 42


def test_total_length_is_none_when_unadvertised():
    assert download._total_length(_FakeResponse(200, b"", {})) is None


# --- splitting a large asset across connections -----------------------------

def test_a_large_asset_is_split_and_reassembled(tmp_path, monkeypatch,
                                                quiet_progress):
    body = bytes(range(256)) * 400                       # 102 400 bytes
    get, asked = _origin(body)
    monkeypatch.setattr(download.requests, "get", get)
    monkeypatch.setattr(download, "_SPLIT_THRESHOLD", 1024)
    monkeypatch.setenv("GALACTICUS_DOWNLOAD_CONNECTIONS", "4")

    dest = tmp_path / "asset.bin"
    assert download._download("https://example/asset", dest,
                              progress=quiet_progress) is True
    assert dest.read_bytes() == body
    # One probe, then one request per span; together the spans cover the file
    # exactly once, with no overlap and nothing left out.
    spans = [request for request in asked if request != "bytes=0-"]
    assert len(spans) == 4
    bounds = sorted(tuple(int(value) for value in span.split("=")[1].split("-"))
                    for span in spans)
    assert bounds[0][0] == 0 and bounds[-1][1] == len(body) - 1
    for (_, end), (start, _) in zip(bounds, bounds[1:]):
        assert start == end + 1


def test_a_small_asset_is_not_split(tmp_path, monkeypatch, quiet_progress):
    body = b"small payload"
    get, asked = _origin(body)
    monkeypatch.setattr(download.requests, "get", get)
    monkeypatch.setattr(download, "_SPLIT_THRESHOLD", 1024)

    dest = tmp_path / "asset.bin"
    download._download("https://example/asset", dest, progress=quiet_progress)
    assert dest.read_bytes() == body
    assert asked == ["bytes=0-"]      # the probe response was the transfer


def test_a_server_without_ranges_is_streamed(tmp_path, monkeypatch,
                                             quiet_progress):
    """The repository archives are generated on the fly and answer 200 to a
    range request; the body of that answer has to be used, not discarded."""
    body = bytes(range(256)) * 400
    get, asked = _origin(body, ranges=False)
    monkeypatch.setattr(download.requests, "get", get)
    monkeypatch.setattr(download, "_SPLIT_THRESHOLD", 1024)

    dest = tmp_path / "asset.bin"
    download._download("https://example/asset", dest, progress=quiet_progress)
    assert dest.read_bytes() == body
    assert asked == ["bytes=0-"]


def test_splitting_can_be_turned_off(tmp_path, monkeypatch, quiet_progress):
    body = bytes(range(256)) * 400
    get, asked = _origin(body)
    monkeypatch.setattr(download.requests, "get", get)
    monkeypatch.setattr(download, "_SPLIT_THRESHOLD", 1024)
    monkeypatch.setenv("GALACTICUS_DOWNLOAD_CONNECTIONS", "1")

    dest = tmp_path / "asset.bin"
    download._download("https://example/asset", dest, progress=quiet_progress)
    assert dest.read_bytes() == body
    assert asked == ["bytes=0-"]


def test_connections_are_clamped(monkeypatch):
    monkeypatch.setenv("GALACTICUS_DOWNLOAD_CONNECTIONS", "0")
    assert download._connections() == 1
    monkeypatch.setenv("GALACTICUS_DOWNLOAD_CONNECTIONS", "1000")
    assert download._connections() == download._MAX_CONNECTIONS
    monkeypatch.setenv("GALACTICUS_DOWNLOAD_CONNECTIONS", "not a number")
    assert download._connections() == download._CONNECTIONS
    monkeypatch.delenv("GALACTICUS_DOWNLOAD_CONNECTIONS")
    assert download._connections() == download._CONNECTIONS


def test_a_missing_asset_is_reported_rather_than_retried(tmp_path, monkeypatch,
                                                         quiet_progress):
    get, asked = _origin(b"", missing=True)
    monkeypatch.setattr(download.requests, "get", get)
    assert download._download("https://example/absent", tmp_path / "x",
                              missing_ok=True, progress=quiet_progress) is False
    assert len(asked) == 1            # a 404 is an answer, not a failure to retry


# --- running the components' downloads together -----------------------------

def test_fetcher_runs_its_jobs_concurrently(tmp_path, monkeypatch):
    """The components are downloaded at the same time; the endpoint throttles
    each connection, so overlapping them is what makes an install quick."""
    started = threading.Barrier(3, timeout=10)

    def fake_download(url, dest, *, log=print, missing_ok=False, progress=None):
        started.wait()                # only returns once all three are in flight
        dest.write_text(url)
        return True

    monkeypatch.setattr(download, "_download", fake_download)
    fetcher = download._Fetcher(log=lambda message: None)
    jobs = [fetcher.enqueue(f"https://example/{index}", tmp_path / str(index),
                            label=str(index)) for index in range(3)]
    fetcher.run()                     # would time out if the jobs were serial
    assert all(job.present for job in jobs)
    assert (tmp_path / "1").read_text() == "https://example/1"


def test_fetcher_reports_a_failed_job(tmp_path, monkeypatch):
    def fake_download(url, dest, *, log=print, missing_ok=False, progress=None):
        if url.endswith("2"):
            raise RuntimeError("no such asset")
        return True

    monkeypatch.setattr(download, "_download", fake_download)
    fetcher = download._Fetcher(log=lambda message: None)
    for index in range(3):
        fetcher.enqueue(f"https://example/{index}", tmp_path / str(index),
                        label=str(index))
    with pytest.raises(RuntimeError, match="no such asset"):
        fetcher.run()


def test_an_optional_job_that_fails_does_not_sink_the_others(tmp_path,
                                                             monkeypatch):
    """The components are fetched together into staging directories which are
    torn down if the fetch phase raises. An artefact the install can rebuild for
    itself must therefore not raise -- otherwise one small failed download
    discards the gigabytes which did arrive."""
    def fake_download(url, dest, *, log=print, missing_ok=False, progress=None):
        if url.endswith("optional"):
            raise RuntimeError("connection reset")
        dest.write_text(url)
        return True

    monkeypatch.setattr(download, "_download", fake_download)
    messages = []
    fetcher = download._Fetcher(log=messages.append)
    needed = fetcher.enqueue("https://example/needed", tmp_path / "needed",
                             label="needed")
    spare = fetcher.enqueue("https://example/optional", tmp_path / "spare",
                            label="optional", optional=True)
    fetcher.run()
    assert needed.present and (tmp_path / "needed").is_file()
    assert not spare.present and isinstance(spare.error, RuntimeError)
    assert any("could not fetch optional" in message for message in messages)
