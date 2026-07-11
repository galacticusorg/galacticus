"""Tests for Galacticus.Build.ScanCache — the shared per-file blob-cache
helpers: cache keys must be stable (they key every incremental rescan) and
cache loading must degrade to a full rescan on any problem, never fail.
"""

import os
import pickle

from Galacticus.Build.ScanCache import file_identifier, load_cache


def test_identifier_replaces_slashes():
    assert file_identifier('a/b/c.F90') == 'a_b_c.F90'


def test_identifier_strips_leading_dot_slash():
    # `./`-relative invocations must produce the same key as bare-relative
    # ones (the strip is a no-op for the absolute paths make passes).
    assert file_identifier('./a/b.F90') == 'a_b.F90'
    assert file_identifier('a/b.F90') == file_identifier('./a/b.F90')


def test_identifier_absolute_path():
    assert file_identifier('/x/y.F90') == '_x_y.F90'


def test_load_missing(tmp_path):
    cache, mtime = load_cache(str(tmp_path / 'nope.blob'))
    assert cache == {} and mtime is None


def test_load_corrupt(tmp_path):
    blob = tmp_path / 'x.blob'
    blob.write_bytes(b'not a pickle \x00\x01')
    cache, mtime = load_cache(str(blob))
    assert cache == {} and mtime is None


def test_load_non_dict(tmp_path):
    blob = tmp_path / 'x.blob'
    blob.write_bytes(pickle.dumps(['a', 'list']))
    cache, mtime = load_cache(str(blob))
    assert cache == {} and mtime is None


def test_load_valid_returns_dict_and_mtime(tmp_path):
    blob = tmp_path / 'x.blob'
    blob.write_bytes(pickle.dumps({'k': 1}))
    cache, mtime = load_cache(str(blob))
    assert cache == {'k': 1}
    assert mtime == os.stat(blob).st_mtime
