"""Tests for Galacticus.Build.FileChanges.update — the only-if-changed file
replacement that (with its ``.up`` sentinels) is the backbone of the build's
incremental minimality: a generated file's mtime must advance only when its
content actually changes.
"""

import os

from Galacticus.Build.FileChanges import update


def _write(path, text):
    path.write_text(text)
    return str(path)


def test_missing_old_moves_new_into_place(tmp_path):
    new = _write(tmp_path / 'x.tmp', 'content')
    old = tmp_path / 'x'
    update(str(old), new)
    assert old.read_text() == 'content'
    assert not os.path.exists(new)


def test_identical_preserves_mtime_and_removes_new(tmp_path):
    old = _write(tmp_path / 'x', 'same')
    os.utime(old, (1000000, 1000000))          # ancient mtime
    new = _write(tmp_path / 'x.tmp', 'same')
    update(old, new)
    assert os.stat(old).st_mtime == 1000000    # untouched
    assert not os.path.exists(new)


def test_different_replaces_and_advances_mtime(tmp_path):
    old = _write(tmp_path / 'x', 'before')
    os.utime(old, (1000000, 1000000))
    new = _write(tmp_path / 'x.tmp', 'after')
    update(old, new)
    assert (tmp_path / 'x').read_text() == 'after'
    assert os.stat(old).st_mtime > 1000000
    assert not os.path.exists(new)


def test_prove_update_touches_sentinel_even_when_unchanged(tmp_path):
    # The `.up` sentinel is make's evidence that the generator ran; it must
    # advance on EVERY run, including no-op ones where the output's own
    # mtime is deliberately preserved.
    old = _write(tmp_path / 'x', 'same')
    sentinel = tmp_path / 'x.up'
    _write(sentinel, '')
    os.utime(old, (1000000, 1000000))
    os.utime(sentinel, (1000000, 1000000))
    new = _write(tmp_path / 'x.tmp', 'same')
    update(old, new, prove_update=True)
    assert os.stat(old).st_mtime == 1000000          # output untouched
    assert os.stat(sentinel).st_mtime > 1000000      # sentinel advanced


def test_prove_update_creates_sentinel(tmp_path):
    new = _write(tmp_path / 'x.tmp', 'content')
    update(str(tmp_path / 'x'), new, prove_update=True)
    assert (tmp_path / 'x.up').exists()
