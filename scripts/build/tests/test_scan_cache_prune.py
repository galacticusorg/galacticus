"""Regression test for stale-entry pruning in the per-file blob caches.

Deleting a source file that carries a directive used to leave its entry in
`codeDirectives.blob` forever: the removal was *detected* (forcing a global
rescan) but a rescan only ever replaces entries for files that still exist,
so the reduction re-emitted the dead file's directives into
`directiveLocations.xml` and `Makefile_Directives` — a phantom make
dependency on a file that no longer exists ("No rule to make target ...").
Deleting the generated fragments did not help: they were regenerated from
the same stale cache. `ScanCache.prune` (called by every catalog script
after it computes its current file set) drops such entries.

This drives codeDirectivesParse.main() end-to-end over a fixture tree:
scan, delete one directive-bearing file, rescan, and require the deleted
file to be gone from the blob, the directive catalog, and the Makefile
fragment.
"""

import os
import pickle
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir))
import codeDirectivesParse as cdp                       # noqa: E402
from Galacticus.Build.ScanCache import (                # noqa: E402
    file_identifier,
    prune,
)


_PROBE = '''module {name}_Module
  !![
  <outputFileClose function="{name}_Report"/>
  !!]
end module {name}_Module
'''


def _run(install):
    cdp.main(['codeDirectivesParse.py', str(install)])


def _setup_tree(tmp_path, monkeypatch, names):
    source = tmp_path / 'source'
    build  = tmp_path / 'build'
    source.mkdir()
    build.mkdir()
    monkeypatch.setenv('BUILDPATH', str(build))
    monkeypatch.setenv('GALACTICUS_BUILD_JOBS', '1')   # deterministic, no fork
    for name in names:
        (source / f'{name}.F90').write_text(_PROBE.format(name=name))
    return source, build


def _blob(build):
    with open(build / 'codeDirectives.blob', 'rb') as fh:
        return pickle.load(fh)


def test_deleted_file_is_pruned_everywhere(tmp_path, monkeypatch):
    source, build = _setup_tree(tmp_path, monkeypatch, ('keep', 'doomed'))
    _run(tmp_path)

    doomed_path = source / 'doomed.F90'
    doomed_id   = file_identifier(str(doomed_path))
    assert doomed_id in _blob(build)          # cached while the file existed

    doomed_path.unlink()
    _run(tmp_path)

    # Gone from the cache, the directive catalog, and the Makefile fragment.
    assert doomed_id not in _blob(build)
    assert 'doomed' not in (build / 'directiveLocations.xml').read_text()
    assert 'doomed' not in (build / 'Makefile_Directives').read_text()

    # The surviving file is untouched.
    keep_id = file_identifier(str(source / 'keep.F90'))
    assert keep_id in _blob(build)
    assert 'keep.F90' in (build / 'directiveLocations.xml').read_text()


def test_unchanged_rescan_prunes_nothing(tmp_path, monkeypatch):
    source, build = _setup_tree(tmp_path, monkeypatch, ('one', 'two'))
    _run(tmp_path)
    before = _blob(build)
    _run(tmp_path)
    assert _blob(build) == before


def test_prune_helper_respects_reserved_keys():
    cache = {'a': 1, 'b': 2, 'blobVersion': 3}
    stale = prune(cache, {'a'}, reserved={'blobVersion'})
    assert cache == {'a': 1, 'blobVersion': 3}
    assert stale == ['b']
