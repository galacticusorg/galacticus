"""Pin codeDirectivesParse's `_canonicalize` pickle-identity invariant.

The per-file blob cache is rewritten only-if-changed, and "changed" is a
byte comparison of the pickle. Entries built in-process share the interned
literal strings (`'files'`, `'includeDirectives'`, ...) across entries, so
pickle writes them once with back-references; entries returned from forked
pool workers are unpickled into fresh, de-interned strings, which would
make the merged blob pickle differently from a serial run — bumping the
blob mtime on every scan and defeating the freshness window.
`_canonicalize` re-interns exactly the literal set `_ENTRY_LITERALS` to
restore byte-identity. That set must be kept in sync with the literals
`_scan_one` produces — an obligation this test enforces: it drives the real
`_scan_one` over fixtures exercising every top-level entry shape (a
functionClass directive with its implicit task dependencies, plus an
include directive), simulates the worker round-trip, and requires the
canonicalized blob to pickle byte-identically to the serial one.
"""

import os
import pickle
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir))
import codeDirectivesParse as cdp  # noqa: E402


_PROBE = '''!![
<functionClass docformat="rst">
 <name>{name}Widget</name>
 <descriptiveName>Widgets</descriptiveName>
</functionClass>
!!]
module {name}_Widgets
  !![
  <include directive="{name}Inc" type="function">
  !!]
  include '{name}.widgets.inc'
  !![
  </include>
  !!]
end module {name}_Widgets
'''


def _entries(tmp_path, roundtrip):
    source = tmp_path / 'source'
    source.mkdir(exist_ok=True)
    (tmp_path / 'build').mkdir(exist_ok=True)
    cdp._WORKER['source_directory'] = str(source)
    cdp._WORKER['build_path']       = str(tmp_path / 'build')
    blob = {}
    for name in ('test', 'other'):
        path = source / f'{name}.F90'
        path.write_text(_PROBE.format(name=name))
        identifier, entry = cdp._scan_one((f'{name}.F90', str(path),
                                           f'id_{name}'))
        if roundtrip:
            # Simulate the per-task pickle round-trip through the process
            # pool: strings lose their interned identity.
            entry = pickle.loads(pickle.dumps(entry))
        blob[identifier] = entry
    return blob


def test_fixture_exercises_full_entry_shape(tmp_path):
    entry = next(iter(_entries(tmp_path, roundtrip=False).values()))
    assert sorted(entry) == ['files', 'functionClasses',
                             'includeDirectives', 'nonIncludeDirectives']


def test_canonicalized_worker_blob_pickles_identically(tmp_path):
    serial = _entries(tmp_path, roundtrip=False)
    worker = _entries(tmp_path, roundtrip=True)
    # Guard: the round-trip really does de-intern (otherwise this test —
    # and _canonicalize itself — would be vacuous; if this ever fails,
    # pickle semantics changed and _canonicalize may be removable).
    assert pickle.dumps(worker) != pickle.dumps(serial)
    canonical = {k: cdp._canonicalize(v) for k, v in worker.items()}
    assert pickle.dumps(canonical) == pickle.dumps(serial)
