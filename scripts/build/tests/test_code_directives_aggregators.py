"""Regression test for the Makefile dependencies of aggregating directives.

Some directives are expanded by the preprocessor from the *content* of every
file that carries another directive: `functionsGlobal` generates a pointer for
each `functionGlobal`, `eventHookManager` reads every `eventHook`,
`metaPropertyDatabase` every `addMetaProperty`, and each `eventHookStatic`
calls every function attached by a directive named after the hook. The
generic `%.p.F90.up` rule only depends on `directiveLocations.xml`, which
changes when a file gains or loses a directive, so adding a second
`functionGlobal` to a file which already had one left `Functions_Global`
stale, and the build failed with "Explicit interface required".

This drives codeDirectivesParse.main() over a fixture tree and requires
`Makefile_Directives` to make each consumer depend on the files it reads.
"""

import os
import re
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir))
import codeDirectivesParse as cdp                       # noqa: E402


_FILES = {
    'pointers.F90': '''module Pointers
  !![
  <functionsGlobal type="pointers"/>
  !!]
end module Pointers
''',
    'producer.F90': '''module Producer
contains
  !![
  <functionGlobal>
   <unitName>Producer_Function</unitName>
   <type>void</type>
  </functionGlobal>
  !!]
  subroutine Producer_Function()
  end subroutine Producer_Function
end module Producer
''',
    'hookCaller.F90': '''module Hook_Caller
contains
  subroutine Hook_Caller_Call()
    !![
    <eventHookStatic name="probeTask"/>
    !!]
  end subroutine Hook_Caller_Call
end module Hook_Caller
''',
    'hooked.F90': '''module Hooked
  !![
  <probeTask function="Hooked_Task"/>
  !!]
end module Hooked
''',
    'unrelated.F90': '''module Unrelated
  !![
  <outputFileClose function="Unrelated_Report"/>
  !!]
end module Unrelated
''',
}


def _prerequisites(makefile, target):
    """Return the set of prerequisites given for `target` in `makefile`."""
    prerequisites = set()
    for line in makefile.splitlines():
        match = re.match(r'^(\S+):\s*(.*)$', line)
        if match and os.path.normpath(match.group(1)) == os.path.normpath(target):
            prerequisites.update(
                os.path.basename(p) for p in match.group(2).split())
    return prerequisites


def test_aggregator_consumers_depend_on_producers(tmp_path, monkeypatch):
    source = tmp_path / 'source'
    build  = tmp_path / 'build'
    source.mkdir()
    build.mkdir()
    monkeypatch.setenv('BUILDPATH', str(build))
    monkeypatch.setenv('GALACTICUS_BUILD_JOBS', '1')   # deterministic, no fork
    for name, content in _FILES.items():
        (source / name).write_text(content)
    cdp.main(['codeDirectivesParse.py', str(tmp_path)])
    makefile = (build / 'Makefile_Directives').read_text()

    # The functionsGlobal consumer depends on the functionGlobal producer only.
    assert _prerequisites(
        makefile, str(build / 'pointers.p.F90.up')) == {'producer.F90'}
    # A static event hook depends on the files attached to that hook only.
    assert _prerequisites(
        makefile, str(build / 'hookCaller.p.F90.up')) == {'hooked.F90'}
    # Files which consume nothing gain no such rule.
    assert _prerequisites(makefile, str(build / 'unrelated.p.F90.up')) == set()
