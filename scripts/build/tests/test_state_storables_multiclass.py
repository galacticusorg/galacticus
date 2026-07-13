"""Regression test for stateStorables.py's Pass 2 merge.

Pass 2 builds one scan task per (instance-file, functionClass) pair. A
single source file that declares concrete instances of *two* different
functionClasses yields two tasks for that file. The merge previously
popped the file's cache entry before every task, so the second task's pop
discarded the first task's instance names — silently dropping that file's
instances of every class but the last from the emitted
`functionClassInstances` list (which drives state store/restore code
generation). The fix pops each instance file's entry only once per run, so
sibling tasks accumulate.

Latent in the current tree (no source file implements instances of two
functionClasses), so this drives stateStorables.main() end-to-end over a
fixture tree that does.
"""

import os
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir))
import stateStorables as ss  # noqa: E402


def _emitted_instances(tmp_path, monkeypatch, *, share_file):
    """Run stateStorables over a fixture tree declaring two functionClasses
    (classA, classB) whose instances live either in one shared file or in
    two separate files, and return the emitted functionClassInstances list.
    """
    install = tmp_path
    source  = install / 'source'
    build   = install / 'build'
    source.mkdir()
    build.mkdir()
    monkeypatch.setenv('BUILDPATH', str(build))
    monkeypatch.setenv('GALACTICUS_BUILD_JOBS', '1')   # deterministic, no fork

    base = source / 'base.F90'
    base.write_text(
        'module Bases\n'
        '  !![\n  <functionClass>\n   <name>classA</name>\n  </functionClass>\n  !!]\n'
        '  !![\n  <functionClass>\n   <name>classB</name>\n  </functionClass>\n  !!]\n'
        'end module Bases\n'
    )

    a_block = '  !![\n  <classA name="classAImpl"/>\n  !!]\n'
    b_block = '  !![\n  <classB name="classBImpl"/>\n  !!]\n'
    if share_file:
        (source / 'impls.F90').write_text(
            'module Impls\n' + a_block + b_block + 'end module Impls\n')
        a_file = b_file = source / 'impls.F90'
    else:
        (source / 'implsA.F90').write_text(
            'module ImplsA\n' + a_block + 'end module ImplsA\n')
        (source / 'implsB.F90').write_text(
            'module ImplsB\n' + b_block + 'end module ImplsB\n')
        a_file = source / 'implsA.F90'
        b_file = source / 'implsB.F90'

    (build / 'directiveLocations.xml').write_text(
        '<directiveLocations>\n'
        f'  <functionClass><file>{base}</file></functionClass>\n'
        f'  <classA><file>{a_file}</file></classA>\n'
        f'  <classB><file>{b_file}</file></classB>\n'
        '</directiveLocations>\n'
    )

    ss.main(['stateStorables.py', str(install)])
    root = ET.parse(build / 'stateStorables.xml').getroot()
    return [e.text for e in root.findall('functionClassInstances')]


def test_two_classes_in_one_file_keep_both_instances(tmp_path, monkeypatch):
    names = _emitted_instances(tmp_path, monkeypatch, share_file=True)
    assert 'classAImpl' in names
    assert 'classBImpl' in names


def test_two_classes_in_separate_files_unaffected(tmp_path, monkeypatch):
    # The one-class-per-file case (every real file today) must be unchanged.
    names = _emitted_instances(tmp_path, monkeypatch, share_file=False)
    assert sorted(names) == ['classAImpl', 'classBImpl']
