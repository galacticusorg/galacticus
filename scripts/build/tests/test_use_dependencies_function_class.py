"""Tests for useDependencies.py's tracking of functionClass instance files.

A file declaring a `functionClass` has the `use` statements of every instance
(implementation) file attributed to it. The scan cache decides whether to
rescan a file from the mtimes of the files its entry tracks, so the instance
files must be tracked too, or an edit to a `use` in an existing instance leaves
stale dependencies in Makefile_Use_Dependencies. A change in which files
implement the class (an existing file gaining or losing an instance directive,
which no mtime of a tracked file reveals) must likewise trigger a rescan.
"""

import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir))
import useDependencies as ud  # noqa: E402


def _new_entry(path):
    return {
        'files':                [path],
        'modulesUsed':          [],
        'dependenciesExplicit': [],
        'modulesProvided':      {},
        'submodules':           [],
        'submodulesProvided':   [],
        'libraryDependencies':  {},
    }


def _apply(locations, class_names):
    """Run _apply_directive_requirements for a file declaring the named
    functionClasses; return (entry, file_names_to_process)."""
    entry      = _new_entry('/src/_class.F90')
    directives = {name: [] for name in ud._DIRECTIVE_NAMES}
    directives['functionClass'] = [{'name': name} for name in class_names]
    to_process = ['/src/_class.F90']
    ud._apply_directive_requirements(
        entry, {'fileName': '_class.F90', 'subDirectoryName': ''},
        directives, locations, None, './work/build/', 'id', [], {},
        to_process,
    )
    return entry, to_process


def test_single_instance_tracked():
    # xml_to_dict gives a scalar `file` for a directive found in one file.
    locations = {'widget': {'file': '/src/a.F90'}}
    entry, to_process = _apply(locations, ['widget'])
    assert to_process == ['/src/_class.F90', '/src/a.F90']
    assert entry['files'] == ['/src/_class.F90', '/src/a.F90']
    assert entry['functionClasses'] == ['widget']
    assert entry['functionClassInstances'] == ['/src/a.F90']
    assert './work/build/function_classes.mod' in entry['modulesUsed']


def test_multiple_instances_and_classes_tracked():
    locations = {
        'widget': {'file': ['/src/b.F90', '/src/a.F90']},
        'gadget': {'file': '/src/c.F90'},
    }
    entry, to_process = _apply(locations, ['widget', 'gadget'])
    # Scan order follows directiveLocations.xml, as before.
    assert to_process[1:] == ['/src/b.F90', '/src/a.F90', '/src/c.F90']
    assert entry['files'][1:] == ['/src/b.F90', '/src/a.F90', '/src/c.F90']
    assert entry['functionClassInstances'] == [
        '/src/a.F90', '/src/b.F90', '/src/c.F90',
    ]


def test_no_locations():
    entry, to_process = _apply(None, ['widget'])
    assert to_process == ['/src/_class.F90']
    assert entry['files'] == ['/src/_class.F90']
    assert entry['functionClassInstances'] == []


def test_non_function_class_entry_never_changed():
    assert not ud._instances_changed(_new_entry('/src/x.F90'), {})


def test_instances_unchanged():
    locations = {'widget': {'file': ['/src/a.F90', '/src/b.F90']}}
    entry, _ = _apply(locations, ['widget'])
    # Order in directiveLocations.xml may vary between regenerations.
    reordered = {'widget': {'file': ['/src/b.F90', '/src/a.F90']}}
    assert not ud._instances_changed(entry, reordered)


def test_instance_added():
    entry, _ = _apply({'widget': {'file': '/src/a.F90'}}, ['widget'])
    added = {'widget': {'file': ['/src/a.F90', '/src/b.F90']}}
    assert ud._instances_changed(entry, added)


def test_instance_removed():
    locations = {'widget': {'file': ['/src/a.F90', '/src/b.F90']}}
    entry, _ = _apply(locations, ['widget'])
    assert ud._instances_changed(entry, {'widget': {'file': '/src/a.F90'}})
