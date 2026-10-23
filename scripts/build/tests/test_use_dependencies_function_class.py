"""Tests for useDependencies.py's tracking of functionClass instance files.

A file declaring a `functionClass` has the `use` statements of every instance
(implementation) file attributed to it. The scan cache decides whether to
rescan a file from the mtimes of the files its entry tracks, so the instance
files must be tracked too, or an edit to a `use` in an existing instance leaves
stale dependencies in Makefile_Use_Dependencies. A change in which files
implement the class (an existing file gaining or losing an instance directive,
which no mtime of a tracked file reveals) must likewise trigger a rescan. The
same applies to the files containing `functionGlobal` directives (whose modules
`functionsGlobal` users import) and to the files declaring static event hooks
(whose containing modules `eventHookStatic` users import).
"""

import os
import subprocess
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
    assert entry['directiveFiles'] == {'widget': ['/src/a.F90']}
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
    assert entry['directiveFiles'] == {
        'widget': ['/src/a.F90', '/src/b.F90'],
        'gadget': ['/src/c.F90'],
    }


def test_no_locations():
    entry, to_process = _apply(None, ['widget'])
    assert to_process == ['/src/_class.F90']
    assert entry['files'] == ['/src/_class.F90']
    assert entry['directiveFiles'] == {'widget': []}


def test_entry_without_directive_files_never_changed():
    assert not ud._directive_files_changed(_new_entry('/src/x.F90'), {})


def test_directive_files_unchanged():
    locations = {'widget': {'file': ['/src/a.F90', '/src/b.F90']}}
    entry, _ = _apply(locations, ['widget'])
    # Order in directiveLocations.xml may vary between regenerations.
    reordered = {'widget': {'file': ['/src/b.F90', '/src/a.F90']}}
    assert not ud._directive_files_changed(entry, reordered)


def test_directive_file_added():
    entry, _ = _apply({'widget': {'file': '/src/a.F90'}}, ['widget'])
    added = {'widget': {'file': ['/src/a.F90', '/src/b.F90']}}
    assert ud._directive_files_changed(entry, added)


def test_directive_file_removed():
    locations = {'widget': {'file': ['/src/a.F90', '/src/b.F90']}}
    entry, _ = _apply(locations, ['widget'])
    assert ud._directive_files_changed(
        entry, {'widget': {'file': '/src/a.F90'}})


def test_directive_file_change_in_second_class_detected():
    locations = {'widget': {'file': '/src/a.F90'},
                 'gadget': {'file': '/src/c.F90'}}
    entry, _ = _apply(locations, ['widget', 'gadget'])
    locations['gadget'] = {'file': ['/src/c.F90', '/src/d.F90']}
    assert ud._directive_files_changed(entry, locations)


# ---------------------------------------------------------------------------
# End-to-end: incremental runs of the script itself
# ---------------------------------------------------------------------------

_SCRIPT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, 'useDependencies.py')
# The repository's python/ directory, which the Makefile puts on PYTHONPATH for
# the build scripts (it is on sys.path here only via pytest, not in a child).
_PYTHON_DIR = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, os.pardir,
    'python'))

_CLASS_SOURCE = """module Widgets
  !![
  <functionClass>
   <name>widget</name>
   <descriptiveName>Widget</descriptiveName>
   <description>A widget.</description>
   <default>a</default>
  </functionClass>
  !!]
end module Widgets
"""

_POINTERS_SOURCE = """module Function_Global_Pointers
  !![
  <functionsGlobal type="pointers"/>
  !!]
end module Function_Global_Pointers
"""


def _function_global_source(module):
    return f"""module Provider
  !![
  <functionGlobal>
   <unitName>Provider_Function</unitName>
   <type>void</type>
   <module>{module}, only : x</module>
  </functionGlobal>
  !!]
end module Provider
"""


class _Tree:
    """A minimal source tree and build directory on which useDependencies.py
    can be run repeatedly, as make would run it."""

    def __init__(self, tmp_path):
        self.root  = tmp_path
        self.build = tmp_path / 'work' / 'build'
        self.build.mkdir(parents=True)
        (tmp_path / 'source').mkdir()
        self.locations = {}
        self.clock     = 1_000_000_000

    def write(self, name, text):
        """Write a source file, with an mtime later than any earlier run."""
        path = self.root / 'source' / name
        path.write_text(text)
        self.clock += 10
        os.utime(path, (self.clock, self.clock))
        return str(path)

    def run(self):
        xml = ['<directives>']
        for name, files in self.locations.items():
            xml.append(f'  <{name}>')
            xml.extend(f'    <file>{f}</file>' for f in files)
            xml.append(f'  </{name}>')
        xml.append('</directives>')
        (self.build / 'directiveLocations.xml').write_text('\n'.join(xml))
        env = dict(os.environ, BUILDPATH=str(self.build))
        env['PYTHONPATH'] = os.pathsep.join(
            [_PYTHON_DIR] + ([env['PYTHONPATH']] if env.get('PYTHONPATH') else []))
        env.pop('GALACTICUS_FCFLAGS', None)
        result = subprocess.run(
            [sys.executable, _SCRIPT, str(self.root)], cwd=str(self.root),
            env=env, capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        # The script stamps the blob when it writes it; move that stamp
        # past every source file written so far, and before any written next.
        self.clock += 10
        os.utime(self.build / 'Makefile_Use_Dependencies.blob',
                 (self.clock, self.clock))

    def modules(self, name):
        """The .mod files in the Makefile rule for source/<name>'s object."""
        target = str(self.build) + '/' + name.replace('.F90', '.o') + ':'
        for line in (self.build / 'Makefile_Use_Dependencies').read_text() \
                .splitlines():
            if line.startswith(target):
                return {os.path.basename(d) for d in line.split()[1:]
                        if d.endswith('.mod')}
        return set()


def test_main_rescans_class_when_instance_edited(tmp_path):
    tree = _Tree(tmp_path)
    tree.write('widget.F90', _CLASS_SOURCE)
    a = tree.write('a.F90', '  use :: Module_Old\n')
    tree.locations = {'widget': [a]}
    tree.run()
    assert 'module_old.mod' in tree.modules('widget.F90')
    tree.write('a.F90', '  use :: Module_New\n')
    tree.run()
    assert 'module_new.mod' in tree.modules('widget.F90')
    assert 'module_old.mod' not in tree.modules('widget.F90')


def test_main_rescans_class_when_existing_file_gains_instance(tmp_path):
    tree = _Tree(tmp_path)
    tree.write('widget.F90', _CLASS_SOURCE)
    a = tree.write('a.F90', '  use :: Module_A\n')
    b = tree.write('b.F90', '  use :: Module_B\n')
    tree.locations = {'widget': [a]}
    tree.run()
    assert 'module_b.mod' not in tree.modules('widget.F90')
    # b.F90 becomes an instance without the class file or `a.F90` changing
    # (in reality its new directive would change its own mtime, but b.F90 is
    # not tracked by the class entry).
    tree.locations = {'widget': [a, b]}
    tree.run()
    assert 'module_b.mod' in tree.modules('widget.F90')


def test_main_rescans_pointers_when_function_global_module_edited(tmp_path):
    tree = _Tree(tmp_path)
    tree.write('pointers.F90', _POINTERS_SOURCE)
    p = tree.write('provider.F90', _function_global_source('Module_Old'))
    tree.locations = {'functionGlobal': [p]}
    tree.run()
    assert 'module_old.mod' in tree.modules('pointers.F90')
    tree.write('provider.F90', _function_global_source('Module_New'))
    tree.run()
    assert 'module_new.mod' in tree.modules('pointers.F90')
    assert 'module_old.mod' not in tree.modules('pointers.F90')


def test_main_rescans_pointers_when_existing_file_gains_function_global(
        tmp_path):
    tree = _Tree(tmp_path)
    tree.write('pointers.F90', _POINTERS_SOURCE)
    tree.write('provider.F90', 'module Provider\nend module Provider\n')
    tree.locations = {}
    tree.run()
    assert 'module_new.mod' not in tree.modules('pointers.F90')
    p = tree.write('provider.F90', _function_global_source('Module_New'))
    tree.locations = {'functionGlobal': [p]}
    tree.run()
    assert 'module_new.mod' in tree.modules('pointers.F90')


_HOOK_USER_SOURCE = """module Hook_User
  !![
  <eventHookStatic name="widgetHook">
  </eventHookStatic>
  !!]
end module Hook_User
"""


def test_main_rescans_static_hook_user_when_hook_moves(tmp_path):
    tree = _Tree(tmp_path)
    tree.write('user.F90', _HOOK_USER_SOURCE)
    hook = """module {}
  !![
  <widgetHook>
  </widgetHook>
  !!]
end module {}
"""
    one = tree.write('one.F90', hook.format('Hook_One', 'Hook_One'))
    two = tree.write('two.F90', 'module Hook_Two\nend module Hook_Two\n')
    tree.locations = {'widgetHook': [one]}
    tree.run()
    assert 'hook_one.mod' in tree.modules('user.F90')
    # The hook moves to another existing file.
    tree.write('one.F90', 'module Hook_One\nend module Hook_One\n')
    tree.write('two.F90', hook.format('Hook_Two', 'Hook_Two'))
    tree.locations = {'widgetHook': [two]}
    tree.run()
    assert 'hook_two.mod' in tree.modules('user.F90')
    assert 'hook_one.mod' not in tree.modules('user.F90')
