"""Tests for scripts/build/preprocessorSources.py.

Regression tests for issue #1395: the preprocessing rules in the Makefile
listed the Fortran input and the directive catalogs among their prerequisites,
but nothing about the Python program that does the preprocessing. Editing
`preprocess.py`, a `SourceTree/Process` hook, or any other module in the
import closure therefore regenerated no `.p.F90` at all — make reported
success and the stale generated Fortran was compiled.

The fix hangs those rules off a digest of the preprocessor's own sources, so
these tests pin the two properties the digest has to have: the closure is
complete (nothing the preprocessor imports is missed, or the staleness comes
straight back), and it is not the whole of `python/` (or every unrelated
Python edit re-preprocesses the tree).
"""

import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir))
from preprocessorSources import ENTRY_POINTS, closure, digest  # noqa: E402


REPO_ROOT = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, os.pardir))


def _write(path, text):
    os.makedirs(os.path.dirname(str(path)), exist_ok=True)
    with open(str(path), 'w', encoding='utf-8') as fh:
        fh.write(text)


def _relative(paths, root):
    return {os.path.relpath(path, str(root)) for path in paths}


def _tree(tmp_path):
    """A miniature `python/` tree plus an entry-point script outside it."""
    root = tmp_path / 'python'
    _write(root / 'Pkg' / '__init__.py',    "VERSION = 1\n")
    _write(root / 'Pkg' / 'Head.py',        "import Pkg.Tail\nimport numpy\n")
    _write(root / 'Pkg' / 'Tail.py',        "from Other import helper\n")
    _write(root / 'Pkg' / 'Unused.py',      "raise RuntimeError('never imported')\n")
    _write(root / 'Other' / '__init__.py',  "def helper():\n    pass\n")
    _write(root / 'Lazy.py',                "OK = True\n")
    _write(tmp_path / 'entry.py',           "from Pkg import Head\n\n"
                                            "def run():\n    import Lazy\n    return Lazy.OK\n")
    return [str(tmp_path / 'entry.py')], [str(root)]


def test_closure_follows_imports_transitively(tmp_path):
    entry_points, roots = _tree(tmp_path)
    found = _relative(closure(entry_points, roots), tmp_path)
    assert found == {
        'entry.py',
        os.path.join('python', 'Pkg', '__init__.py'),
        os.path.join('python', 'Pkg', 'Head.py'),
        os.path.join('python', 'Pkg', 'Tail.py'),
        os.path.join('python', 'Other', '__init__.py'),
        # Imported inside a function: the hooks do this to break import cycles,
        # and the module is every bit as much part of the preprocessor.
        os.path.join('python', 'Lazy.py'),
    }


def test_closure_excludes_unimported_and_third_party(tmp_path):
    entry_points, roots = _tree(tmp_path)
    found = _relative(closure(entry_points, roots), tmp_path)
    # `Pkg.Unused` is in the tree but nothing imports it; `numpy` is not ours.
    assert os.path.join('python', 'Pkg', 'Unused.py') not in found
    assert not any('numpy' in path for path in found)


def test_closure_resolves_relative_imports(tmp_path):
    root = tmp_path / 'python'
    _write(root / 'Pkg' / '__init__.py',      "from . import Sibling\n")
    _write(root / 'Pkg' / 'Sibling.py',       "from ..Neighbor import thing\n")
    _write(root / 'Neighbor' / '__init__.py', "thing = 1\n")
    _write(tmp_path / 'entry.py',             "import Pkg\n")
    found = _relative(closure([str(tmp_path / 'entry.py')], [str(root)]), tmp_path)
    assert os.path.join('python', 'Pkg', 'Sibling.py')       in found
    assert os.path.join('python', 'Neighbor', '__init__.py') in found


def test_closure_follows_sibling_scripts(tmp_path):
    """A script importing a sibling script in its own directory is followed."""
    _write(tmp_path / 'scripts' / 'entry.py',   "import sidekick\n")
    _write(tmp_path / 'scripts' / 'sidekick.py', "X = 1\n")
    found = _relative(closure([str(tmp_path / 'scripts' / 'entry.py')], []), tmp_path)
    assert found == {os.path.join('scripts', 'entry.py'),
                     os.path.join('scripts', 'sidekick.py')}


def test_digest_tracks_content_not_mtime(tmp_path):
    entry_points, roots = _tree(tmp_path)
    paths    = closure(entry_points, roots)
    original = digest(paths, str(tmp_path))

    # Touching a file without changing it must not change the digest: an
    # unchanged digest is what keeps a spurious re-preprocess sweep away.
    os.utime(str(tmp_path / 'python' / 'Pkg' / 'Tail.py'), (0, 0))
    assert digest(closure(entry_points, roots), str(tmp_path)) == original

    # Changing one does.
    _write(tmp_path / 'python' / 'Pkg' / 'Tail.py',
           "from Other import helper\n# altered\n")
    assert digest(closure(entry_points, roots), str(tmp_path)) != original


def test_digest_ignores_files_outside_the_closure(tmp_path):
    entry_points, roots = _tree(tmp_path)
    original = digest(closure(entry_points, roots), str(tmp_path))
    _write(tmp_path / 'python' / 'Pkg' / 'Unused.py', "raise RuntimeError('still unused')\n")
    assert digest(closure(entry_points, roots), str(tmp_path)) == original


def test_real_preprocessor_closure():
    """The closure of the real preprocessor, against the real tree."""
    entry_points = [os.path.join(REPO_ROOT, entry_point) for entry_point in ENTRY_POINTS]
    found        = _relative(closure(entry_points, [os.path.join(REPO_ROOT, 'python')]), REPO_ROOT)

    # The driver, the parser, and the hook registry -- editing any of these
    # changes the generated Fortran.
    for expected in [
        os.path.join('scripts', 'build', 'preprocess.py'),
        os.path.join('python', 'Galacticus', 'Build', 'SourceTree', '__init__.py'),
        os.path.join('python', 'Galacticus', 'Build', 'SourceTree', 'Process', 'all.py'),
        os.path.join('python', 'Fortran', 'Utils.py'),
    ]:
        assert expected in found

    # Every registered Process hook must be present: `Process/all.py` imports
    # the full set, and a hook missing from the closure is a hook whose edits
    # would silently not be applied to the already-preprocessed tree.
    process_directory = os.path.join('python', 'Galacticus', 'Build', 'SourceTree', 'Process')
    for file_name in os.listdir(os.path.join(REPO_ROOT, process_directory)):
        if file_name.endswith('.py'):
            assert os.path.join(process_directory, file_name) in found

    # The node-component generator is reached only through a function-level
    # import in the componentBuilder hook, but it emits Fortran all the same.
    assert os.path.join('python', 'Galacticus', 'Build', 'Components', 'Components.py') in found

    # ... and the closure is genuinely narrower than "all of python/": Python
    # that has nothing to do with preprocessing must not drag the whole tree
    # through the preprocessor again.
    assert os.path.join('python', 'cloudy.py')                                          not in found
    assert os.path.join('python', 'LibraryInterfaces', 'Emitters.py')                   not in found
    assert os.path.join('python', 'Galacticus', 'Parameters', 'catalog.py')             not in found
