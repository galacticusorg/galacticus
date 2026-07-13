"""Tests for scripts/build/sourceDigests.py dependency-file tracking.

Regression tests for the hierarchical-source-tree handling of
`_append_dep_files`: the `.d` sidecar and the object entries it lists both
carry sub-directory components (`$BUILDPATH/<subdir>/<name>.{d,o}`), which an
earlier version discarded (it kept only basenames), silently tracking no
dependencies for any file in a subdirectory — i.e. almost every file.
"""

import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir))
from sourceDigests import _append_dep_files, _object_files_from_dep  # noqa: E402


def _make_tree(tmp_path):
    """Create a minimal hierarchical source/build tree:

    source/interface/fortran_interop.F90 depends (per its .d sidecar) on
    itself, on source/utility/hashes/md5.F90, and on an include
    source/objects/nodes.Inc.
    """
    source_root = tmp_path / 'source'
    build_path  = tmp_path / 'build'
    (source_root / 'interface').mkdir(parents=True)
    (source_root / 'utility' / 'hashes').mkdir(parents=True)
    (source_root / 'objects').mkdir(parents=True)
    (build_path / 'interface').mkdir(parents=True)

    (source_root / 'interface' / 'fortran_interop.F90').write_text('program p\nend\n')
    (source_root / 'utility' / 'hashes' / 'md5.F90').write_text('module m\nend\n')
    (source_root / 'objects' / 'nodes.Inc').write_text('! include\n')

    dep_file = build_path / 'interface' / 'fortran_interop.d'
    dep_file.write_text(
        f"{build_path}/interface/fortran_interop.o\n"
        f"{build_path}/utility/hashes/md5.o\n"
        f"{build_path}/objects/nodes.o\n"
        f"{build_path}/no/such/file.o\n"
    )
    return str(source_root), str(build_path)


def test_object_files_from_dep_preserves_subdirectories(tmp_path):
    source_root, build_path = _make_tree(tmp_path)
    objs = _object_files_from_dep(
        os.path.join(build_path, 'interface', 'fortran_interop.d'), build_path)
    assert 'interface/fortran_interop.o' in objs
    assert 'utility/hashes/md5.o' in objs


def test_append_dep_files_hierarchical(tmp_path):
    source_root, build_path = _make_tree(tmp_path)
    files = []
    _append_dep_files(files, 'interface/fortran_interop.F90',
                      build_path, source_root)
    assert os.path.join(source_root, 'interface', 'fortran_interop.F90') in files
    assert os.path.join(source_root, 'utility', 'hashes', 'md5.F90') in files
    # `.Inc` sources are found via the suffix search.
    assert os.path.join(source_root, 'objects', 'nodes.Inc') in files
    # Objects with no corresponding source are skipped silently.
    assert not any('no/such' in f for f in files)


def test_append_dep_files_missing_sidecar_is_noop(tmp_path):
    source_root, build_path = _make_tree(tmp_path)
    files = []
    _append_dep_files(files, 'interface/no_sidecar.F90',
                      build_path, source_root)
    assert files == []


def test_append_dep_files_non_f90_is_noop(tmp_path):
    source_root, build_path = _make_tree(tmp_path)
    files = []
    _append_dep_files(files, 'interface/fortran_interop.cpp',
                      build_path, source_root)
    assert files == []
