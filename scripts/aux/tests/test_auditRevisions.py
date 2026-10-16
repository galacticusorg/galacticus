"""Tests for scripts/aux/auditRevisions.py.

The audit fails the `Validate-Revisions` PR check when a parameter file's
`lastModified` revision, or a commit registered in `migrations.xml`, is not in
the history of `HEAD`. Both kinds of hash are typically lost when the branch
which recorded them is rebased or amended before it is merged, and each breaks
migration: a missing `lastModified` revision makes `parametersMigrate.py`
refuse to run, and a missing migration commit makes that migration silently
never fire. The cases below build a small throwaway repository and check that

  * a clean repository passes,
  * a `lastModified` revision which does not exist, or which exists only on a
    branch not merged into `HEAD`, is reported,
  * a migration commit which is abbreviated, does not exist, or is not an
    ancestor of `HEAD` is reported, and
  * a shallow clone is refused rather than reported as full of missing commits.
"""

import importlib.util
import os
import subprocess

import pytest

_AUDIT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, 'auditRevisions.py')


def _load():
    """Import `auditRevisions.py` as a module (it has no `.py`-importable name
    on the path, and executes nothing at import time beyond its definitions)."""
    spec = importlib.util.spec_from_file_location('auditRevisions', _AUDIT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


audit = _load()


def _git(root, *arguments):
    # An empty template and no hooks, so that the throwaway repository is unaffected by any local git configuration.
    return subprocess.run(
        ['git', '-C', str(root), '-c', 'core.hooksPath=/dev/null', '-c', 'user.name=Test', '-c', 'user.email=test@example.com',
         '-c', 'commit.gpgSign=false', *arguments],
        capture_output=True, text=True, check=True).stdout.strip()


def _write(root, path, content):
    full = os.path.join(root, path)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, 'w') as file:
        file.write(content)


def _parameters(revision):
    return f'<parameters>\n  <lastModified revision="{revision}" time="2026-01-01T00:00:00"/>\n</parameters>\n'


def _migrations(*commits):
    entries = ''.join(f'  <migration commit="{commit}"/>\n' for commit in commits)
    return f'<migrations>\n{entries}</migrations>\n'


@pytest.fixture
def repository(tmp_path):
    """A repository with two commits on its main line; returns (root, first, second)."""
    root = tmp_path / 'repo'
    root.mkdir()
    _git(root, 'init', '--quiet', '--template=', '-b', 'main')
    _write(root, 'README', 'one\n')
    _git(root, 'add', '-A')
    _git(root, 'commit', '--quiet', '-m', 'first')
    first = _git(root, 'rev-parse', 'HEAD')
    _write(root, 'README', 'two\n')
    _git(root, 'add', '-A')
    _git(root, 'commit', '--quiet', '-m', 'second')
    second = _git(root, 'rev-parse', 'HEAD')
    return root, first, second


def _commit(root, files):
    for path, content in files.items():
        _write(root, path, content)
    _git(root, 'add', '-A')
    _git(root, 'commit', '--quiet', '-m', 'parameters')


def test_clean(repository):
    root, first, second = repository
    _commit(root, {'parameters/a.xml': _parameters(first), 'scripts/aux/migrations.xml': _migrations(second)})
    assert audit.audit(str(root)) == []


def test_missing_last_modified(repository):
    root, first, second = repository
    missing = 'f' * 40
    _commit(root, {
        'parameters/a.xml': _parameters(missing),
        'parameters/b.xml': _parameters(missing),
        'parameters/c.xml': _parameters(first),
        'scripts/aux/migrations.xml': _migrations(second),
    })
    issues = audit.audit(str(root))
    assert len(issues) == 1
    assert missing in issues[0] and 'does not resolve' in issues[0] and '2 file(s)' in issues[0]


def test_last_modified_attributes_split_across_lines(repository):
    root, first, second = repository
    missing = 'e' * 40
    content = f'<parameters>\n  <lastModified time="2026-01-01T00:00:00"\n                revision="{missing}"/>\n</parameters>\n'
    _commit(root, {'parameters/a.xml': content, 'scripts/aux/migrations.xml': _migrations(second)})
    issues = audit.audit(str(root))
    assert len(issues) == 1 and missing in issues[0]


def test_unmerged_branch(repository):
    root, first, second = repository
    # A commit which exists, but only on a branch which HEAD does not contain.
    _git(root, 'checkout', '--quiet', '-b', 'side', first)
    _write(root, 'README', 'side\n')
    _git(root, 'add', '-A')
    _git(root, 'commit', '--quiet', '-m', 'side')
    side = _git(root, 'rev-parse', 'HEAD')
    _git(root, 'checkout', '--quiet', 'main')
    _commit(root, {'parameters/a.xml': _parameters(side), 'scripts/aux/migrations.xml': _migrations(side)})
    issues = audit.audit(str(root))
    assert len(issues) == 2
    assert all('is not an ancestor of HEAD' in issue for issue in issues)


def test_bad_migration_commits(repository):
    root, first, second = repository
    _commit(root, {
        'parameters/a.xml': _parameters(first),
        'scripts/aux/migrations.xml': _migrations(second[:10], second.upper(), 'd' * 40, second),
    })
    issues = audit.audit(str(root))
    assert len(issues) == 3
    assert 'full, lower-case' in issues[0] and 'full, lower-case' in issues[1]
    assert 'd' * 40 in issues[2] and 'does not resolve' in issues[2]


def test_shallow_clone_refused(repository, tmp_path, monkeypatch, capsys):
    root, first, second = repository
    _commit(root, {'parameters/a.xml': _parameters(first), 'scripts/aux/migrations.xml': _migrations(second)})
    shallow = tmp_path / 'shallow'
    subprocess.run(['git', 'clone', '--quiet', '--template=', '--depth', '1', f'file://{root}', str(shallow)],
                   capture_output=True, check=True)
    monkeypatch.setattr('sys.argv', ['auditRevisions.py', '--check', str(shallow)])
    assert audit.main() == 2
    assert 'shallow clone' in capsys.readouterr().err
