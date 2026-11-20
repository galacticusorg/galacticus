"""Tests that the documentation and test-marker audits examine the tree in which they live.

`docBlockAudit.py` and `auditTestMarkers.py` are run before committing, often from a worktree in a shell whose
`GALACTICUS_EXEC_PATH` names a different checkout. If they audited the tree named by that variable, they would report success
about code they never opened. Each is run here with `GALACTICUS_EXEC_PATH` (and the working directory) pointing at an empty
decoy tree, and must still audit the files of its own tree.

Andrew Benson, Claude (19-November-2026).
"""

import os
import re
import subprocess
import sys

_AUX = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))


def _decoy(tmp_path):
    """Build an empty tree with the few files the audits read, and return the environment pointing at it."""
    for directory in ('source', os.path.join('docs', 'manuals'), 'testSuite'):
        os.makedirs(os.path.join(tmp_path, directory), exist_ok=True)
    for path in (os.path.join('docs', 'conf.py'), os.path.join('docs', 'Galacticus.bib')):
        open(os.path.join(tmp_path, path), 'w').close()
    environment = dict(os.environ)
    environment['GALACTICUS_EXEC_PATH'] = str(tmp_path)
    return environment


def test_docBlockAudit_audits_own_tree(tmp_path):
    environment = _decoy(tmp_path)
    result = subprocess.run([sys.executable, os.path.join(_AUX, 'docBlockAudit.py')], cwd=tmp_path, env=environment,
                            capture_output=True, text=True)
    match = re.search(r'in (\d+) source files', result.stdout)
    assert match, result.stdout + result.stderr
    # The decoy tree holds no source files; the tree holding the script holds well over a thousand.
    assert int(match.group(1)) > 1000, result.stdout


def test_auditTestMarkers_audits_own_tree(tmp_path):
    environment = _decoy(tmp_path)
    result = subprocess.run([sys.executable, os.path.join(_AUX, 'auditTestMarkers.py')], cwd=tmp_path, env=environment,
                            capture_output=True, text=True)
    match = re.search(r'(\d+) scripts audited', result.stdout)
    assert match, result.stdout + result.stderr
    # The decoy tree holds no test scripts; the tree holding the script holds over a hundred.
    assert int(match.group(1)) > 100, result.stdout
