"""Integration tests for the `++Task:` contract between profiler.sh
(producer, the SHELL wrapper of `GALACTICUS_BUILD_OPTION=compileprof`
builds) and buildProfiler.py (consumer, the report generator) — an
undocumented-by-code format pair where neither side knows the other exists.

Also pins profiler.sh's exit-status propagation: as make's SHELL it must
forward the wrapped command's status, or make cannot detect any failed
recipe during a profiled build.
"""

import os
import re
import subprocess
import sys

import pytest

_REPO      = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          os.pardir, os.pardir, os.pardir)
_PROFILER  = os.path.join(_REPO, 'scripts', 'build', 'profiler.sh')
_REPORTER  = os.path.join(_REPO, 'scripts', 'build', 'buildProfiler.py')

# The consumer's parse regex (mirrors buildProfiler.py).
_TASK_RE = re.compile(
    r"^\+\+Task: \{([\s\d\+\-:]+)\|([\s\d\+\-:]+)(?:\|(-?\d+))?\} '(.*)")

# profiler.sh uses GNU `date --rfc-3339=seconds`; skip where unsupported
# (e.g. BSD date on macOS).
_HAVE_GNU_DATE = subprocess.run(
    ['date', '--rfc-3339=seconds'], capture_output=True).returncode == 0

pytestmark = pytest.mark.skipif(
    not _HAVE_GNU_DATE, reason='requires GNU date (--rfc-3339)')


def _run(command):
    """Invoke profiler.sh the way make does: `$(SHELL) -c '<command>'`."""
    return subprocess.run([_PROFILER, '-c', command],
                          capture_output=True, text=True)


def test_task_line_parses_with_consumer_regex():
    result = _run('true')
    task_lines = [l for l in result.stdout.splitlines()
                  if l.startswith('++Task:')]
    assert len(task_lines) == 1
    m = _TASK_RE.match(task_lines[0])
    assert m, f'consumer regex rejected: {task_lines[0]!r}'
    assert m.group(4) == "true'"          # command, with closing quote
    assert int(m.group(3)) >= -1          # maxRSS in kB, or -1 sentinel


def test_exit_status_propagated():
    assert _run('true').returncode == 0
    assert _run('false').returncode == 1
    assert _run('exit 7').returncode == 7


def test_command_output_passes_through():
    result = _run('echo recipe-output')
    assert 'recipe-output' in result.stdout


def test_report_generated_from_task_lines(tmp_path):
    # End-to-end: a profiler-produced log renders to an HTML report.
    log = tmp_path / 'build.log'
    log.write_text(_run('sleep 0').stdout)
    html = tmp_path / 'profile.html'
    result = subprocess.run(
        [sys.executable, _REPORTER, str(log), str(html),
         '--durationMinimum', '0'],
        capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert html.exists() and html.stat().st_size > 0
