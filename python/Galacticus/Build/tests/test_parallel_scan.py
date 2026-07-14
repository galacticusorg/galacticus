"""Tests for Galacticus.Build.ParallelScan.

Pins the two contracts every cataloging script relies on:

* results come back in task order, identical to a serial run (the emitted
  makefiles/catalogs must be byte-identical whatever the parallelism);
* the worker count follows GALACTICUS_BUILD_JOBS, then make's ``-j`` from
  MAKEFLAGS (bare ``-j`` = unlimited, absent ``-j`` under make = serial),
  then the standalone default, clamped to the task count.
"""

import pytest

from Galacticus.Build.ParallelScan import resolve_jobs, scan, _DEFAULT_STANDALONE


# Workers must be module-level so the process pool can pickle them.

def _square(task):
    return task * task


def _fail_on_three(task):
    if task == 3:
        raise ValueError("worker exploded")
    return task


# ---------------------------------------------------------------------------
# scan(): ordering and serial equivalence
# ---------------------------------------------------------------------------

def test_results_in_task_order_parallel():
    tasks = list(range(50))
    assert scan(tasks, _square, 'test', jobs=8) == [t * t for t in tasks]


def test_parallel_matches_serial():
    tasks = list(range(23))
    serial   = scan(tasks, _square, 'test', jobs=1)
    parallel = scan(tasks, _square, 'test', jobs=4)
    assert parallel == serial


def test_empty_tasks():
    assert scan([], _square, 'test') == []


def test_single_task_runs_serially():
    assert scan([7], _square, 'test', jobs=16) == [49]


def test_worker_failure_exits_with_serial_hint(capsys):
    with pytest.raises(SystemExit) as exc:
        scan(list(range(8)), _fail_on_three, 'myScript.py', jobs=4)
    message = str(exc.value)
    assert 'myScript.py' in message
    assert 'GALACTICUS_BUILD_JOBS=1' in message


# ---------------------------------------------------------------------------
# resolve_jobs(): precedence and clamping
# ---------------------------------------------------------------------------

def test_override_wins(monkeypatch):
    monkeypatch.setenv('GALACTICUS_BUILD_JOBS', '3')
    monkeypatch.setenv('MAKEFLAGS', ' -j16')
    assert resolve_jobs(100) == 3


def test_override_zero_means_cpu_count(monkeypatch):
    monkeypatch.setenv('GALACTICUS_BUILD_JOBS', '0')
    assert resolve_jobs(10000) >= 1


def test_makeflags_jN(monkeypatch):
    monkeypatch.delenv('GALACTICUS_BUILD_JOBS', raising=False)
    monkeypatch.setenv('MAKEFLAGS', ' -j4 --jobserver-auth=fifo:/tmp/x')
    assert resolve_jobs(100) == 4


def test_makeflags_without_j_is_serial(monkeypatch):
    # Under make with no -j the build itself is serial; forking a pool of
    # workers on top of it would oversubscribe unasked.
    monkeypatch.delenv('GALACTICUS_BUILD_JOBS', raising=False)
    monkeypatch.setenv('MAKEFLAGS', '')
    assert resolve_jobs(100) == 1


def test_makeflags_bare_j_means_cpu_count(monkeypatch):
    monkeypatch.delenv('GALACTICUS_BUILD_JOBS', raising=False)
    monkeypatch.setenv('MAKEFLAGS', ' -j')
    assert resolve_jobs(100000) >= 1


def test_standalone_default(monkeypatch):
    monkeypatch.delenv('GALACTICUS_BUILD_JOBS', raising=False)
    monkeypatch.delenv('MAKEFLAGS', raising=False)
    assert resolve_jobs(10000) == _DEFAULT_STANDALONE


def test_clamped_to_task_count(monkeypatch):
    monkeypatch.setenv('GALACTICUS_BUILD_JOBS', '64')
    assert resolve_jobs(2) == 2


def test_at_least_one(monkeypatch):
    monkeypatch.setenv('GALACTICUS_BUILD_JOBS', '-5')
    assert resolve_jobs(10) >= 1
