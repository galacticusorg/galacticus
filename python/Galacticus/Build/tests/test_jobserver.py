"""Tests for Galacticus.Build.Jobserver.

The contract that matters to the build: never take more slots than make has
free, and never lose a token. A leaked token permanently shrinks the whole
build's parallelism, so the release path is pinned here for both the `fifo:`
(make >= 4.4) and legacy `R,W` descriptor transports, including when the body
raises.
"""

import os
import subprocess
import sys
import textwrap

import pytest

from Galacticus.Build.Jobserver import slots, jobserver_available


@pytest.fixture
def fifo_jobserver(tmp_path, monkeypatch):
    """A `fifo:`-style jobserver holding `n` tokens, advertised via MAKEFLAGS."""
    path = str(tmp_path / 'fifo')
    os.mkfifo(path)
    fd = os.open(path, os.O_RDWR | os.O_NONBLOCK)

    def fill(n):
        os.write(fd, b'+' * n)
        monkeypatch.setenv('MAKEFLAGS', f' -j8 --jobserver-auth=fifo:{path}')
        return fd

    yield fill
    os.close(fd)


@pytest.fixture
def fds_jobserver(monkeypatch):
    """A legacy `R,W` descriptor-pair jobserver holding `n` tokens."""
    read_fd, write_fd = os.pipe()

    def fill(n):
        os.write(write_fd, b'+' * n)
        monkeypatch.setenv('MAKEFLAGS',
                           f' -j8 --jobserver-auth={read_fd},{write_fd}')
        return read_fd, write_fd

    yield fill
    os.close(read_fd)
    os.close(write_fd)


# ---------------------------------------------------------------------------
# No jobserver -> caller's own policy is untouched
# ---------------------------------------------------------------------------

def test_no_makeflags_yields_desired(monkeypatch):
    monkeypatch.delenv('MAKEFLAGS', raising=False)
    with slots(8) as granted:
        assert granted == 8


def test_makeflags_without_auth_yields_desired(monkeypatch):
    monkeypatch.setenv('MAKEFLAGS', ' -j8')
    with slots(8) as granted:
        assert granted == 8
    assert not jobserver_available()


def test_unreachable_fifo_yields_desired(monkeypatch, tmp_path):
    monkeypatch.setenv('MAKEFLAGS',
                       f" -j8 --jobserver-auth=fifo:{tmp_path / 'absent'}")
    with slots(8) as granted:
        assert granted == 8


def test_stale_fds_yield_desired(monkeypatch):
    # Descriptors make did not pass down (recipe not marked `+`).
    monkeypatch.setenv('MAKEFLAGS', ' -j8 --jobserver-auth=71,72')
    with slots(8) as granted:
        assert granted == 8


# ---------------------------------------------------------------------------
# Token accounting
# ---------------------------------------------------------------------------

def test_fifo_grants_tokens_plus_implicit_slot(fifo_jobserver):
    fifo_jobserver(3)
    with slots(8) as granted:
        # 3 tokens + the slot we implicitly own.
        assert granted == 4


def test_fifo_never_exceeds_desired(fifo_jobserver):
    fifo_jobserver(32)
    with slots(4) as granted:
        assert granted == 4


def test_fifo_empty_pool_falls_back_to_serial(fifo_jobserver):
    fifo_jobserver(0)
    with slots(8) as granted:
        assert granted == 1


def test_fds_transport_grants_tokens(fds_jobserver):
    fds_jobserver(2)
    with slots(8) as granted:
        assert granted == 3


def test_tokens_returned_to_pool(fifo_jobserver):
    fd = fifo_jobserver(3)
    with slots(8) as granted:
        assert granted == 4
    # All three must be back, or the build loses them for good.
    assert os.read(fd, 16) == b'+++'


def test_tokens_returned_on_exception(fifo_jobserver):
    fd = fifo_jobserver(3)
    with pytest.raises(RuntimeError):
        with slots(8):
            raise RuntimeError("boom")
    assert os.read(fd, 16) == b'+++'


def test_token_bytes_preserved_exactly(fifo_jobserver):
    """make distinguishes its tokens by value; they must go back unaltered."""
    fd = fifo_jobserver(0)
    os.write(fd, b'abc')
    with slots(4) as granted:
        assert granted == 4      # 3 distinct tokens + our implicit slot
    assert sorted(os.read(fd, 16)) == sorted(b'abc')


def test_desired_one_takes_no_tokens(fifo_jobserver):
    fd = fifo_jobserver(2)
    with slots(1) as granted:
        assert granted == 1
    # Nothing was taken, so the pool is untouched.
    assert os.read(fd, 16) == b'++'


def test_jobserver_available(fifo_jobserver):
    fifo_jobserver(1)
    assert jobserver_available()


# ---------------------------------------------------------------------------
# End-to-end under a real `make`
# ---------------------------------------------------------------------------

def _run_probe(tmp_path, make_args):
    """Run a `+`-marked recipe that asks `slots` for 64 workers under make."""
    package_root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))))
    probe = tmp_path / 'probe.py'
    probe.write_text(textwrap.dedent(f"""
        import sys
        sys.path.insert(0, {package_root!r})
        from Galacticus.Build.Jobserver import slots
        with slots(64) as granted:
            print("GRANTED", granted)
    """))
    makefile = tmp_path / 'Makefile'
    makefile.write_text(f"probe:\n\t+@{sys.executable} {probe}\n")
    result = subprocess.run(['make', '-f', str(makefile), *make_args, 'probe'],
                            capture_output=True, text=True, cwd=str(tmp_path))
    assert result.returncode == 0, result.stderr
    return int(result.stdout.split('GRANTED')[1].split()[0])


@pytest.mark.parametrize("jobs", [2, 4])
def test_under_real_make(tmp_path, jobs):
    """A `+`-marked recipe asking for far more than make allows is capped to it.

    This is the whole point of the module, so it is pinned against a real make
    rather than a synthetic pipe: make -jN publishes N-1 tokens, and with our
    probe as the only running job we should get exactly those plus our own slot.
    """
    if subprocess.run(['make', '--version'], capture_output=True).returncode != 0:
        pytest.skip("make unavailable")
    assert _run_probe(tmp_path, [f'-j{jobs}']) == jobs


def test_under_real_make_j1_has_no_jobserver(tmp_path):
    """`make -j1` is serial and publishes no jobserver, so `slots` cannot cap.

    Nothing is over-subscribed in practice: `ParallelScan.resolve_jobs` reads the
    `-j1` from MAKEFLAGS and asks for a single worker in the first place. Pinned
    so the division of labour between the two stays deliberate.
    """
    if subprocess.run(['make', '--version'], capture_output=True).returncode != 0:
        pytest.skip("make unavailable")
    assert _run_probe(tmp_path, ['-j1']) == 64
