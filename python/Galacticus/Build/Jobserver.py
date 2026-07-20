"""Client for GNU make's jobserver, used to size a build script's own pool.

`MAKEFLAGS` carries make's `-jN`, but N is the limit for the build *as a whole*:
every recipe make is running concurrently shares it. A script that reads `-j8`
and forks eight workers while make is already compiling seven files
oversubscribes the machine eightfold. The jobserver exists to arbitrate this --
make publishes a fixed pool of tokens, and a job may only run as many parallel
tasks as it holds tokens for.

This module implements the client side of that protocol (documented in the GNU
make manual, "Sharing Job Slots with GNU make"):

* make advertises the pool in `MAKEFLAGS` as `--jobserver-auth=<spec>`
  (`--jobserver-fds=<spec>` before make 4.2);
* `<spec>` is either `fifo:PATH` -- a named pipe, the default since make 4.4 --
  or `R,W`, a pair of inherited pipe file descriptors;
* the pool holds N-1 tokens; every job implicitly owns one slot for itself, so
  it reads a byte per *additional* task it wants to run in parallel;
* a read that would block means the pool is empty: someone else holds the
  tokens, and we run with fewer workers rather than waiting;
* every token read MUST be written back, byte-for-byte -- the bytes differ and
  make checks them. Losing one shrinks the build's parallelism permanently.

Tokens are taken once, up front, and held for the scan's duration. That is
coarser than acquiring per task, but a scan is a single short burst of work and
it keeps the failure mode safe: we can only ever under-subscribe, never over.

Whether a recipe can reach the pool depends on the transport, so recipes running
a scan should be marked `+` in the Makefile:

* with the `R,W` transport make only passes the descriptors down to `+`-marked
  recipes, so an unmarked recipe finds them closed and gets no pool (reported
  here as simply having no jobserver, leaving the caller's own `-j` policy);
* with `fifo:` the path is in `MAKEFLAGS` regardless, so an unmarked recipe can
  still reach the pool. Honouring it is the conservative outcome either way --
  we only ever take tokens that are genuinely free, and always give them back --
  but `+` is what the manual requires and the only portable way to get the
  descriptors.

Andrew Benson (2026).
"""

import contextlib
import errno
import fcntl
import os
import re
import select

__all__ = ['slots', 'jobserver_available']

_AUTH_RE = re.compile(r'--jobserver-(?:auth|fds)=(\S+)')


def _spec():
    """Return the `--jobserver-auth` spec from `MAKEFLAGS`, or None."""
    makeflags = os.environ.get('MAKEFLAGS')
    if not makeflags:
        return None
    match = _AUTH_RE.search(makeflags)
    return match.group(1) if match else None


def _open_fifo(path):
    """Open make's named-pipe jobserver.

    Opened O_RDWR: a read-only open would see EOF whenever no writer happens to
    hold the pipe, and O_NONBLOCK keeps an empty pool from blocking us. The file
    description is ours alone, so setting O_NONBLOCK cannot disturb make.
    """
    try:
        fd = os.open(path, os.O_RDWR | os.O_NONBLOCK)
    except OSError:
        return None
    return fd, fd, [fd]


def _open_fds(spec):
    """Adopt make's inherited pipe file descriptors (`R,W`, pre-4.4 style)."""
    try:
        read_text, write_text = spec.split(',', 1)
        read_fd, write_fd = int(read_text), int(write_text)
    except ValueError:
        return None
    for fd in (read_fd, write_fd):
        try:
            fcntl.fcntl(fd, fcntl.F_GETFD)
        except OSError:
            # Not inherited -- the recipe is not marked `+`, so make deliberately
            # withheld the pool.
            return None
    # These descriptors are shared with make: do NOT set O_NONBLOCK on them, as
    # the flag lives on the open file description and would leak into make
    # itself. `_take` polls instead.
    return read_fd, write_fd, []


def _connect():
    """Return `(read_fd, write_fd, fds_to_close)` for make's jobserver, or None."""
    spec = _spec()
    if spec is None:
        return None
    if spec.startswith('fifo:'):
        return _open_fifo(spec[len('fifo:'):])
    return _open_fds(spec)


def _take(read_fd, own_fd):
    """Read one token, or return None if the pool is empty.

    When the descriptor is ours (`fifo:`) it is already non-blocking. When it is
    make's own we must not touch its flags, so readiness is polled first; the
    poll/read gap is a benign race -- another job may take the token between the
    two, and we simply get fewer workers.
    """
    if not own_fd:
        if not select.select([read_fd], [], [], 0)[0]:
            return None
    try:
        token = os.read(read_fd, 1)
    except OSError as exc:
        if exc.errno in (errno.EAGAIN, errno.EWOULDBLOCK, errno.EBADF):
            return None
        raise
    return token or None


def jobserver_available():
    """True when a usable make jobserver is reachable (diagnostics/tests)."""
    connection = _connect()
    if connection is None:
        return False
    for fd in connection[2]:
        os.close(fd)
    return True


@contextlib.contextmanager
def slots(desired):
    """Yield the number of tasks we may actually run in parallel, <= `desired`.

    With no jobserver, `desired` is yielded unchanged (the caller's own policy
    applies). Otherwise up to `desired - 1` tokens are acquired -- the caller
    always owns one slot implicitly -- and released on exit, including on error.
    """
    connection = _connect()
    if connection is None or desired <= 1:
        yield desired
        return

    read_fd, write_fd, close_fds = connection
    own_fd = bool(close_fds)
    tokens = bytearray()
    try:
        while len(tokens) < desired - 1:
            token = _take(read_fd, own_fd)
            if token is None:
                break
            tokens += token
        yield 1 + len(tokens)
    finally:
        # Returning tokens is not optional: a dropped token is a slot the whole
        # build loses until make exits. Push back byte-for-byte and keep going
        # even if one write fails.
        for index in range(len(tokens)):
            try:
                os.write(write_fd, tokens[index:index + 1])
            except OSError:
                pass
        for fd in close_fds:
            try:
                os.close(fd)
            except OSError:
                pass
