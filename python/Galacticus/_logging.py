"""Default logging configuration for Galacticus library modules.

Convention used across the build pipeline:

* Library modules (under ``Galacticus.Build.*`` and ``Galacticus.Constraints.*``)
  declare a module-level logger via::

      import logging
      logger = logging.getLogger(__name__)

  and emit ``logger.info()`` / ``logger.warning()`` / ``logger.error()`` /
  ``logger.debug()`` for diagnostic output, in place of bare ``print()``.

* Script entrypoints under ``scripts/build/`` and ``scripts/aux/`` that
  drive the library code call :func:`configure_default` once at startup.
  Without that call, Python's root logger sits at ``WARNING`` and every
  ``logger.info`` from the library is silently dropped -- preserving the
  current visible-by-default behaviour requires this one-line opt-in at
  the entry point.

* Tests need not configure anything.  Pytest's ``caplog`` fixture captures
  log records regardless of the root configuration, and the default
  ``WARNING`` level means tests run quietly without the ``-->`` build
  progress output that ``print`` used to produce.
"""
from __future__ import annotations

import logging
import sys

__all__ = ['configure_default']

_DEFAULT_FORMAT = "%(message)s"


def configure_default(level: int = logging.INFO) -> None:
    """Set up a single stdout handler at *level* using a message-only
    format that mirrors the historical ``print()`` output style.

    Idempotent: subsequent calls are no-ops, which lets multiple
    entry-point scripts call it without doubling handlers when they
    invoke each other.
    """
    root = logging.getLogger()
    if any(getattr(h, "_galacticus_default", False) for h in root.handlers):
        # Already configured by an earlier call.
        return
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter(_DEFAULT_FORMAT))
    handler._galacticus_default = True  # type: ignore[attr-defined]
    root.addHandler(handler)
    root.setLevel(level)
