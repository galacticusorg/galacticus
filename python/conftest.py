# Top-level conftest for tests under python/.
#
# Sets GALACTICUS_EXEC_PATH (used by modules that do
# `sys.path.insert(0, ...)` at import time, see Group 2 audit) and prepends
# python/ to sys.path so tests can do `from <Module> import ...` regardless
# of the cwd they are launched from.
#
# A narrower conftest at python/Galacticus/Build/SourceTree/tests/conftest.py
# predates this one and remains in place; both are idempotent.

import os
import sys


_HERE           = os.path.dirname(os.path.abspath(__file__))
GALACTICUS_ROOT = os.path.abspath(os.path.join(_HERE, '..'))

os.environ.setdefault('GALACTICUS_EXEC_PATH', GALACTICUS_ROOT)

if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
