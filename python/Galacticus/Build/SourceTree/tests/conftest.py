"""Pytest configuration for SourceTree regression tests.

Each module under test does
    sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'python'))
at import time.  We honour that by setting GALACTICUS_EXEC_PATH to the
Galacticus root before any test imports run, and by prepending the same
`python/` directory to sys.path so the imports resolve when pytest is
invoked from a different cwd.

The conftest does NOT auto-load every Process submodule — individual tests
import only what they exercise so failures are localised and so we don't
pull in hooks (functionClass, eventHooksStatic, …) that need a populated
`$BUILDPATH` to import their static fixtures.
"""

import os
import sys


_HERE            = os.path.dirname(os.path.abspath(__file__))
GALACTICUS_ROOT  = os.path.abspath(os.path.join(_HERE, '..', '..', '..', '..', '..'))

os.environ.setdefault('GALACTICUS_EXEC_PATH', GALACTICUS_ROOT)

_PYTHON_DIR = os.path.join(GALACTICUS_ROOT, 'python')
if _PYTHON_DIR not in sys.path:
    sys.path.insert(0, _PYTHON_DIR)
