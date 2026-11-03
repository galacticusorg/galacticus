"""Pytest configuration for SourceTree regression tests.

Each module under test does
    sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'python'))
at import time.  We honor that by setting GALACTICUS_EXEC_PATH to the
Galacticus root before any test imports run, and by prepending the same
`python/` directory to sys.path so the imports resolve when pytest is
invoked from a different cwd.

The conftest does NOT auto-load every Process submodule — individual tests
import only what they exercise so failures are localized and so we don't
pull in hooks (functionClass, eventHooksStatic, …) that need a populated
`$BUILDPATH` to import their static fixtures.
"""

import os
import sys


_HERE            = os.path.dirname(os.path.abspath(__file__))
GALACTICUS_ROOT  = os.path.abspath(os.path.join(_HERE, '..', '..', '..', '..', '..'))

# Set, not setdefault. These tests exercise *this* tree, so the variable must
# name it: modules under test read it at import to find the schemas they
# validate directives against, and deferring to the environment means a shell
# configured for another checkout silently changes what the tests check. On
# master, with it naming a different checkout,
# `test_harvest_records_capability_markers` fails because that tree's
# `objectBuilder.xsd` predates the `withholds` attribute. Only CI, which sets
# the variable to its workspace, was unaffected.
os.environ['GALACTICUS_EXEC_PATH'] = GALACTICUS_ROOT

_PYTHON_DIR = os.path.join(GALACTICUS_ROOT, 'python')
if _PYTHON_DIR not in sys.path:
    sys.path.insert(0, _PYTHON_DIR)
