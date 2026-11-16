"""Top-level conftest for tests under python/.

Sets GALACTICUS_EXEC_PATH (used by modules that do
`sys.path.insert(0, ...)` at import time, see Group 2 audit) and prepends
python/ to sys.path so tests can do `from <Module> import ...` regardless
of the cwd they are launched from.

A narrower conftest at python/Galacticus/Build/SourceTree/tests/conftest.py
predates this one and remains in place; both are idempotent.
"""

import os
import sys


_HERE           = os.path.dirname(os.path.abspath(__file__))
GALACTICUS_ROOT = os.path.abspath(os.path.join(_HERE, '..'))

# Set, not setdefault. These tests exercise *this* tree, so the variable must
# name it: modules under test read it at import to find the schemas they
# validate directives against, and deferring to the environment means a shell
# configured for another checkout silently changes what the tests check. On
# master, with it naming a different checkout,
# `test_harvest_records_capability_markers` fails because that tree's
# `objectBuilder.xsd` predates the `withholds` attribute. Only CI, which sets
# the variable to its workspace, was unaffected.
os.environ['GALACTICUS_EXEC_PATH'] = GALACTICUS_ROOT

if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
