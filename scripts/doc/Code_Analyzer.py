#!/usr/bin/env python3
"""Scan Fortran90 source code and output a LaTeX file describing it.

Outputs lists of modules used/used-by with hyperlinks, descriptions, and code
line counts.  Collects much more information on subroutine calls, function
calls, module procedures, etc. — could be made to output that data if required.

Andrew Benson (01-May-2010 original Perl); Python port 2026.

Usage: Code_Analyzer.py <sourceDir> <outputFile>
"""

import os
import re
import sys

_exec_path = os.environ.get(
    'GALACTICUS_EXEC_PATH',
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')),
)
sys.path.insert(0, os.path.join(_exec_path, 'python'))
from latex_utils import latex_encode                                   # noqa: E402
from build.fortran_utils import (get_fortran_line, extract_bracketed,  # noqa: E402
                                  extract_variables)

# ---------------------------------------------------------------------------
# Global data structures (populated by process_file, consumed by
# build_modules_hash and output_data).
# ---------------------------------------------------------------------------

# unitID → {unitType, unitName, belongsTo?, codeLines?, contains?,
#            modulesUsed?, usedBy?, subroutinesCalled?, functionCalls?,
#            methods?, comments?, <intrinsic-type lists>…}
units = {}

# moduleName → unitID  (filled by build_modules_hash)
modules = {}

# ---------------------------------------------------------------------------
# Regex constants
# ---------------------------------------------------------------------------

# use statements:  group(3) is the module name.
USE_RE = re.compile(
    r'^\s*use\s*(,\s*intrinsic\s*)*(::)??\s*([a-z0-9_]+)',
    re.IGNORECASE,
)

# Direct subroutine calls:  group(2) is the subroutine name.
CALL_RE = re.compile(
    r'(^|\W)call\s+([a-z0-9_]+)\s*(\(|\s*$)',
    re.IGNORECASE,
)

# Type-bound subroutine calls:  group(2) variable, group(3) method name.
CALL_TYPE_BOUND_RE = re.compile(
    r'(^|\W)call\s+([a-z0-9_]+)\s*%\s*([a-z0-9_]+)\s*(\(|\s*$)',
    re.IGNORECASE,
)

# Module procedure declarations:  group(1) is the procedure name.
MODULE_PROCEDURE_RE = re.compile(
    r'^\s*module\s+procedure\s+([a-z0-9_]+)\s*$',
    re.IGNORECASE,
)

# Derived-type variable declarations:  group(1) type name, group(3) variable list.
DERIVED_TYPE_RE = re.compile(
    r'^\s*type\s*\(\s*([a-z0-9_]+)\s*\)([\sa-z0-9_,:\+\-\*\/\(\)]*::)*\s*([a-z0-9_,:\+\-\*\/\(\)]+)\s*$',
    re.IGNORECASE,
)

# Type-bound procedure declarations:  group(2) method name, group(3) procedure list.
TYPE_BOUND_RE = re.compile(
    r'^\s*(procedure|generic)\s*::\s*([a-z0-9_]+)\s*=>\s*([a-z0-9_,\s]+)$',
    re.IGNORECASE,
)

# Intrinsic variable declarations.
# Each entry: intrinsic label, compiled regex, and 1-based group indices for
# the kind spec (type_grp), attribute list (attrs_grp), and variable list (vars_grp).
#
# character is special: its outer kind/len spec is group 1, but there are two
# nested capturing groups inside it (groups 2 and 3), which pushes attributes
# to group 4 and variables to group 5.
INTRINSIC_DECLARATIONS = {
    'integer': {
        'intrinsic': 'integer',
        'type_grp':  1,
        'attrs_grp': 2,
        'vars_grp':  3,
        'regex': re.compile(
            r'^\s*integer\s*(\(\s*kind\s*=\s*[a-zA-Z0-9_]+\s*\))*([\sa-zA-Z0-9_,:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9_,:=>\+\-\*\/\(\)\[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'real': {
        'intrinsic': 'real',
        'type_grp':  1,
        'attrs_grp': 2,
        'vars_grp':  3,
        'regex': re.compile(
            r'^\s*real\s*(\(\s*kind\s*=\s*[a-zA-Z0-9_]+\s*\))*([\sa-zA-Z0-9_,:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9\._,:=>\+\-\*\/\(\)\[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'double': {
        'intrinsic': 'double precision',
        'type_grp':  1,
        'attrs_grp': 2,
        'vars_grp':  3,
        'regex': re.compile(
            r'^\s*double\s+precision\s*(\(\s*kind\s*=\s*[a-zA-Z0-9_]+\s*\))*([\sa-zA-Z0-9_,:=\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9\._,:=>\+\-\*\/\(\)\[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'logical': {
        'intrinsic': 'logical',
        'type_grp':  1,
        'attrs_grp': 2,
        'vars_grp':  3,
        'regex': re.compile(
            r'^\s*logical\s*(\(\s*kind\s*=\s*[a-zA-Z0-9_]+\s*\))*([\sa-zA-Z0-9_,:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9_,:=>\+\-\*\/\(\)\[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'character': {
        'intrinsic': 'character',
        'type_grp':  1,
        'attrs_grp': 4,  # two nested groups inside the kind/len spec shift attrs to 4
        'vars_grp':  5,
        'regex': re.compile(
            r'^\s*character\s*(\((\s*(len|kind)\s*=\s*[a-zA-Z0-9_,\+\-\*\(\)]+\s*)+\))*([\sa-zA-Z0-9_,:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9_,:=>\+\-\*\/\(\)\[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'procedure': {
        'intrinsic': 'procedure',
        'type_grp':  1,
        'attrs_grp': 2,
        'vars_grp':  3,
        'regex': re.compile(
            r'^\s*procedure\s*(\(\s*[a-zA-Z0-9_]*\s*\))*([\sa-zA-Z0-9_,:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9_,:=>\+\-\*\/\(\)]+)\s*$',
            re.IGNORECASE,
        ),
    },
}

# Unit-opening patterns.  unit_name_grp is the 1-based capturing group that
# holds the unit name.
UNIT_OPENERS = {
    'module': {
        'unit_name_grp': 1,
        'regex': re.compile(
            r'^\s*module\s+(?!procedure\s)([a-z0-9_]+)$',
            re.IGNORECASE,
        ),
    },
    'program': {
        'unit_name_grp': 1,
        'regex': re.compile(
            r'^\s*program\s+([a-z0-9_]+)',
            re.IGNORECASE,
        ),
    },
    'subroutine': {
        'unit_name_grp': 2,
        'regex': re.compile(
            r'^\s*(impure\s+|pure\s+|elemental\s+|recursive\s+|module\s+)*\s*subroutine\s+([a-z0-9_]+)',
            re.IGNORECASE,
        ),
    },
    # Groups: 1=modifiers, 2=return-type keyword, 3=kind-spec, 4=kind/len=, 5=kind|len, 6=name
    'function': {
        'unit_name_grp': 6,
        'regex': re.compile(
            r'^\s*(impure\s+|pure\s+|elemental\s+|recursive\s+)*\s*'
            r'(real\s*|integer\s*|double\s+precision\s*|double\s+complex\s*|character\s*|logical\s*|module\s*)*\s*'
            r'(\(((kind|len)=)??\w*\))*\s*function\s+([a-z0-9_]+)',
            re.IGNORECASE,
        ),
    },
    # Groups: 1=abstract?, 2=interface name (may be empty)
    'interface': {
        'unit_name_grp': 2,
        'regex': re.compile(
            r'^\s*(abstract\s+)??interface\s+([a-z0-9_()/\+\-\*\.=]*)',
            re.IGNORECASE,
        ),
    },
    # Groups: 1=attributes?, 2=::?, 3=type name
    'type': {
        'unit_name_grp': 3,
        'regex': re.compile(
            r'^\s*type\s*(,\s*abstract\s*|,\s*public\s*|,\s*private\s*|,\s*extends\s*\([a-zA-Z0-9_]+\)\s*|,\s*bind\(c\)\s*)*(::)??\s*([a-z0-9_{}\u00a6]+)\s*$',
            re.IGNORECASE,
        ),
    },
}

# Unit-closing patterns.  unit_name_grp is the 1-based capturing group for the
# unit name; interface has no capturing group so unit_name_grp=0 signals that.
UNIT_CLOSERS = {
    'module': {
        'unit_name_grp': 1,
        'regex': re.compile(r'^\s*end\s+module\s+([a-z0-9_]+)', re.IGNORECASE),
    },
    'program': {
        'unit_name_grp': 1,
        'regex': re.compile(r'^\s*end\s+program\s+([a-z0-9_]+)', re.IGNORECASE),
    },
    'subroutine': {
        'unit_name_grp': 1,
        'regex': re.compile(r'^\s*end\s+subroutine\s+([a-z0-9_]+)', re.IGNORECASE),
    },
    'function': {
        'unit_name_grp': 1,
        'regex': re.compile(r'^\s*end\s+function\s+([a-z0-9_]+)', re.IGNORECASE),
    },
    'interface': {
        'unit_name_grp': 0,  # no capturing group; name check skipped via unitType match
        'regex': re.compile(r'^\s*end\s+interface', re.IGNORECASE),
    },
    'type': {
        'unit_name_grp': 1,
        'regex': re.compile(r'^\s*end\s+type\s+([a-z0-9_{}\u00a6]+)', re.IGNORECASE),
    },
}

# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def process_file(file_path):
    """Scan one source file and populate the global `units` dict.

    TODO: implementation added step by step in subsequent tasks.
    """
    pass  # stub


def build_modules_hash():
    """Build `modules` dict and record usedBy back-references in `units`.

    TODO: implementation added in a later task.
    """
    pass  # stub


def output_data(output_file):
    """Write the LaTeX output file.

    TODO: implementation added in a later task.
    """
    pass  # stub


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 3:
        print('Usage: Code_Analyzer.py <sourceDir> <outputFile>', file=sys.stderr)
        sys.exit(1)

    source_dir  = sys.argv[1]
    output_file = sys.argv[2]

    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        print('Error: BUILDPATH environment variable is not set.', file=sys.stderr)
        sys.exit(1)

    # Walk the source tree and process every eligible file.
    for root, _dirs, files in os.walk(source_dir):
        for fname in files:
            process_file(os.path.join(root, fname))

    build_modules_hash()
    output_data(output_file)


if __name__ == '__main__':
    main()
