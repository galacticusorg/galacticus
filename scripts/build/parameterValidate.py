#!/usr/bin/env python3
"""Validate Galacticus parameter file(s) against `parameters.catalog.json`.

Usage:
    parameterValidate.py [--catalog <file>] [--exceptions <file>] [--structural] \
                         <paramFile-or-dir> [...]

Reports selector errors (a `value` that names no real implementation), unknown
parameter names (within a functionClass scope), type/enumeration/range
violations, idRef-reference errors, and -- with ``--structural`` -- structural
errors (duplicate top-level parameters, missing values, multiple `<value>`
children).  With a directory argument, recurses over ``*.xml`` and prints a
per-file and aggregate summary -- useful as an audit of `testSuite/parameters/`.

``--catalog <file>`` enables the catalog-aware checks; omit it to run only the
catalog-independent checks (XInclude, references, and structural).  ``--structural``
enables the structural checks (appropriate only for curated, pure parameter
files; they false-positive on non-parameter or tree-data XML).

Exit status is non-zero if any non-excepted error-level finding is reported, so
this doubles as a CI gate.  ``--exceptions <file>`` points at a JSON file listing
findings that are intentional (e.g. files that test the validation/migration
machinery itself); matched findings are reported as "expected" and do not fail.

Andrew Benson (2026).
"""

import json
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, 'python'))

from Galacticus.Parameters.validate import validate_file


def _iter_param_files(arguments):
    for argument in arguments:
        if os.path.isdir(argument):
            for dirpath, dirnames, filenames in os.walk(argument):
                dirnames[:] = sorted(d for d in dirnames if not d.startswith('.'))
                for name in sorted(filenames):
                    if name.endswith('.xml'):
                        yield os.path.join(dirpath, name)
        else:
            yield argument


def _load_exceptions(path):
    """Load the JSON exceptions list, or [] if no path given."""
    if not path:
        return []
    with open(path) as fh:
        return json.load(fh).get('exceptions', [])


def _is_excepted(path, finding, exceptions):
    """True if `finding` matches a documented exception (by file basename,
    finding kind, and the trailing parameter/element name of its path)."""
    base = os.path.basename(path)
    name = finding.path.rsplit('/', 1)[-1]
    return any(
        e.get('file') == base and e.get('kind') == finding.kind
        and e.get('name') == name
        for e in exceptions
    )


def main(argv):
    paths = []
    catalog_path = None
    exceptions_path = None
    structural = False
    i = 1
    while i < len(argv):
        if argv[i] == '--catalog':
            catalog_path = argv[i + 1]
            i += 2
            continue
        if argv[i] == '--exceptions':
            exceptions_path = argv[i + 1]
            i += 2
            continue
        if argv[i] == '--structural':
            structural = True
            i += 1
            continue
        paths.append(argv[i])
        i += 1
    if not paths:
        print("Usage: parameterValidate.py [--catalog <file>] [--exceptions <file>] "
              "[--structural] <paramFile-or-dir> [...]", file=sys.stderr)
        return 2
    catalog = None
    if catalog_path:
        with open(catalog_path) as fh:
            catalog = json.load(fh)
    exceptions = _load_exceptions(exceptions_path)

    files = list(_iter_param_files(paths))
    totals = {'error': 0, 'warning': 0}
    kind_totals = {}
    files_with_errors = 0
    parse_errors = 0
    expected = 0

    for path in files:
        findings, parse_error = validate_file(path, catalog, structural=structural)
        if parse_error:
            parse_errors += 1
            print(f"\n{path}\n  [parse] {parse_error}")
            continue
        reportable = [f for f in findings if not _is_excepted(path, f, exceptions)]
        expected += len(findings) - len(reportable)
        if not reportable:
            continue
        if any(f.level == 'error' for f in reportable):
            files_with_errors += 1
        print(f"\n{path}")
        for finding in reportable:
            totals[finding.level] = totals.get(finding.level, 0) + 1
            kind_totals[finding.kind] = kind_totals.get(finding.kind, 0) + 1
            print(f"  [{finding.level}/{finding.kind}] {finding.path}: {finding.message}")

    print("\n" + "=" * 60)
    print(f"files checked      : {len(files)}")
    print(f"parse errors       : {parse_errors}")
    print(f"files with errors  : {files_with_errors}")
    print(f"errors             : {totals.get('error', 0)}")
    print(f"warnings           : {totals.get('warning', 0)}")
    print(f"expected (excepted): {expected}")
    print("by kind            : " + (", ".join(
        f"{k}={v}" for k, v in sorted(kind_totals.items())) or "none"))
    return 1 if (totals.get('error', 0) or parse_errors) else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
