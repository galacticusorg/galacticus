#!/usr/bin/env python3
"""Emit `parameters.catalog.json` -- a static, machine-readable catalog of the
input parameters accepted by every `functionClass` implementation.

For each implementation the catalog records the parameters it reads (with an
inferred type and the provenance of that inference) and the nested objects it
builds from the parameter tree.  It is the foundation for typed validation of
parameter files and, longer term, a typed Python configuration layer.

Usage:
    parameterCatalog.py <sourceDirectory> [<outputFile>]

`<sourceDirectory>` is the repository root (the directory containing `source/`).
`<outputFile>` defaults to `<sourceDirectory>/parameters.catalog.json`.

Andrew Benson (2026).
"""

import json
import os
import sys

# Make the in-tree Python packages importable when run directly.
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, 'python'))

from Galacticus.Parameters.catalog import build_catalog


def _coverage_summary(catalog):
    """Return per-provenance parameter counts for a quick coverage report."""
    counts = {}
    total = 0
    for implementation in catalog['implementations'].values():
        for parameter in implementation['parameters']:
            counts[parameter['provenance']] = counts.get(parameter['provenance'], 0) + 1
            total += 1
    return total, counts


def main(argv):
    if len(argv) not in (2, 3):
        print("Usage: parameterCatalog.py <sourceDirectory> [<outputFile>]",
              file=sys.stderr)
        return 1
    source_directory = argv[1]
    output_path = argv[2] if len(argv) == 3 else os.path.join(
        source_directory, 'parameters.catalog.json')

    os.environ.setdefault('GALACTICUS_EXEC_PATH', source_directory)
    source_root = os.path.join(source_directory, 'source')
    if not os.path.isdir(source_root):
        print(f"parameterCatalog.py: no source directory at {source_root}",
              file=sys.stderr)
        return 1

    catalog = build_catalog(source_root, log=lambda m: print(f"  {m}", file=sys.stderr))

    with open(output_path, 'w') as fh:
        json.dump(catalog, fh, indent=1, sort_keys=True)
        fh.write('\n')

    total, counts = _coverage_summary(catalog)
    print(f"wrote {output_path}", file=sys.stderr)
    print(f"  functionClasses : {len(catalog['functionClasses'])}", file=sys.stderr)
    print(f"  implementations : {len(catalog['implementations'])}", file=sys.stderr)
    print(f"  parameters      : {total}", file=sys.stderr)
    for provenance in ('explicit', 'default', 'declaration', 'unknown'):
        n = counts.get(provenance, 0)
        pct = (100.0 * n / total) if total else 0.0
        print(f"    {provenance:12s}: {n:5d}  ({pct:5.1f}%)", file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
