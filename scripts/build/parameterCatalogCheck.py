#!/usr/bin/env python3
"""Sanity-check a generated `parameters.catalog.json`.

Run after `parameterCatalog.py` to guard against a source change that breaks
catalog generation or regresses type-inference coverage.  Exits non-zero (with a
diagnostic) when an invariant fails; intended for CI.

Usage:
    parameterCatalogCheck.py <catalogFile>

Andrew Benson (2026).
"""

import json
import sys

# Invariant thresholds.  These are deliberately loose floors/ceilings -- they
# catch gross regressions (generation broke, coverage collapsed) without being
# brittle to ordinary growth of the source tree.
MIN_FUNCTION_CLASSES = 100
MIN_IMPLEMENTATIONS  = 500
MAX_UNKNOWN_FRACTION = 0.02


def main(argv):
    if len(argv) != 2:
        print("Usage: parameterCatalogCheck.py <catalogFile>", file=sys.stderr)
        return 1
    with open(argv[1]) as fh:
        catalog = json.load(fh)

    function_classes = catalog.get('functionClasses', {})
    implementations  = catalog.get('implementations', {})
    parameters = [p for impl in implementations.values()
                  for p in impl.get('parameters', [])]
    unknown = sum(1 for p in parameters if p.get('provenance') == 'unknown')
    fraction = unknown / len(parameters) if parameters else 1.0

    print(f"functionClasses={len(function_classes)} "
          f"implementations={len(implementations)} "
          f"parameters={len(parameters)} "
          f"unknown={unknown} ({fraction:.1%})")

    failures = []
    if len(function_classes) < MIN_FUNCTION_CLASSES:
        failures.append(f"too few functionClasses: {len(function_classes)} "
                        f"< {MIN_FUNCTION_CLASSES}")
    if len(implementations) < MIN_IMPLEMENTATIONS:
        failures.append(f"too few implementations: {len(implementations)} "
                        f"< {MIN_IMPLEMENTATIONS}")
    if not parameters:
        failures.append("no parameters catalogued")
    elif fraction >= MAX_UNKNOWN_FRACTION:
        failures.append(f"type-inference coverage regressed: {fraction:.1%} "
                        f"unknown >= {MAX_UNKNOWN_FRACTION:.0%}")

    if failures:
        for failure in failures:
            print(f"parameterCatalogCheck: FAIL: {failure}", file=sys.stderr)
        return 1
    print("parameterCatalogCheck: OK")
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
