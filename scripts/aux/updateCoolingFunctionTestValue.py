#!/usr/bin/env python3
"""
Update the expected cooling function value in tests.cooling_functions.F90
when the Cloudy CIE cooling function test fails due to table changes.

Usage:
    python3 updateCoolingFunctionTestValue.py <test_log_file> <source_file>
"""

import re
import sys


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <test_log_file> <source_file>")
        sys.exit(1)

    test_log_file = sys.argv[1]
    source_file = sys.argv[2]

    with open(test_log_file) as f:
        log = f.read()

    # The Assert framework note for a failing double-precision scalar uses format e16.8:
    #   FAILED: Cloudy CIE cooling function {<actual_e16.8> ≉ <expected_e16.8>}
    # Extract the actual computed value (the portion before the ≉ symbol).
    m = re.search(r'FAILED: Cloudy CIE cooling function \{(.+?)≉', log)
    if not m:
        print("Could not parse actual value from test log; manual update required")
        sys.exit(0)

    raw = m.group(1).strip()
    # Convert Fortran-style scientific notation to a Python float, then to a
    # Fortran double-precision literal (d-notation) with 16 significant figures.
    val = float(raw.replace('D', 'e').replace('d', 'e'))
    fortran_val = f"{val:.16e}".replace('e', 'd')

    with open(source_file) as f:
        src = f.read()

    updated = re.sub(
        r"(call Assert\('Cloudy CIE cooling function',coolantAtomicCIECloudy,)[0-9.Ded+\-]+(,relTol=1\.0d-6\))",
        rf'\g<1>{fortran_val}\g<2>',
        src
    )
    if updated == src:
        print("WARNING: pattern not found in test source; manual update required")
    else:
        with open(source_file, 'w') as f:
            f.write(updated)
        print(f"Updated Cloudy CIE cooling function expected value to: {fortran_val}")


if __name__ == '__main__':
    main()
