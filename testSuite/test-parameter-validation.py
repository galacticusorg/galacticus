#!/usr/bin/env python3
import subprocess
import os
import re

# Validate a set of test parameter files using the structural checks of
# scripts/build/parameterValidate.py (the catalog-free `--structural` mode, which
# subsumes the former scripts/aux/validateParameters.py).
# Andrew Benson (ported to Python)

failures = 0


def validate(path):
    return subprocess.run(
        f"cd ..; ./scripts/build/parameterValidate.py --structural {path}",
        shell=True,
    ).returncode


# The default parameter file should be valid.
if validate("parameters.xml") == 0:
    print("PASSED: validation of default parameter file")
else:
    print("FAILED: validation of default parameter file")
    failures += 1

# Validation fixtures encode their expected validity in their file name
# (`-valid.xml` / `-invalid.xml`).
validationDir = "parameters/validation"
if os.path.isdir(validationDir):
    for fileName in sorted(os.listdir(validationDir)):
        match = re.search(r'-(valid|invalid)\.xml$', fileName)
        if not match:
            continue
        validity = match.group(1)
        returncode = validate(f"testSuite/parameters/validation/{fileName}")
        passed = ((returncode == 0 and validity == "valid")
                  or (returncode != 0 and validity == "invalid"))
        if passed:
            print(f"PASSED: validation of '{fileName}'")
        else:
            print(f"FAILED: validation of '{fileName}'")
            failures += 1

if failures > 0:
    raise SystemExit(1)
