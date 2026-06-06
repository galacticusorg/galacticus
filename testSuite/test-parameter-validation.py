#!/usr/bin/env python3
import subprocess
import sys
import os

# Validate a set of test parameter files.
# Andrew Benson (ported to Python)

# Validate the default parameter file.
status = subprocess.run("cd ..; ./scripts/aux/validateParameters.py parameters.xml", shell=True)
if status.returncode == 0:
    print("PASSED: validation of default parameter file")
else:
    print("FAILED: validation of default parameter file")

# Find all validation parameter files and run them.
validationDir = "parameters/validation"
if os.path.isdir(validationDir):
    for fileName in sorted(os.listdir(validationDir)):
        import re
        m = re.search(r'-(valid|invalid)\.xml$', fileName)
        if m:
            validity = m.group(1)
            status = subprocess.run(
                f"cd ..; ./scripts/aux/validateParameters.py testSuite/parameters/validation/{fileName}",
                shell=True
            )
            if (status.returncode == 0 and validity == "valid") or (status.returncode != 0 and validity == "invalid"):
                print(f"PASSED: validation of '{fileName}'")
            else:
                print(f"FAILED: validation of '{fileName}'")
