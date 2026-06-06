#!/usr/bin/env python3
import subprocess
import sys
import os

# Check that duplicated output property names are caught.
# Andrew Benson (ported to Python)

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Run the model.
with open("outputs/duplicatedOutputPropertyName.log", "w") as logFile:
    status = subprocess.run(
        "cd ..; export OMP_NUM_THREADS=1; ./Galacticus.exe testSuite/parameters/duplicatedOutputPropertyName.xml",
        shell=True, stdout=logFile, stderr=subprocess.STDOUT
    )

if status.returncode == 0:
    print("FAILED: duplicated parameter name was not detected")
else:
    # Check that the error message was given.
    result = subprocess.run("grep -q 'duplicate property name' outputs/duplicatedOutputPropertyName.log", shell=True)
    if result.returncode == 0:
        print("SUCCESS: duplicated parameter name was detected")
    else:
        print("FAILED: duplicated parameter name error message was not given")
