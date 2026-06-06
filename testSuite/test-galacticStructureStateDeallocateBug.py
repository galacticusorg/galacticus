#!/usr/bin/env python3
import subprocess
import sys

# Test for galactic structure state deallocate bug.
# Andrew Benson (ported to Python)

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Run the model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/galacticStructureStateDeallocateBug.xml", shell=True)
if status.returncode != 0:
    print("FAILED: model failed to run")
    sys.exit(0)

print("SUCCESS: galactic structure state deallocate bug test")
