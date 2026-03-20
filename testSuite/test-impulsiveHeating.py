#!/usr/bin/env python3
import subprocess
import sys

# Test impulsive heating model.
# Andrew Benson (ported to Python)

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Run the model.
status = subprocess.run("export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe testSuite/parameters/impulsiveHeating.xml", shell=True)
if status.returncode == 0:
    print("success: impulsive heating model ran successfully")
else:
    print("FAIL: impulsive heating model failed to run")
