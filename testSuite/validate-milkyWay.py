#!/usr/bin/env python3
import sys
import os
import subprocess
import validate

# Run models to validate a Milky Way model.
# Andrew Benson (10-August-2022)

# Create output path.
try:
    os.mkdir("outputs")
except FileExistsError:
    pass

# Run the validation model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/validate_milkyWay.xml",shell=True)
if status.returncode != 0:
    print("FAILED: Milky Way validation model failed to run")
    sys.exit()

# Extract and validate the likelihoods.
validate.extract("outputs/validate_milkyWay.hdf5","Milky Way model","milkyWayModel","testSuite/parameters/validate_milkyWay.xml")

print("SUCCESS: Milky Way validation model")
