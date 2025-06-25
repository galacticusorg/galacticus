#!/usr/bin/env python3
import sys
import os
import subprocess
import validate

# Run models to validate a dark matter only subhalo evolution model.
# Andrew Benson (05-August-2022)

# Create output path.
try:
    os.mkdir("outputs")
except FileExistsError:
    pass

# Run the validation model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/validate_darkMatterOnlySubHalos.xml",shell=True)
if status.returncode != 0:
    print("FAILED: dark matter-only subhalos validation model failed to run")
    sys.exit()

# Extract and validate the likelihoods.
validate.extract("outputs/validate_darkMatterOnlySubHalos.hdf5","Dark Matter Only Subhalos","darkMatterOnlySubhalos","testSuite/parameters/validate_darkMatterOnlySubHalos.xml")

print("SUCCESS: dark matter-only subhalos validation model")
