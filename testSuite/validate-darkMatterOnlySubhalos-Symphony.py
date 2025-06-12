#!/usr/bin/env python3
import sys
import os
import subprocess
import validate

# Run models to validate a dark matter only subhalo evolution model against data from the Symphony suite.
# Andrew Benson (09-June-2025)

# Create output path.
try:
    os.mkdir("outputs")
except FileExistsError:
    pass

# Run the validation model for Symphony CDM.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/validate_darkMatterOnlySubhalos_Symphony_CDM.xml parameters/reference/changeSymphony.xml",shell=True)
if status.returncode != 0:
    print("FAILED: dark matter-only subhalos validation model (CDM) failed to run")
    sys.exit()

# Extract and validate the likelihoods.
validate.extract("outputs/validate_darkMatterOnlySubhalos_Symphony_CDM.hdf5","Dark Matter Only Subhalos (Symphony Milky Way)","darkMatterOnlySubhalosSymphonyMilkyWay","testSuite/parameters/validate_darkMatterOnlySubhalos_Symphony_CDM.xml")

print("SUCCESS: dark matter-only subhalos validation model (CDM)")
