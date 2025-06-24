#!/usr/bin/env python3
import sys
import os
import subprocess
import validate

# Run models to validate a dark matter only subhalo evolution model against data from the COZMIC suite.
# Andrew Benson (09-June-2025)

# Create output path.
try:
    os.mkdir("outputs")
except FileExistsError:
    pass

# Run the validation model COZMIC WDM 3keV.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/validate_darkMatterOnlySubhalos_Symphony_CDM.xml parameters/reference/changeSymphony.xml parameters/reference/powerSpectraSuppressed.xml parameters/reference/warmDarkMatter.xml testSuite/parameters/validate_darkMatterOnlySubhalos_COZMIC_WDM:3keV.xml",shell=True)
if status.returncode != 0:
    print("FAILED: dark matter-only subhalos validation model (WDM 3keV) failed to run")
    sys.exit()

# Extract and validate the likelihoods.
validate.extract("outputs/validate_darkMatterOnlySubhalos_COZMIC_WDM:3keV.hdf5","Dark Matter Only Subhalos (COZMIC WDM 3keV Milky Way)","darkMatterOnlySubhalosCOZMICWDM3keVMilkyWay","testSuite/parameters/validate_darkMatterOnlySubhalos_COZMIC_WDM:3keV.xml")

print("SUCCESS: dark matter-only subhalos validation model (WDM 3keV)")
