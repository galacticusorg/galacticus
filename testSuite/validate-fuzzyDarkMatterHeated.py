#!/usr/bin/env python3
import sys
import os
import subprocess
import validate

# Run models to validate a Fuzzy Dark Matter model using solitonNFWHeated.
# Yu Zhao (10-Angust-2026)

# Create output path.
try:
    os.mkdir("outputs")
except FileExistsError:
    pass

# Run the validation model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/validate_darkMatterOnlySubhalos_Symphony_resolutionX1_CDM.xml parameters/reference/powerSpectraSuppressed.xml parameters/reference/fuzzyDarkMatterHeated.xml testSuite/parameters/resolutionM1e9.xml",shell=True)
if status.returncode != 0:
    print("FAILED: Fuzzy Dark Matter validation model using solitonNFWHeated failed to run")
    sys.exit(0)

print("SUCCESS: Fuzzy Dark Matter validation model using solitonNFWHeated")
