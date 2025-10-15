#!/usr/bin/env python3
import sys
import os
import subprocess
import validate
import re

# Run models to validate a dark matter only subhalo evolution model against data from the Symphony/COZMIC suite.
# Andrew Benson (09-June-2025)

# Get argument.
if len(sys.argv) != 4:
    sys.exit("FAIL: usage: validate-darkMatterOnlySubhalos-SymphonyCOZMIC.py <suite> <resolution> <simulation>")
suite      = sys.argv[1]
resolution = sys.argv[2]
simulation = sys.argv[3]

# By default the "true" simulation is just the same as the simulation argument.
simulationTrue = simulation

# Create output path.
try:
    os.mkdir("outputs")
except FileExistsError:
    pass

# Construct the parameter files to run.
## Use the Symphony resolutionX1 file as our base.
parameterFiles = "testSuite/parameters/validate_darkMatterOnlySubhalos_Symphony_resolutionX1_CDM.xml parameters/reference/changeSymphony.xml"
## For non-CDM models, add the relevant change files.
matchWDM = re.search(r'^WDM([0-9\.]+)keV',simulation)
if suite == "COZMIC" and matchWDM:
    parameterFiles += " parameters/reference/powerSpectraSuppressed.xml parameters/reference/warmDarkMatter.xml"
    simulationTrue  = "WDM:"+matchWDM.group(1)+"keV"
## For higher resolution models, add the relevant change file.
if   resolution == "X8":
        parameterFiles += " testSuite/parameters/resolutionSymphonyX8.xml"
elif resolution == "X64":
        parameterFiles += " testSuite/parameters/resolutionSymphonyX64.xml"
## Add the specific model file.
if not (suite == "Symphony" and resolution == "X1" and simulation == "CDM"):
    parameterFiles += " testSuite/parameters/validate_darkMatterOnlySubhalos_"+suite+"_resolution"+resolution+"_"+simulationTrue+".xml"

# Run the validation model.
status = subprocess.run("cd ..; mpirun --n 1 ./Galacticus.exe "+parameterFiles,shell=True)
if status.returncode != 0:
    print("FAILED: dark matter-only subhalos validation model ("+suite+"; resolution"+resolution+"; "+simulationTrue+") failed to run")
    sys.exit()

# Extract and validate the likelihoods.
validate.extract("outputs/validate_darkMatterOnlySubhalos_"+suite+"_MilkyWay_resolution"+resolution+"_"+simulationTrue+":MPI0000.hdf5","Dark Matter Only Subhalos ("+suite+" "+simulationTrue+" resolution "+resolution+" Milky Way)","darkMatterOnlySubhalos"+suite+resolution+simulation+"MilkyWay","testSuite/parameters/validate_darkMatterOnlySubhalos_"+suite+"_resolution"+resolution+"_"+simulationTrue+".xml")

print("SUCCESS: dark matter-only subhalos validation model ("+suite+" "+simulationTrue+" resolution "+resolution+" Milky Way)")
