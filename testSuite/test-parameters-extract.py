#!/usr/bin/env python3
import subprocess
import h5py
import re
import filecmp
import lxml.etree as ET

# Check that parameter extraction produces consistent results. A simple parameter file is run with Galacticus. A new parameter
# file is then extracted from the HDF5 output. That new parameter file is then run with Galacticus, and another parameter file is
# extracted from this second HDF5 file. The two extracted parameter files are then compared - they should be identical.
# Andrew Benson (06-September-2024)

# Run the model and check for completion.
print("Running model...")
status = subprocess.run("mkdir -p outputs",shell=True)
log = open("outputs/test-parameters-extract.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/parametersExtract.xml; ./scripts/parameters/parametersExtract.py testSuite/outputs/parametersExtract.hdf5 testSuite/outputs/parametersExtract.xml; ./Galacticus.exe testSuite/outputs/parametersExtract.xml; ./scripts/parameters/parametersExtract.py testSuite/outputs/parametersExtract.hdf5 testSuite/outputs/parametersExtractSecond.xml; ./Galacticus.exe testSuite/parameters/parametersExtract.xml --dry-run;  ./scripts/parameters/parametersExtract.py testSuite/outputs/parametersExtract.hdf5 testSuite/outputs/parametersExtractDryRun.xml",stdout=log,stderr=log,shell=True)
log.close()
if status.returncode != 0:
    print("...done ("+str(status)+")")
    print("FAILED: model run:")
    subprocess.run("cat outputs/test-parameters-extract.log",shell=True)
    sys.exit()
else:
    print("...done")
    print("Checking for errors...")
    status = subprocess.run("grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" outputs/test-parameters-extract.log",shell=True)
    if status.returncode == 0:
        print("...done ("+str(status)+")")
        print("FAILED: model run (errors):")
        subprocess.run("cat outputs/test-parameters-extract.log",shell=True)
        sys.exit()
    else:
        print("...done")
        print("SUCCESS: model run")


# Compare the two extracted parameter files - they should be identical.
status = "SUCCESS" if filecmp.cmp("outputs/parametersExtract.xml","outputs/parametersExtractSecond.xml") else "FAILED"
print(status+": consistent extracted parameters")

# Parse the output parameter file.
parameters = ET.parse("outputs/parametersExtractDryRun.xml")

# Check that default parameters have correct values.
linearGrowth  = parameters.findall("./linearGrowth[@value='collisionlessMatter']")
if len(linearGrowth) == 1:
    print('SUCCESS: default `<linearGrowth value="collisionlessMatter"/>` is present')
else:
    print('FAILED: default `<linearGrowth value="collisionlessMatter"/>` is not present')

# Check that references and targets are correctly included.
cosmologyTarget  = parameters.findall("./cosmologyParameters[@id='refCosmo']")
cosmologyPointer = parameters.findall("./transferFunction/cosmologyParameters[@idRef='refCosmo']")
if len(cosmologyTarget ) == 1:
    print('SUCCESS: target `<cosmologyParameters/>` is present')
else:
    print('FAILED: target `<cosmologyParameters/>` is not present')
if len(cosmologyPointer) == 1:
    print('SUCCESS: reference `<cosmologyParameters/>` is present')
else:
    print('FAILED: reference `<cosmologyParameters/>` is not present')
