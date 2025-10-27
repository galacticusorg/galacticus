#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check that impossible tree constrains can not be met.
# Andrew Benson (20-October-2025)

# Run the model and check for completion.
print("Running model...")
status = subprocess.run("mkdir -p outputs",shell=True)
log = open("testSuite/outputs/constrainedTreeImpossible.log","w")
status = subprocess.run("./Galacticus.exe testSuite/regressions/constrainedTreeImpossible.xml",stdout=log,stderr=log,shell=True)
log.close()
print("...done ("+str(status)+")")
if status.returncode == 0:
    print("FAILED: model run succeeded (but expected failure):")
    subprocess.run("cat testSuite/outputs/constrainedTreeImpossible.log",shell=True)
    sys.exit()
print("Checking for errors...")
status = subprocess.run("grep -q -e \"maximum number of trials exceeded\" testSuite/outputs/constrainedTreeImpossible.log",shell=True)
print("...done ("+str(status)+")")
if status.returncode != 0:
    print("FAILED: expected error not found:")
    subprocess.run("cat testSuite/outputs/constrainedTreeImpossible.log",shell=True)
    sys.exit()
print("SUCCESS: model failed as expected")
