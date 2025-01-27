#!/usr/bin/env python3
import subprocess
import sys
import h5py
import math
import numpy as np

# Test for shell-crossing capture in decaying dark matter models.
# Andrew Benson (18-November-2024)

# Run the model and check for completion
print("Running model...")
status = subprocess.run("mkdir -p outputs",shell=True)
log = open("outputs/test-decayingDM-shell-crossing-capture.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/decayingDM-shell-crossing-capture.xml",stdout=log,stderr=log,shell=True)
log.close()
print("...done ("+str(status)+")")
if status.returncode != 0:
    print("FAILED: model run:")
    subprocess.run("cat outputs/test-decayingDM-shell-crossing-capture.log",shell=True)
    sys.exit()
print("Checking for errors...")
status = subprocess.run("grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" outputs/test-decayingDM-shell-crossing-capture.log",shell=True)
print("...done ("+str(status)+")")
if status.returncode == 0:
    print("FAILED: model run (errors):")
    subprocess.run("cat outputs/test-decayingDM-shell-crossing-capture.log",shell=True)
    sys.exit()
print("SUCCESS: model run")
# Open the model and extract the density profile.
model          = h5py.File('outputs/decayingDM-shell-crossing-capture.hdf5','r')
nodes          = model['Outputs/Output1/nodeData']
densityProfile = nodes['densityProfile'][:]
    
# Check for negative density profiles.
if np.any(densityProfile.flatten() < 0.0):
    print("FAIL: caustic detected")
else:
    print("SUCCESS: no caustic detected")
    
