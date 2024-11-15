#!/usr/bin/env python3
import subprocess
import sys
import h5py
import math
import numpy as np

# Test for caustics in decaying dark matter models. These can happen if the shell-crossing radius is approached too closely, and
# can lead to failures in numerical integration of the density profiles resulting in excessively large rotation curves.
# Andrew Benson (28-October-2024)

# Run the model and check for completion
print("Running model...")
status = subprocess.run("mkdir -p outputs",shell=True)
log = open("outputs/test-decayingDM-caustic.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/decayingDM-caustic.xml",stdout=log,stderr=log,shell=True)
log.close()
print("...done ("+str(status)+")")
if status.returncode != 0:
    print("FAILED: model run:")
    subprocess.run("cat outputs/test-decayingDM-caustic.log",shell=True)
    sys.exit()
print("Checking for errors...")
status = subprocess.run("grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" outputs/test-decayingDM-caustic.log",shell=True)
print("...done ("+str(status)+")")
if status.returncode == 0:
    print("FAILED: model run (errors):")
    subprocess.run("cat outputs/test-decayingDM-caustic.log",shell=True)
    sys.exit()
print("SUCCESS: model run")
# Open the model and extract the rotation curve.
model         = h5py.File('outputs/decayingDM-caustic.hdf5','r')
nodes         = model['Outputs/Output1/nodeData']
rotationCurve = nodes['rotationCurve'][:]
    
# Check for excessively large rotation velocities.
if np.any(rotationCurve.flatten() > 30.0):
    print("FAIL: caustics detected")
else:
    print("SUCCESS: no caustics detected")
    
