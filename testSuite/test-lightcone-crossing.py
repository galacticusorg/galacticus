#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check that lightcone crossing times are computed correctly.
# Andrew Benson (30-August-2024)

# Run the model and check for completion.
print("Running model...")
status = subprocess.run("mkdir -p outputs",shell=True)
log = open("outputs/test-lightconeCrossing.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/test-lightconeCrossing.xml",stdout=log,stderr=log,shell=True)
log.close()
print("...done ("+str(status)+")")
if status.returncode != 0:
    print("FAILED: model run:")
    subprocess.run("cat outputs/test-lightconeCrossing.log",shell=True)
    sys.exit()
print("Checking for errors...")
status = subprocess.run("grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" outputs/test-lightconeCrossing.log",shell=True)
print("...done ("+str(status)+")")
if status.returncode == 0:
    print("FAILED: model run (errors):")
    subprocess.run("cat outputs/test-lightconeCrossing.log",shell=True)
    sys.exit()
print("SUCCESS: model run")

# Open the model and extract the required properties.
model                         = h5py.File('outputs/test-lightconeCrossing.hdf5','r')
output                        = model['Lightcone/Output1']
nodes                         = output['nodeData']
time                          = nodes['time'                         ][:]
lightconePositionX            = nodes['lightconePositionX'           ][:]
positionPositionX             = nodes['positionPositionX'            ][:]
lightconeRedshiftCosmological = nodes['lightconeRedshiftCosmological'][:]

# Verify that the lightcone and non-lightcone positions are consistent.
if np.allclose(lightconePositionX,positionPositionX*(1.0+lightconeRedshiftCosmological),rtol=1.0e-4):
    print("SUCCESS: lightcone and non-lightcone positions agree")
else:
    print("FAIL: lightcone and non-lightcone positions disagree")

# Verify that the lightcone position is consistent with the expected position at the crossing time.
## Position interpolation is a cubic polynomial is comoving coordinates. The coefficients of that polynomial were computed for our
## test galaxy using Mathematica.
coefficients     = np.array([0.0953563 ,-2.11698  ,3.94789   ,998.074]   )
powers           = np.array([time[0]**3,time[0]**2,time[0]**1,time[0]**0])
distanceExpected = np.sum(coefficients*powers)
if np.allclose(lightconePositionX,distanceExpected,rtol=1.0e-4):
    print("SUCCESS: lightcone position is correct for crossing time")
else:
    print("FAIL: lightcone position is incorrect for crossing time")
