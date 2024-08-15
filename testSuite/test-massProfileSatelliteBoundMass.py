#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check internal self-consistency of the mass profile output at the satellite bound mass.
# Andrew Benson (09-August-2024)

# Run the model and check for completion
print("Running model...")
status = subprocess.run("mkdir -p outputs",shell=True)
log = open("outputs/test-mass-profile-satellite-bound-mass.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/massProfileSatelliteBoundMass.xml",stdout=log,stderr=log,shell=True)
log.close()
print("...done ("+str(status)+")")
if status.returncode != 0:
    print("FAILED: model run:")
    subprocess.run("cat outputs/test-mass-profile-satellite-bound-mass.log",shell=True)
    sys.exit()
print("Checking for errors...")
status = subprocess.run("grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" outputs/test-mass-profile-satellite-bound-mass.log",shell=True)
print("...done ("+str(status)+")")
if status.returncode == 0:
    print("FAILED: model run (errors):")
    subprocess.run("cat outputs/test-mass-profile-satellite-bound-mass.log",shell=True)
    sys.exit()
print("SUCCESS: model run")

# Open the model and extract the mass profile and satellite bound mass.
model              = h5py.File('outputs/massProfileSatelliteBoundMass.hdf5','r')
nodes              =            model['Outputs/Output1/nodeData']
nodeIndex          =            nodes['nodeIndex'               ][:]
satelliteBoundMass =            nodes['satelliteBoundMass'      ][:]
massProfile        = np.squeeze(nodes['massProfile'             ][:])
massProfileRadius  = np.squeeze(nodes['massProfileRadius'       ][:])

# Check for consistency. Ignore cases where the radius is 1.0e10 (`radiusLarge` in Galacticus) - these correspond to cases where
# the mass shell corresponding to `satelliteBoundMass` is actually unbound and has moved to infinity.
failures = (np.abs(massProfile-satelliteBoundMass) > 1.0e-3*satelliteBoundMass) & (massProfileRadius < 1.0e10)
if np.count_nonzero(failures) > 0:
    print("FAILED: masses are inconsistent")
    print("nodeIndex\tsatelliteBoundMass\tmassProfilet\tmassProfileRadius\t(ratio)")
    for i in range(len(failures)):
        if failures[i]:
            ratio = massProfile[i]/satelliteBoundMass[i]
            print("%d\t%12.6e\t%12.6e\t%12.6e\t(%12.6e)" % ( nodeIndex[i], satelliteBoundMass[i], massProfile[i], massProfileRadius[i], ratio ))
else:
    print("SUCCESS: masses are consistent" )
    
