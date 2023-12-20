#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np


# Check lengths of adaptive star formation histories are constant.
# Andrew Benson (20-December-2023)

# Run the model and check for completion.
print("Running model...")
status = subprocess.run("mkdir -p outputs/regressions",shell=True)
log = open("outputs/regressions/adaptiveSFHLengths.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/regressions/adaptiveSFHLengths.xml",stdout=log,stderr=log,shell=True)
log.close()
print("...done ("+str(status)+")")
if status.returncode != 0:
    print("FAILED: model run:")
    subprocess.run("cat outputs/regressions/adaptiveSFHLengths.log",shell=True)
    sys.exit()

# Check lengths of SFHs.
model     = h5py.File('outputs/regressions/adaptiveSFHLengths.hdf5','r')
nodes     = model['Outputs/Output137/nodeData']
diskSFHs  = nodes['diskStarFormationHistoryMass']
times     = diskSFHs.attrs.get('time')
lengthSFH = -1
for sfh in diskSFHs[:]:
    if sfh[0].size > 0:
        if sfh[0].size != len(times):
            print("FAILED: SFH and time lengths differ")
            sys.exit()
        if lengthSFH == -1:
            lengthSFH = sfh[0].size
        elif lengthSFH != sfh[0].size:
            print("FAILED: SFH lengths differ")
            sys.exit()

# Lengths all agreed - success.
print("success: SFH lengths match")
