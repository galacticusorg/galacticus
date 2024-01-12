#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np


# Check lengths of adaptive star formation histories are constant.
# Andrew Benson (20-December-2023)

# Run the model and check for completion.
print("Running model...")
status = subprocess.run("mkdir -p testSuite/outputs/regressions",shell=True)
log = open("testSuite/outputs/regressions/adaptiveSFHLengths.log","w")
status = subprocess.run("./Galacticus.exe testSuite/regressions/adaptiveSFHLengths.xml",stdout=log,stderr=log,shell=True)
log.close()
print("...done ("+str(status)+")")
if status.returncode != 0:
    print("FAILED: model run:")
    subprocess.run("cat testSuite/outputs/regressions/adaptiveSFHLengths.log",shell=True)
    sys.exit()

# Check lengths of SFHs.
model        = h5py.File('testSuite/outputs/regressions/adaptiveSFHLengths.hdf5','r')
outputs      = model['Outputs']
timeStep     = 0.1
for outputName, output in model['Outputs'].items():
    time         = output.attrs.get('outputTime')
    nodes        = output['nodeData']
    diskSFHs     = nodes['diskStarFormationHistoryMass']
    spheroidSFHs = nodes['spheroidStarFormationHistoryMass']
    times        = diskSFHs.attrs.get('time')
    if outputName != "Ouptut1":
        if times[-1]-time > timeStep:
            print("FAILED: SFH times inconsistent with output time")
            sys.exit()
    lengthSFH    = -1
    for i in range(diskSFHs.size):
        diskSFH     = diskSFHs    [i]
        spheroidSFH = spheroidSFHs[i]
        if diskSFH.size > 0:
            if diskSFH[0].size > 0:
                if diskSFH[0].size != len(times):
                    print("FAILED: diskSFH and time lengths differ")
                    sys.exit()
                if lengthSFH == -1:
                    lengthSFH = diskSFH[0].size
                elif lengthSFH != diskSFH[0].size:
                    print("FAILED: disk SFH lengths differ")
                    sys.exit()
        if spheroidSFH.size > 0:
            if spheroidSFH[0].size > 0:
                if spheroidSFH[0].size != len(times):
                    print("FAILED: spheroidSFH and time lengths differ")
                    sys.exit()
                if lengthSFH == -1:
                    lengthSFH = spheroidSFH[0].size
                elif lengthSFH != spheroidSFH[0].size:
                    print("FAILED: spheroid SFH lengths differ")
                    sys.exit()
        if diskSFH.size > 0 and spheroidSFH.size > 0:
            if diskSFH[0].size != spheroidSFH[0].size:
                print("FAILED: disk and spheroid SFH lengths differ")
                sys.exit()
                    
# Lengths all agreed - success.
print("success: SFH lengths match")
