#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Test the subhalo promotion works correctly when we elect to always promote the most massive progenitor (even if it is a subhalo
# and non-subhalo progenitors are present).
# Andrew Benson (28-January-2024)
print("Running model...")
status = subprocess.run("mkdir -p outputs/",shell=True)
log = open("outputs/mostMassiveProgenitorIsSubhalo.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/mostMassiveProgenitorIsSubhalo.xml",stdout=log,stderr=log,shell=True)
log.close()
print("...done ("+str(status)+")")
if status.returncode != 0:
    print("FAILED: model run:")
    subprocess.run("cat outputs/mostMassiveProgenitorIsSubhalo.log",shell=True)
    sys.exit()
print("Checking for errors...")
status = subprocess.run("grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" outputs/mostMassiveProgenitorIsSubhalo.log",shell=True)
print("...done ("+str(status)+")")
if status.returncode == 0:
    print("FAILED: model run (errors):")
    subprocess.run("cat outputs/mostMassiveProgenitorIsSubhalo.log",shell=True)
    sys.exit()
print("SUCCESS: model run")

# Check for jumps in the mass accretion history - these would occur if the most massive progenitor was not correctly promoted.
model   = h5py.File('outputs/mostMassiveProgenitorIsSubhalo.hdf5','r')
outputs = model['Outputs']
for output in outputs.keys():
    print(output+":")
    # Get the nodeData group.
    nodes  = outputs[output+"/nodeData"]
    idx    = nodes['nodeIndex'               ][:]
    iso    = nodes['nodeIsIsolated'          ][:]
    times  = nodes['haloAccretionHistoryTime'][:]
    masses = nodes['haloAccretionHistoryMass'][:]
    for i in range(len(iso)):
        ratio = masses[i][1:]/masses[i][0:-1]
        if len(ratio) > 0:
            rmax = np.max(ratio)
            if rmax > 5.0:
                print("FAIL: Large jump in MAH of: "+str(rmax)+" "+str(idx[i]))
                sys.exit()
print("SUCCESS: most massive subhalo promotion")
