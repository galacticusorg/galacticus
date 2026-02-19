#!/usr/bin/env python3
# Test that pruning for lightcone output does not result in spurious removal of galaxies that should be in the lightcone.
import numpy as np
import h5py
import subprocess
import sys
import os

# Create output path.
try:
    os.mkdir("outputs")
except FileExistsError:
    pass

# Run the models.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/testPruneLightconePruned.xml"  ,shell=True)
if status.returncode != 0:
    print("FAILED: pruned model failed to run"  )
    sys.exit()
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/testPruneLightconeUnpruned.xml",shell=True)
if status.returncode != 0:
    print("FAILED: unpruned model failed to run")
    sys.exit()

# Read data.
pruned        = h5py.File('outputs/testPruneLightconePruned.hdf5'  ,'r')
unpruned      = h5py.File('outputs/testPruneLightconeUnpruned.hdf5','r')
prunedIDs     = pruned  ['Lightcone/Output1/nodeData/nodeIndex'][:]
unprunedIDs   = unpruned['Lightcone/Output1/nodeData/nodeIndex'][:]
prunedTimes   = pruned  ['Lightcone/Output1/nodeData/time'     ][:]
unprunedTimes = unpruned['Lightcone/Output1/nodeData/time'     ][:]

# We expect three galaxies in the outputs.
if prunedIDs  .shape[0] != 3:
    print("FAILED: pruned model has incorrect number of galaxies in output"  )
    sys.exit()
if unprunedIDs.shape[0] != 3:
    print("FAILED: unpruned model has incorrect number of galaxies in output")
    sys.exit()
    
# Node indices and times must match.
prunedOrder   = np.argsort(prunedIDs  )
unprunedOrder = np.argsort(unprunedIDs)
if np.any(prunedIDs  [prunedOrder] != unprunedIDs  [unprunedOrder]):
    print("FAILED: pruned and unpruned node IDs do not match")
if np.any(prunedTimes[prunedOrder] != unprunedTimes[unprunedOrder]):
    print("FAILED: pruned and unpruned times do not match"   )
    
# If we have not yet failed, we have succeeded.
print("SUCCESS: pruned and unpruned lightcone content matches")
