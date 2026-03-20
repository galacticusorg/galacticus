#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check calculations of maximum host halo mass.
# Andrew Benson (ported to Python)

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Run the model.
with open("outputs/massHostMaximum.log", "w") as logFile:
    status = subprocess.run(
        "cd ..; ./Galacticus.exe testSuite/parameters/massHostMaximum.xml",
        shell=True, stdout=logFile, stderr=subprocess.STDOUT
    )
if status.returncode != 0:
    print("FAILED:  model run:")
    with open("outputs/massHostMaximum.log") as f:
        print(f.read())
else:
    print("SUCCESS: model run")

# Read the model data and check for consistency.
with h5py.File("outputs/massHostMaximum.hdf5", "r") as model:
    outputs = model["Outputs"]
    for outputName in outputs.keys():
        output          = outputs[outputName]
        expansionFactor = output.attrs["outputExpansionFactor"]
        nodeData        = output["nodeData"]
        isIsolated      = nodeData["nodeIsIsolated"][:]
        massHost        = nodeData["satelliteHostHaloMass"][:]
        massHostMaximum = nodeData["massHostMaximum"][:]
        satellites      = np.where(isIsolated == 0)[0]
        redshift        = f"{1.0/expansionFactor - 1.0:3.1f}"
        status_str      = "SUCCESS" if np.all(massHost[satellites] == massHostMaximum[satellites]) else "FAILED"
        print(f"{status_str}: maximum host mass at z={redshift}")
