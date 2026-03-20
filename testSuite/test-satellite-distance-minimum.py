#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check calculations of minimum satellite distance from host center.
# Andrew Benson (ported to Python)

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Run the model.
with open("outputs/satelliteDistanceMinimum.log", "w") as logFile:
    status = subprocess.run(
        "cd ..; ./Galacticus.exe testSuite/parameters/satelliteDistanceMinimum.xml",
        shell=True, stdout=logFile, stderr=subprocess.STDOUT
    )
if status.returncode != 0:
    print("FAILED:  model run:")
    with open("outputs/satelliteDistanceMinimum.log") as f:
        print(f.read())
else:
    print("SUCCESS: model run")

# Read the model data and check for consistency.
with h5py.File("outputs/satelliteDistanceMinimum.hdf5", "r") as model:
    outputs = model["Outputs"]
    for outputName in outputs.keys():
        output          = outputs[outputName]
        nodeData        = output["nodeData"]
        isIsolated      = nodeData["nodeIsIsolated"][:]
        distanceMinimum = nodeData["satelliteDistanceMinimum"][:]
        positionX       = nodeData["satellitePositionX"][:]
        positionY       = nodeData["satellitePositionY"][:]
        positionZ       = nodeData["satellitePositionZ"][:]
        distance        = np.sqrt(positionX**2 + positionY**2 + positionZ**2)
        satellites      = np.where(isIsolated == 0)[0]
        tolerance       = 1.0e-6
        status_str      = "SUCCESS" if np.all(distanceMinimum[satellites] <= distance[satellites] * (1.0 + tolerance)) else "FAILED"
        offset          = distanceMinimum[satellites] / distance[satellites] - 1.0
        print(f"Maximum fractional offset: {offset.max()}")
        print(f"{status_str}: minimum distance from host center")
