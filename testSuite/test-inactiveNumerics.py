#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Run models that test that inactive properties are evolved when using an ODE solver
# that does not support inactive property evaluation.
# Andrew Benson (ported to Python)

# Make output directory.
subprocess.run("mkdir -p outputs/", shell=True)

# Run the model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/inactiveNumerics.xml", shell=True)
if status.returncode != 0:
    print("FAIL: inactiveNumerics model failed to run")

# Check that luminosities are non-zero.
nonZeroLuminosities = False
with h5py.File("outputs/inactiveNumerics.hdf5", "r") as model:
    nodeData = model["Outputs/Output1/nodeData"]
    for datasetName in nodeData.keys():
        if "Luminosities" not in datasetName:
            continue
        dataset = nodeData[datasetName][:]
        if np.any(dataset > 0.0):
            nonZeroLuminosities = True
            break

if nonZeroLuminosities:
    print("SUCCESS: non-zero luminosities are present")
else:
    print("FAIL: all luminosities are zero")
