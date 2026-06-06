#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Run models that test that the merger tree branch subsampling algorithm.
# Andrew Benson (ported to Python)

# Make output directory.
subprocess.run("mkdir -p outputs/", shell=True)

# Run the subsampled model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/mergerTreeBranchSubsampled.xml", shell=True)
if status.returncode != 0:
    print("FAIL: merger tree branch subsampling model failed to run")
    sys.exit(0)

# Run the not subsampled model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/mergerTreeBranchNotSubsampled.xml", shell=True)
if status.returncode != 0:
    print("FAIL: merger tree branch no subsampling model failed to run")
    sys.exit(0)

# Read data and construct counts of subhalos.
models = [
    {"label": "subsampled",    "fileName": "outputs/mergerTreeBranchSubsampled.hdf5"   },
    {"label": "notSubsampled", "fileName": "outputs/mergerTreeBranchNotSubsampled.hdf5"}
]
for model in models:
    with h5py.File(model["fileName"], "r") as f:
        nodeData = f["Outputs/Output1/nodeData"]
        for key in ("basicMass", "nodeIsIsolated", "nodeSubsamplingWeight"):
            model[key] = nodeData[key][:]

massesLogarithmic = [11, 10, 9, 8]
for model in models:
    model["countSubhalos"]      = np.zeros(len(massesLogarithmic))
    model["countSubhalosError"] = np.zeros(len(massesLogarithmic))
    for i, massLog in enumerate(massesLogarithmic):
        selectionSubhalos = np.where((model["nodeIsIsolated"] == 0) & (np.log10(model["basicMass"]) >= massLog))[0]
        selectionHalos    = np.where( model["nodeIsIsolated"] == 1)[0]
        weightHalos = model["nodeSubsamplingWeight"][selectionHalos].sum()
        model["countSubhalos"][i]      = model["nodeSubsamplingWeight"][selectionSubhalos].sum() / weightHalos
        model["countSubhalosError"][i] = np.sqrt((model["nodeSubsamplingWeight"][selectionSubhalos]**2).sum()) / weightHalos

# Compare counts.
offsetScaled = np.abs(models[0]["countSubhalos"] - models[1]["countSubhalos"]) / np.sqrt(
    models[0]["countSubhalosError"]**2 + models[1]["countSubhalosError"]**2
)
if np.any(offsetScaled > 3.0):
    print("FAIL: merger tree branch subsampling changes subhalo mass function at > 3σ")
else:
    print("SUCCESS: merger tree branch subsampling does not change subhalo mass function at > 3σ")
