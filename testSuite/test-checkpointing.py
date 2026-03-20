#!/usr/bin/env python3
import subprocess
import sys
import os
import h5py
import numpy as np

# Run a set of Galacticus models to test checkpointing functionality.
# Andrew Benson (ported to Python)

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Run model without checkpointing.
status = subprocess.run("export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe testSuite/parameters/checkpointingNoCheckpoints.xml", shell=True)
if status.returncode != 0:
    print("FAILED: failed to run model without checkpointing")
    sys.exit(0)

# Run model with checkpointing - interrupted.
subprocess.run("cd ..; rm -f checkpoint.chk", shell=True)
subprocess.run("export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe testSuite/parameters/checkpointingCheckpoints.xml", shell=True)
if not os.path.exists("../checkpoint.chk"):
    print("FAILED: failed to produce a checkpoint file")
    sys.exit(0)

# Run model with checkpointing - resuming.
status = subprocess.run("export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe testSuite/parameters/checkpointingResume.xml", shell=True)
if status.returncode != 0:
    print("FAILED: failed to run model resuming from checkpoint")
    sys.exit(0)

# Read model data.
data = {}
for modelName in ("checkpointingNoCheckpoints", "checkpointingCheckpoints"):
    with h5py.File(f"outputs/{modelName}.hdf5", "r") as model:
        nodes = model["Outputs/Output1/nodeData"]
        data[modelName] = {
            "nodeIndex": nodes["nodeIndex"][:],
            "basicMass": nodes["basicMass"][:],
        }

# Compare results.
status_str        = "SUCCESS"
toleranceRelative = 0.01
noChk = data["checkpointingNoCheckpoints"]
chk   = data["checkpointingCheckpoints"]
for i in range(len(noChk["nodeIndex"])):
    match = np.where(chk["nodeIndex"] == noChk["nodeIndex"][i])[0]
    if len(match) != 1:
        continue
    massNoCheckpoints = noChk["basicMass"][i]
    massCheckpoints   = chk["basicMass"][match[0]]
    agrees            = abs(massNoCheckpoints - massCheckpoints) < toleranceRelative * massNoCheckpoints
    if not agrees:
        status_str = "FAIL"
        print(f"\t({i}) {massNoCheckpoints} == {massCheckpoints} ? : FAILURE")

print(f"{status_str}: resume from checkpoint file")
