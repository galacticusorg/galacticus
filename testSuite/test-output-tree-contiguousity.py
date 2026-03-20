#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Run models that test that merger trees are output in contiguous blocks.
# Andrew Benson (ported to Python)

# Make output directory.
subprocess.run("mkdir -p outputs/", shell=True)

# Run the model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/outputTreeContiguosity.xml", shell=True)
if status.returncode != 0:
    print("FAIL: output tree contiguosity model failed to run")
    sys.exit(0)

# Check for contiguous output.
success = True
with h5py.File("outputs/outputTreeContiguosity.hdf5", "r") as model:
    outputs = model["Outputs"]
    for outputName in outputs.keys():
        output          = outputs[outputName]
        nodeData        = output["nodeData"]
        treeIndex       = output["mergerTreeIndex"][:]
        treeStart       = output["mergerTreeStartIndex"][:]
        treeCount       = output["mergerTreeCount"][:]
        treeIndexByNode = nodeData["mergerTreeIndex"][:]
        for i in range(len(treeIndex)):
            indexTarget = treeIndex[i]
            indexStart  = treeStart[i]
            indexEnd    = indexStart + treeCount[i]
            if not np.all(treeIndexByNode[indexStart:indexEnd] == indexTarget):
                success = False
                break

status_str = "SUCCESS" if success else "FAILED"
print(f"{status_str}: merger tree output contiguosity")
