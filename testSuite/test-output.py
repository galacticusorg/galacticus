#!/usr/bin/env python3
import subprocess
import sys
import os
import h5py

# Run a single Galacticus model and ensure that outputs are created.
# Andrew Benson (ported to Python)

# Simply run the models.
subprocess.run("cd ..; mkdir -p testSuite/outputs; ./Galacticus.exe testSuite/parameters/test-output.xml", shell=True)

# Check for outputs.
if not os.path.exists("outputs/test-output.hdf5"):
    print("test-output.py: FAILED to run Galacticus model")
    sys.exit(0)

with h5py.File("outputs/test-output.hdf5", "r") as f:
    if "Outputs" not in f:
        print("test-output.py FAIL - no Outputs group exists")
        sys.exit(0)
    if "Output1" not in f["Outputs"]:
        print("test-output.py FAIL - no Output1 group exists")
        sys.exit(0)
    if "nodeData" not in f["Outputs/Output1"]:
        print("test-output.py FAIL - no nodeData group exists")
        sys.exit(0)
    if "nodeIndex" not in f["Outputs/Output1/nodeData"]:
        print("test-output.py FAIL - no nodeIndex dataset exists")
        sys.exit(0)

print("test-output.py: SUCCESS")
