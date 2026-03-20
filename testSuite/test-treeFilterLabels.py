#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Run a Galacticus model to test node labeling functionality.
# Andrew Benson (ported to Python)

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Run model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/treeFilterLabels.xml", shell=True)
if status.returncode != 0:
    print("FAILED: failed to run model")
    sys.exit(0)

# Read model data.
with h5py.File("outputs/treeFilterLabels.hdf5", "r") as model:
    nodes = model["Outputs/Output1/nodeData"]
    if "nodeLabelLMC" not in nodes:
        print("FAILED: label not output")
        sys.exit(0)
    nodeIsIsolated      = nodes["nodeIsIsolated"][:]
    nodeLabelLMC        = nodes["nodeLabelLMC"][:]
    darkMatterVmax      = nodes["darkMatterProfileDMOVelocityMaximum"][:]
    basicTimeLastIsolated = nodes["basicTimeLastIsolated"][:]

hosts = np.where(nodeIsIsolated == 1)[0]
lmcs  = np.where(nodeLabelLMC  == 1)[0]
print(f"Found {len(lmcs)} LMCs in {len(hosts)} trees")

if len(lmcs) < len(hosts):
    print("FAILED: LMCs not found in all trees")
    sys.exit(0)

if np.any(basicTimeLastIsolated[lmcs] < 11.8):
    print("FAILED: LMCs labelled prior to infall time")
    sys.exit(0)

# Allow some tolerance in Vmax test as peak velocities are determined numerically for tidally-heated subhalos.
if np.any(darkMatterVmax[lmcs] < 54.0):
    print("FAILED: LMCs labelled at low Vmax")
    sys.exit(0)

print("SUCCESS: filtered tree labeling")
