#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check construction of constrained merger trees.
# Andrew Benson (ported to Python)

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Run the model.
with open("outputs/constrainedMergerTrees.log", "w") as logFile:
    status = subprocess.run(
        "cd ..; ./Galacticus.exe testSuite/parameters/constrainedMergerTrees.xml",
        shell=True, stdout=logFile, stderr=subprocess.STDOUT
    )
if status.returncode != 0:
    print("FAILED:  model run:")
    with open("outputs/constrainedMergerTrees.log") as f:
        print(f.read())
else:
    print("SUCCESS: model run")

# Define the constraint.
redshiftConstraint = 8.0
massConstraint     = 1.0e12

# Read the model data and check constraint is met.
with h5py.File("outputs/constrainedMergerTrees.hdf5", "r") as model:
    trees = model["mergerTreeStructures"]
    massesClosest = []
    for treeName in trees.keys():
        tree                  = trees[treeName]
        isOnConstrainedBranch = tree["nodeIsConstrained"][:]
        redshift              = tree["redshift"][:]
        mass                  = tree["massBasic"][:]
        constrainedBranch     = np.where((isOnConstrainedBranch == 1) & (redshift >= redshiftConstraint))[0]
        if len(constrainedBranch) > 0:
            redshiftDelta   = np.abs(redshift[constrainedBranch] - redshiftConstraint)
            redshiftClosest = np.argmin(redshiftDelta)
            massesClosest.append(mass[constrainedBranch[redshiftClosest]])

    massesClosest = np.array(massesClosest)
    tolerance     = 1.0e-3
    status_str    = "SUCCESS" if np.all(np.abs(massesClosest - massConstraint) < massConstraint * tolerance) else "FAILED"
    print(f"{status_str}: constrained merger trees")

    # Check propagation of constrained branch indicator.
    output        = model["Outputs/Output1"]
    nodes         = output["nodeData"]
    treeIndices   = output["mergerTreeIndex"][:]
    treeIndex     = nodes["mergerTreeIndex"][:]
    isConstrained = nodes["nodeIsConstrained"][:]
    statusIndicator = "SUCCESS"
    for i in range(len(treeIndices)):
        selection = np.where((treeIndex == treeIndices[i]) & (isConstrained == 1))[0]
        if len(selection) != 1:
            statusIndicator = "FAILED"
            break
    print(f"{statusIndicator}: constrained branch indicator")
