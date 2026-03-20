#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Run models that test the branchless merger tree algorithm.
# Andrew Benson (ported to Python)

# Make output directory.
subprocess.run("mkdir -p outputs/", shell=True)

# Run the branchless model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/mergerTreeBranchless.xml", shell=True)
if status.returncode != 0:
    print("FAIL: merger tree branchless model failed to run")
    sys.exit(0)

success = True
with h5py.File("outputs/mergerTreeBranchless.hdf5", "r") as model:
    structure = model["mergerTreeStructures"]
    for treeName in structure.keys():
        tree = structure[treeName]
        nodeIndex          = tree["nodeIndex"][:]
        parentIndex        = tree["parentIndex"][:]
        nodeIsOnMainBranch = tree["nodeIsOnMainBranch"][:]
        # If this tree is truly branchless then any node which is not on the main branch should
        # have its parent in the subsequent entry in the arrays, and that parent should be on the main branch.
        sideBranch = np.where(nodeIsOnMainBranch == 0)[0]
        mainBranch = sideBranch + 1
        if len(sideBranch) > 0:
            if not (np.all(parentIndex[sideBranch] == nodeIndex[mainBranch]) and
                    np.all(nodeIsOnMainBranch[mainBranch] == 1)):
                success = False
                break

status_str = "SUCCESS" if success else "FAILED"
print(f"{status_str}: branchless merger tree construction")
