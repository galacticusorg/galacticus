#!/usr/bin/env python3
"""Check that collapse times from a tree-dependent critical overdensity belong to the tree being evolved.

The `environmental` critical overdensity depends on each tree's environment. When no mass or node dependence is
present, `criticalOverdensityTimeOfCollapse` inverts δ_c(t) using a table built for one tree, and must rebuild it
when the tree changes. A calculation reset on a node of a new tree used to adopt that tree's ID as the one the table
belonged to, so every tree evolved by a thread reused the table of the first tree it tabulated (fixed October 2026).

The model builds four trees in different environments on one thread and outputs each node's collapse redshift. It is
then rerun starting at the second tree. The collapse redshifts of trees 2-4 must be the same in both runs, and the
root halos (which have the same mass, so differ only through their environment) must not all share one collapse
redshift.

Andrew Benson, Claude
"""
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

import h5py
import numpy as np

# Relative tolerance for agreement of collapse redshifts between the two runs.
TOLERANCE_AGREEMENT = 1.0e-6
# Minimum relative spread required in the root halo collapse redshifts across trees.
SPREAD_MINIMUM      = 1.0e-3

# This script runs with the working directory set to testSuite/, while the models are run from the repository root.
parameterFile = "parameters/criticalOverdensityTreeDependent.xml"
os.makedirs("outputs", exist_ok=True)

# Write a variant of the parameter file which begins at the second tree.
tree = ET.parse(parameterFile)
root = tree.getroot()
constructor = root.find("mergerTreeConstructor")
ET.SubElement(constructor, "treeBeginAt", value="2")
root.find("outputFileName").set("value", "testSuite/outputs/criticalOverdensityTreeDependentBegin2.hdf5")
tree.write("outputs/criticalOverdensityTreeDependentBegin2.xml")

runs = {
    "all"   : ("testSuite/parameters/criticalOverdensityTreeDependent.xml", "outputs/criticalOverdensityTreeDependent.hdf5"           ),
    "begin2": ("testSuite/outputs/criticalOverdensityTreeDependentBegin2.xml", "outputs/criticalOverdensityTreeDependentBegin2.hdf5"),
}

# Run the models on a single thread, so that one thread evolves every tree.
environment = dict(os.environ, OMP_NUM_THREADS="1")
redshifts = {}
for name, (parameters, output) in runs.items():
    if os.path.exists(output):
        os.remove(output)
    status = subprocess.run(f"cd ..; ./Galacticus.exe {parameters}", shell=True, capture_output=True, env=environment)
    log = status.stdout.decode() + status.stderr.decode()
    if status.returncode != 0 or not os.path.exists(output):
        print(f"FAILED: the tree-dependent critical overdensity model ('{name}') did not run")
        print(log[-2000:])
        sys.exit(0)
    # Collect collapse redshifts, keyed by tree and node index.
    redshifts[name] = {}
    with h5py.File(output, "r") as f:
        try:
            output_ = f["Outputs/Output1"]
            nodes   = output_["nodeData"]
            for index, start, count in zip(output_["mergerTreeIndex"][:], output_["mergerTreeStartIndex"][:], output_["mergerTreeCount"][:]):
                for i in range(start, start+count):
                    redshifts[name][(int(index), int(nodes["nodeIndex"][i]))] = (float(nodes["haloCollapseRedshift"][i]), bool(nodes["nodeIsIsolated"][i]))
        except KeyError as e:
            print(f"FAILED: expected dataset missing from the '{name}' output: {e}")
            sys.exit(0)
print("SUCCESS: models run")

# Trees 2-4 must give the same collapse redshifts whichever tree a thread began with.
failed = False
compared = 0
for key, (redshift, isolated) in redshifts["begin2"].items():
    if key not in redshifts["all"]:
        print(f"FAILED: node {key[1]} of tree {key[0]} is missing from the run of all trees")
        failed = True
        continue
    redshiftAll = redshifts["all"][key][0]
    compared += 1
    if abs(redshift-redshiftAll) > TOLERANCE_AGREEMENT*abs(redshiftAll):
        print(f"FAILED: tree {key[0]}, node {key[1]}: collapse redshift {redshiftAll:.9f} when all trees are run, but {redshift:.9f} when beginning at tree 2")
        failed = True
if compared == 0:
    print("FAILED: no nodes were compared")
    failed = True
if not failed:
    print(f"SUCCESS: collapse redshifts of {compared} nodes in trees 2-4 do not depend on the first tree evolved")

# Root halos have identical masses, so their collapse redshifts differ only through their environments.
rootRedshifts = np.array([redshift for (index, node), (redshift, isolated) in redshifts["all"].items() if node == 1 and isolated])
if rootRedshifts.size < 4:
    print(f"FAILED: expected four root halos, found {rootRedshifts.size}")
elif (rootRedshifts.max()-rootRedshifts.min()) < SPREAD_MINIMUM*rootRedshifts.mean():
    print(f"FAILED: root halo collapse redshifts are all {rootRedshifts.mean():.9f}, independent of their environments")
else:
    print(f"SUCCESS: root halo collapse redshifts depend on environment ({rootRedshifts.min():.6f} to {rootRedshifts.max():.6f})")
sys.exit(0)
