#!/usr/bin/env python3
import subprocess
import sys
import os
import h5py
import numpy as np
import xml.etree.ElementTree as ET

# Test that a state retrieve reproduces merger trees built using generalized Press-Schechter branching rates from the Farahi
# excursion set solver (issue #1786). Two trees are built and their state is stored. A second model then retrieves that state
# and must rebuild the final tree exactly. This is done for both tree processing orders:
#  * The Farahi solver's tables are allocatable arrays with zero lower bounds. If these are restored with lower bounds of 1
#    (as they once were) the branching rates are wrong after a restore, and the rebuilt tree differs in either order.
#  * When trees are processed in descending order, a state retrieve must resume at the tree whose state was stored. It once
#    resumed at the first tree instead, so rebuilding every tree from the final state.
# Andrew Benson, Claude

outputPath = "outputs/stateRestoreFarahi"
subprocess.run(f"rm -rf {outputPath}; mkdir -p {outputPath}", shell=True)

def readTrees(fileName):
    # Read the node masses of each tree in a merger tree export file, keyed by tree index.
    with h5py.File(fileName, "r") as treesFile:
        forestIndex = treesFile["forestIndex"]
        indices     = forestIndex["forestIndex"  ][:]
        starts      = forestIndex["firstNode"    ][:]
        counts      = forestIndex["numberOfNodes"][:]
        masses      = treesFile["forestHalos/nodeMass"][:]
    return {int(index): masses[start:start+count] for index, start, count in zip(indices, starts, counts)}

def runModel(label, processDescending, store):
    # Build the parameter file for this model from the base parameters, and run it. Returns the name of the tree export file,
    # or None on failure.
    parameters = ET.parse("parameters/state/farahiBase.xml")
    root       = parameters.getroot()
    stateFile  = ET.SubElement(root, "stateFileRoot" if store else "stateRetrieveFileRoot")
    stateFile.set("value", f"testSuite/{outputPath}/{label}State")
    mode = "store" if store else "retrieve"
    root.find("outputFileName").set("value", f"testSuite/{outputPath}/{label}_{mode}.hdf5")
    root.find("mergerTreeConstructor/processDescending").set("value", "true" if processDescending else "false")
    root.find("mergerTreeOperator/outputFileName").set("value", f"testSuite/{outputPath}/{label}_{mode}_trees.hdf5")
    parameterFileName = f"{outputPath}/{label}_{mode}.xml"
    parameters.write(parameterFileName)
    # Run serially, so that the stored state is that from the start of the final tree.
    status = subprocess.run(
        f"export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe testSuite/{parameterFileName} > testSuite/{outputPath}/{label}_{mode}.log 2>&1",
        shell=True
    )
    if status.returncode != 0:
        print(f"FAILED: {label} {mode} model failed (see {outputPath}/{label}_{mode}.log)")
        return None
    return f"{outputPath}/{label}_{mode}_trees.hdf5"

failed = False
for label, processDescending in (("descending", True), ("ascending", False)):
    storeFileName = runModel(label, processDescending, store=True)
    if storeFileName is None:
        failed = True
        continue
    retrieveFileName = runModel(label, processDescending, store=False)
    if retrieveFileName is None:
        failed = True
        continue
    treesStore    = readTrees(storeFileName   )
    treesRetrieve = readTrees(retrieveFileName)
    # Trees are numbered 1 and 2. The final tree processed - which the retrieve must rebuild, and only that tree - is tree 1
    # when processing in descending order, and tree 2 otherwise.
    indexFinal = 1 if processDescending else 2
    if sorted(treesStore.keys()) != [1, 2]:
        print(f"FAILED: {label}: expected trees [1, 2] in the stored model, found {sorted(treesStore.keys())}")
        failed = True
    elif sorted(treesRetrieve.keys()) != [indexFinal]:
        print(f"FAILED: {label}: expected only tree {indexFinal} in the retrieved model, found {sorted(treesRetrieve.keys())}")
        failed = True
    elif len(treesRetrieve[indexFinal]) != len(treesStore[indexFinal]):
        print(f"FAILED: {label}: tree {indexFinal} has {len(treesRetrieve[indexFinal])} nodes after state retrieve, but {len(treesStore[indexFinal])} originally")
        failed = True
    elif not np.array_equal(treesRetrieve[indexFinal], treesStore[indexFinal]):
        print(f"FAILED: {label}: node masses in tree {indexFinal} changed after state retrieve")
        failed = True
    else:
        print(f"{label}: tree {indexFinal} ({len(treesStore[indexFinal])} nodes) is reproduced exactly after state retrieve")

if not failed:
    print("SUCCESS: state retrieve reproduces merger trees built using Farahi excursion set branching rates")
sys.exit(0)
