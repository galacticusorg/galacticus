#!/usr/bin/env python3
"""Check the export and output-structure merger tree operators on read merger trees.

Two failures fixed in October 2026:

* `mergerTreeOperatorExport` with `includeSubhalos=true` (the default) counted the subhalo history entries it would write
  from the first history that was *gettable*, but wrote them from the first history that was *non-empty*. With `preset`
  satellites the bound mass history is gettable but may be empty while another history is not, so entries were written beyond
  the end of the arrays, and the run aborted with heap corruption ("free(): invalid next size").
* `mergerTreeOperatorOutputStructure` wrote each tree of a forest to a group named after the tree index, which all trees of a
  forest share, so the second tree of a forest wrote into the first one's datasets, and the run aborted with an HDF5
  `H5Dset_extent` error. It now writes the whole forest to one group.

Andrew Benson, Claude
"""
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

import h5py
import numpy as np

# This script runs with the working directory set to testSuite/, while models are run from the repository root.
pathOutputs = os.path.abspath("outputs/mergerTreeExportStructure")
os.makedirs(pathOutputs, exist_ok=True)

def run(label, tree):
    fileParameters = os.path.join(pathOutputs, f"{label}.xml")
    tree.write(fileParameters)
    status = subprocess.run(["./Galacticus.exe", fileParameters], cwd="..", capture_output=True)
    log    = status.stdout.decode()+status.stderr.decode()
    with open(os.path.join(pathOutputs, f"{label}.log"), "w") as fileLog:
        fileLog.write(log)
    if status.returncode != 0:
        print(f"FAILED: model '{label}' did not run (log in testSuite/outputs/mergerTreeExportStructure/{label}.log):")
        print(log[-1500:])
        return False
    return True

# Export with subhalos.
fileExport = os.path.join(pathOutputs, "exported.hdf5")
if os.path.exists(fileExport):
    os.remove(fileExport)
tree = ET.parse("regressions/treeWithNoPrimaryProgenitor.xml")
root = tree.getroot()
for element in root.findall("outputFileName"):
    root.remove(element)
operator = ET.SubElement(root, "mergerTreeOperator", value="export")
ET.SubElement(operator, "outputFileName" , value=fileExport)
ET.SubElement(operator, "includeSubhalos", value="true"    )
ET.SubElement(root, "outputFileName", value=os.path.join(pathOutputs, "export.hdf5"))
if run("export", tree):
    with h5py.File(fileExport, "r") as exported:
        halos      = exported["forestHalos" if "forestHalos" in exported else "haloTrees"]
        indexNode  = halos["nodeIndex"      ][:]
        indexDesc  = halos["descendantIndex" if "descendantIndex" in halos else "descendentIndex"][:]
        indexHost  = halos["hostIndex"      ][:]
        known      = set(indexNode.tolist())
        dangling   = sum(1 for i in indexDesc if i != -1 and i not in known)+sum(1 for i in indexHost if i not in known)
        if len(known) == len(indexNode) and dangling == 0 and np.all(halos["nodeMass"][:] > 0.0):
            print(f"SUCCESS: exported {len(indexNode)} nodes (including subhalos), all with unique indices and resolvable hosts and descendants")
        else:
            print(f"FAILED: exported nodes are inconsistent ({len(indexNode)} nodes, {len(known)} unique indices, {dangling} unresolvable host or descendant indices)")

# Output structure for forests of several trees.
fileTrees = "data/mergerTrees/millenniumTestTreesGLC.hdf5"
indexForest = 102000000
tree = ET.parse("parameters/extractSingleTreeRun.xml")
root = tree.getroot()
root.find("mergerTreeConstructor/fileNames").set("value", "testSuite/"+fileTrees)
for name in ("mergerTreeOperator", "outputFileName"):
    for element in root.findall(name):
        root.remove(element)
operator  = ET.SubElement(root, "mergerTreeOperator", value="outputStructure")
extractor = ET.SubElement(operator, "nodePropertyExtractor", value="multi")
ET.SubElement(extractor, "nodePropertyExtractor", value="nodeIndices")
ET.SubElement(extractor, "nodePropertyExtractor", value="massBasic"  )
fileStructure = os.path.join(pathOutputs, "structure.hdf5")
ET.SubElement(root, "outputFileName", value=fileStructure)
if run("structure", tree):
    with h5py.File(fileTrees, "r") as source, h5py.File(fileStructure, "r") as model:
        trees      = source["treeIndex"]
        halos      = source["haloTrees"]
        expected   = set()
        for first, count, index in zip(trees["firstNode"][:], trees["numberOfNodes"][:], trees["treeIndex"][:]):
            if index == indexForest:
                indexNode = halos["nodeIndex"][first:first+count]
                expected |= set(indexNode[indexNode == halos["hostIndex"][first:first+count]].tolist())
        found = model[f"mergerTreeStructures/tree{indexForest}/nodeIndex"][:]
        if len(found) == len(expected) and set(found.tolist()) == expected:
            print(f"SUCCESS: structure of forest {indexForest} contains each of its {len(expected)} isolated nodes once")
        else:
            print(f"FAILED: structure of forest {indexForest} has {len(found)} nodes ({len(set(found.tolist()))} unique), expected the {len(expected)} isolated nodes of the forest")
sys.exit(0)
