#!/usr/bin/env python3
"""Check that the progenitor mass function analysis counts each parent halo once in forests of several trees.

The `progenitorMassFunction` analysis normalizes its function by the summed weight of the parent halos, which it accumulates
by walking the whole forest each time it is told that a new tree is being output. That happens once for each tree of a
forest, so in a forest of N trees each parent used to be counted N times (fixed October 2026).

This script takes the single-tree forests of `testSuite/data/mergerTrees/millenniumTestTreesGLC.hdf5` and writes three tree
files: one with the forests as they are, and two in which each forest also contains a second tree, a single isolated halo
whose mass is below the range of parent masses (so it contributes neither parents nor progenitors). In the first the added
tree follows the original one and exists at z=0. In the second it is the first tree of the forest and exists only at z~2,
so that it is removed from the forest (as a finished remnant) before the progenitor output, which tests that the trees of
the forest are left pointing to the new first tree. The progenitor mass functions measured from the three files must be
identical.

Andrew Benson, Claude
"""
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

import h5py
import numpy as np

# Relative tolerance for agreement of the progenitor mass functions.
TOLERANCE              = 1.0e-9
# Mass (in the file's units, 10¹⁰ M☉/h) of the halo added to each forest, below the minimum parent mass used below.
MASS_ADDED             = 0.1
# Redshift of the halo added as the first tree of each forest.
REDSHIFT_ADDED_FIRST   = 2.0
# Redshifts of parents and progenitors. The progenitor redshift is replaced by that of the nearest snapshot.
REDSHIFT_PARENT        = 0.0
REDSHIFT_PROGENITOR    = 1.0

# This script runs with the working directory set to testSuite/, while the models are run from the repository root.
os.makedirs("outputs", exist_ok=True)

# Build the two tree files.
with h5py.File("data/mergerTrees/millenniumTestTreesGLC.hdf5", "r") as source:
    trees        = source["treeIndex"]
    halos        = source["haloTrees"]
    descendant   = halos["descendentIndex"][:]
    redshifts    = np.unique(halos["redshift"][:])
    redshiftProgenitor = redshifts[np.argmin(np.abs(redshifts-REDSHIFT_PROGENITOR))]
    forests      = [(start, count) for start, count in zip(trees["firstNode"][:], trees["numberOfNodes"][:]) if np.sum(descendant[start:start+count] < 0) == 1]
    if len(forests) < 10:
        print(f"FAILED: expected at least 10 single-tree forests in the test file, found {len(forests)}")
        sys.exit(0)
    indexNodeNext = halos["nodeIndex"][:].max()+1
    variants = (("Single", False, False), ("Multiple", True, False), ("MultipleFirst", True, True))
    for label, addTree, addFirst in variants:
        fileName = f"outputs/progenitorMassFunctionForests{label}.hdf5"
        with h5py.File(fileName, "w") as target:
            for name in source.keys():
                if name not in ("haloTrees", "treeIndex"):
                    source.copy(source[name], target, name=name)
            haloTrees = target.create_group("haloTrees")
            for key, value in source["haloTrees"].attrs.items():
                haloTrees.attrs[key] = value
            datasets   = {name: [] for name in halos.keys() if halos[name].shape != ()}
            firstNode  = []
            countNodes = []
            indexTree  = []
            countTotal = 0
            for i, (start, count) in enumerate(forests):
                firstNode.append(countTotal)
                indexTree.append(trees["treeIndex"][start == trees["firstNode"][:]][0])
                countAdded = 0
                if addTree:
                    # Make a single, isolated, low-mass halo, offset from the forest's root.
                    root = start+np.nonzero(descendant[start:start+count] < 0)[0][0]
                    added = {name: halos[name][root:root+1].copy() for name in datasets}
                    added["nodeIndex"         ][0] = indexNodeNext+i
                    added["descendentIndex"   ][0] = -1
                    added["hostIndex"         ][0] = indexNodeNext+i
                    added["nodeMass"          ][0] = MASS_ADDED
                    added["redshift"          ][0] = redshifts[np.argmin(np.abs(redshifts-REDSHIFT_ADDED_FIRST))] if addFirst else 0.0
                    added["position"          ][0] = halos["position"][root]+1.0
                    added["halfMassRadius"    ][0] = 0.01
                    added["particleIndexCount"][0] = 0
                    countAdded = 1
                if addTree and addFirst:
                    for name in datasets:
                        datasets[name].append(added[name])
                for name in datasets:
                    datasets[name].append(halos[name][start:start+count])
                if addTree and not addFirst:
                    for name in datasets:
                        datasets[name].append(added[name])
                countNodes.append(count+countAdded)
                countTotal += count+countAdded
            for name, values in datasets.items():
                haloTrees.create_dataset(name, data=np.concatenate(values))
                for key, value in halos[name].attrs.items():
                    haloTrees[name].attrs[key] = value
            # In the source file "descendantIndex" is a reference to "descendentIndex" (an older spelling). Write it directly.
            haloTrees.create_dataset("descendantIndex", data=haloTrees["descendentIndex"][:])
            treeIndex = target.create_group("treeIndex")
            treeIndex.create_dataset("firstNode"    , data=np.array(firstNode , dtype=trees["firstNode"    ].dtype))
            treeIndex.create_dataset("numberOfNodes", data=np.array(countNodes, dtype=trees["numberOfNodes"].dtype))
            treeIndex.create_dataset("treeIndex"    , data=np.array(indexTree , dtype=trees["treeIndex"    ].dtype))
            for key, value in source.attrs.items():
                target.attrs[key] = value

# Write a parameter file for the model used by test-Millennium-tree-builder.py (which reads these trees), adding the
# progenitor mass function analysis.
def parameterFile(treeFile, outputFile):
    tree = ET.parse("parameters/extractSingleTreeRun.xml")
    root = tree.getroot()
    def add(parent, name, value=None):
        element = ET.SubElement(parent, name)
        if value is not None:
            element.set("value", value)
        return element
    for name in ("mergerTreeOutputter", "outputTimes", "outputFileName", "outputAnalysis"):
        for element in root.findall(name):
            root.remove(element)
    constructor = root.find("mergerTreeConstructor")
    constructor.find("fileNames").set("value", treeFile)
    add(constructor, "outputTimeSnapTolerance", "1.0e-3")
    add(root, "mergerTreeOutputter", "analyzer")
    add(add(root, "outputTimes", "list"), "redshifts", f"{redshiftProgenitor:.6f} {REDSHIFT_PARENT:.6f}")
    analysis = add(root, "outputAnalysis", "progenitorMassFunction")
    for name, value in (("label", "Forests"), ("comment", "Progenitor mass function"), ("massRatioMinimum", "1.0e-3"),
                        ("massRatioMaximum", "1.0"), ("countMassRatio", "10"), ("massParentMinimum", "5.0e11"),
                        ("massParentMaximum", "1.0e14"), ("redshiftProgenitor", f"{redshiftProgenitor:.6f}"),
                        ("redshiftParent", f"{REDSHIFT_PARENT:.6f}"), ("alwaysIsolatedOnly", "false"),
                        ("virialDensityContrastDefinition", "sphericalCollapseClsnlssMttrCsmlgclCnstnt"),
                        # Synthetic target data. These play no part in the comparison, but the analysis requires them.
                        ("targetLabel", "synthetic"), ("functionValueTarget", " ".join(["1.0"]*10)),
                        ("functionCovarianceTarget", " ".join("0.01" if i == j else "0.0" for i in range(10) for j in range(10))),
                        ("functionCountTarget", " ".join(["10"]*10))):
        add(analysis, name, value)
    add(root, "outputFileName", outputFile)
    path = outputFile.replace(".hdf5", ".xml")
    tree.write(os.path.join("..", path))
    return path

functions = {}
for label, addTree, addFirst in variants:
    treeFile   = f"testSuite/outputs/progenitorMassFunctionForests{label}.hdf5"
    outputFile = f"testSuite/outputs/progenitorMassFunctionForests{label}Model.hdf5"
    path       = parameterFile(treeFile, outputFile)
    if os.path.exists(os.path.join("..", outputFile)):
        os.remove(os.path.join("..", outputFile))
    status = subprocess.run(f"cd ..; ./Galacticus.exe {path}", shell=True, capture_output=True)
    log = status.stdout.decode() + status.stderr.decode()
    if status.returncode != 0 or not os.path.exists(os.path.join("..", outputFile)):
        print(f"FAILED: the progenitor mass function model ({label}) did not run")
        print(log[-2000:])
        sys.exit(0)
    with h5py.File(os.path.join("..", outputFile), "r") as model:
        try:
            functions[label] = model["analyses/progenitorMassFunctionForests/progenitorMassFunction"][:]
        except KeyError as e:
            print(f"FAILED: progenitor mass function missing from the {label} output: {e}")
            sys.exit(0)
print("SUCCESS: models run")

if not np.any(functions["Single"] > 0.0):
    print("FAILED: the progenitor mass function is empty")
    sys.exit(0)
descriptions = {"Multiple": "a second tree (at z=0) is added to", "MultipleFirst": f"a first tree (at z={REDSHIFT_ADDED_FIRST:.0f}, later removed as a remnant) is added to"}
for label, description in descriptions.items():
    if np.allclose(functions[label], functions["Single"], rtol=TOLERANCE, atol=0.0):
        print(f"SUCCESS: progenitor mass function is unchanged when {description} each of {len(forests)} forests")
    else:
        nonZero = functions["Single"] > 0.0
        ratio   = functions[label][nonZero]/functions["Single"][nonZero]
        print(f"FAILED: progenitor mass function changes when {description} each forest (ratio {ratio.min():.4g} to {ratio.max():.4g})")
sys.exit(0)
