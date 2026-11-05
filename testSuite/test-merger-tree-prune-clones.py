#!/usr/bin/env python3
import subprocess
import sys
import os
import h5py
import numpy as np
import lxml.etree as ET

# Test the merger tree operator that prunes cloned nodes. Merger trees read from file which contain clones (nodes inserted by
# the "read" constructor at very slightly earlier times than their parent) are processed with and without the "pruneClones"
# operator. The pre-evolution trees are exported in each case, and are then evolved (which would fail if pruning left dangling
# pointers). We check that:
#  * the trees contain clones without pruning (else the test would be vacuous);
#  * with pruning, no clone remains which is the sole progenitor of its parent (clones with siblings must be retained, as they
#    prevent a sibling which is intended to become a subhalo from being promoted to be the primary progenitor; clones which carry
#    events, such as those that are the target of a subhalo promotion, are also retained);
#  * only clones are removed;
#  * clones in all forests and all trees of each forest (not just the first) are pruned.
# Andrew Benson, Claude (23-October-2026)

# Cases: base parameter file, elements replacing those in it, an operator to apply before pruning clones (in both the pruned and
# unpruned models), and whether sole-progenitor clones must exist before pruning. All clones inserted by the "read" constructor in
# these trees have siblings, and so must be retained. Pruning by mass first (preserving primary progenitors) inserts clones as the
# sole progenitor of a parent whose only progenitor was pruned; these must be removed.
millennium = """
<root>
  <mergerTreeConstructor value="read">
    <fileNames value="testSuite/data/mergerTrees/millenniumTestTreesGLC.hdf5"/>
    <presetOrbits value="false"/>
    <presetMergerTimes value="false"/>
    <presetSubhaloIndices value="false"/>
    <presetMergerNodes value="false"/>
    <presetSubhaloMasses value="false"/>
    <presetPositions value="true"/>
    <presetScaleRadii value="false"/>
    <presetAngularMomenta value="false"/>
  </mergerTreeConstructor>
  <cosmologyParameters value="simple">
    <OmegaMatter value="0.25"/>
    <OmegaDarkEnergy value="0.75"/>
    <OmegaBaryon value="0.045"/>
    <HubbleConstant value="73.0"/>
  </cosmologyParameters>
  <cosmologicalMassVariance value="filteredPower">
    <sigma_8 value="0.9"/>
  </cosmologicalMassVariance>
  <componentPosition value="cartesian"/>
</root>
"""
cases = [
    {"label": "noPrimaryProgenitor"         , "parameterFile": "testSuite/regressions/treeWithNoPrimaryProgenitor.xml"  , "replace": None      , "prePrune": None, "soleClonesExpected": False},
    {"label": "mostMassiveProgenitorSubhalo", "parameterFile": "testSuite/parameters/mostMassiveProgenitorIsSubhalo.xml", "replace": None      , "prePrune": None, "soleClonesExpected": False},
    {"label": "millennium"                  , "parameterFile": "testSuite/parameters/test-merger-tree-write-base.xml"   , "replace": millennium, "prePrune": None, "soleClonesExpected": False},
    {"label": "millenniumPrunedByMass"      , "parameterFile": "testSuite/parameters/test-merger-tree-write-base.xml"   , "replace": millennium, "prePrune": '<mergerTreeOperator value="pruneByMass"><massThreshold value="1.0e11"/></mergerTreeOperator>', "soleClonesExpected": True},
]

def structure(fileName):
    # Read the exported trees, returning per-forest arrays.
    forests = {}
    with h5py.File(fileName,'r') as trees:
        halos = {key: trees['forestHalos'][key][:] for key in ('nodeIndex','descendantIndex','hostIndex','redshift')}
        index = trees['forestIndex']
        for forestIndex, first, count in zip(index['forestIndex'][:],index['firstNode'][:],index['numberOfNodes'][:]):
            forests[forestIndex] = {key: value[first:first+count] for key, value in halos.items()}
    return forests

def clones(forest):
    # Count clones, and those which are the sole progenitor of their parent. A clone is a node at the same time as its descendant.
    # (The tolerance on 1+z corresponds to a fractional tolerance in time of about 1e-5, as used by the operator.)
    position = {index: i for i, index in enumerate(forest['nodeIndex'])}
    count    = 0
    sole     = 0
    for i in range(len(forest['nodeIndex'])):
        j = position.get(forest['descendantIndex'][i])
        if j is not None and abs(forest['redshift'][j]-forest['redshift'][i]) < 1.0e-5*(1.0+forest['redshift'][i]):
            count += 1
            if np.count_nonzero(forest['descendantIndex'] == forest['descendantIndex'][i]) == 1:
                sole += 1
    return count, sole

subprocess.run("mkdir -p outputs/",shell=True)
for case in cases:
    results = {}
    for prune in (False, True):
        variant = case["label"]+("Pruned" if prune else "Unpruned")
        # Construct the parameter file.
        parameters = ET.parse("../"+case["parameterFile"])
        root       = parameters.getroot()
        for element in root.findall("outputFileName")+root.findall("mergerTreeOperator"):
            root.remove(element)
        if case["replace"] is not None:
            for element in ET.fromstring(case["replace"]):
                for existing in root.findall(element.tag):
                    root.remove(existing)
                root.append(element)
        ET.SubElement(root,"outputFileName",value="testSuite/outputs/pruneClones"+variant+".hdf5")
        sequence = ET.SubElement(root,"mergerTreeOperator",value="sequence")
        if case["prePrune"] is not None:
            sequence.append(ET.fromstring(case["prePrune"]))
        if prune:
            ET.SubElement(sequence,"mergerTreeOperator",value="pruneClones")
        export   = ET.SubElement(sequence,"mergerTreeOperator",value="export")
        ET.SubElement(export,"outputFileName" ,value="testSuite/outputs/pruneClones"+variant+"Trees.hdf5")
        ET.SubElement(export,"exportFormat"   ,value="galacticus")
        ET.SubElement(export,"includeSubhalos",value="false")
        parameters.write("outputs/pruneClones"+variant+".xml")
        # Remove any previous outputs (the export operator appends to an existing file).
        for suffix in (".hdf5","Trees.hdf5"):
            if os.path.exists("outputs/pruneClones"+variant+suffix):
                os.remove("outputs/pruneClones"+variant+suffix)
        # Run the model.
        print("Running model '"+variant+"'...")
        with open("outputs/pruneClones"+variant+".log","w") as log:
            status = subprocess.run("cd ..; ./Galacticus.exe testSuite/outputs/pruneClones"+variant+".xml",stdout=log,stderr=log,shell=True)
        print("...done ("+str(status.returncode)+")")
        errors = subprocess.run("grep -q -i -e fatal -e aborted -e segmentation \"outputs/pruneClones"+variant+".log\"",shell=True)
        if status.returncode != 0 or errors.returncode == 0:
            print("FAILED: model '"+variant+"' did not run successfully:")
            subprocess.run("cat \"outputs/pruneClones"+variant+".log\"",shell=True)
            sys.exit(0)
        results[prune] = structure("outputs/pruneClones"+variant+"Trees.hdf5")
    # Compare the structures.
    unpruned = results[False]
    pruned   = results[True ]
    if set(unpruned.keys()) != set(pruned.keys()):
        print("FAILED: ["+case["label"]+"] set of forests differs after pruning")
        sys.exit(0)
    countClonesUnpruned     = 0
    countClonesPruned       = 0
    countSoleUnpruned       = 0
    countSolePruned         = 0
    countForestsWithClones  = 0
    for name in sorted(unpruned.keys()):
        clonesUnpruned, soleUnpruned  = clones(unpruned[name])
        clonesPruned  , solePruned    = clones(pruned  [name])
        countClonesUnpruned          += clonesUnpruned
        countClonesPruned            += clonesPruned
        countSoleUnpruned            += soleUnpruned
        countSolePruned              += solePruned
        if clonesUnpruned > 0:
            countForestsWithClones += 1
        # Only clones must be removed.
        removed = len(unpruned[name]['nodeIndex'])-len(pruned[name]['nodeIndex'])
        if removed != clonesUnpruned-clonesPruned:
            print("FAILED: ["+case["label"]+"; forest "+str(name)+"] "+str(removed)+" nodes removed, but "+str(clonesUnpruned-clonesPruned)+" clones pruned")
            sys.exit(0)
    print("["+case["label"]+"] forests: "+str(len(unpruned))+"; forests with clones: "+str(countForestsWithClones)+"; clones (sole progenitors) before/after pruning: "+str(countClonesUnpruned)+" ("+str(countSoleUnpruned)+")/"+str(countClonesPruned)+" ("+str(countSolePruned)+")")
    if countClonesUnpruned == 0:
        print("FAILED: ["+case["label"]+"] no clones found in the unpruned forests - test is not exercising the operator")
        sys.exit(0)
    if case["soleClonesExpected"] and countSoleUnpruned == 0:
        print("FAILED: ["+case["label"]+"] no sole-progenitor clones found in the unpruned forests - test is not exercising pruning")
        sys.exit(0)
    if countSolePruned > 0:
        print("FAILED: ["+case["label"]+"] "+str(countSolePruned)+" sole-progenitor clones remain after pruning")
        sys.exit(0)
    if countClonesUnpruned-countClonesPruned != countSoleUnpruned:
        print("FAILED: ["+case["label"]+"] "+str(countClonesUnpruned-countClonesPruned)+" clones pruned, but "+str(countSoleUnpruned)+" were sole progenitors")
        sys.exit(0)
    print("SUCCESS: ["+case["label"]+"] clones pruned")
