#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check the reference output of the merger tree outputter. With `[outputReferences]` set, the
# outputter creates a group per merger tree holding an HDF5 region reference, per property, to
# that tree's block of rows in the corresponding `nodeData` dataset. This verifies that each
# reference resolves to the correct dataset and selects exactly the rows recorded for that tree
# in `mergerTreeStartIndex` and `mergerTreeCount`.
#
# No other parameter file in the repository enables `[outputReferences]`, so without this test
# the reference path is never exercised.
# Claude (13-October-2026)

subprocess.run("mkdir -p outputs",shell=True)

# Run the model and check for completion.
print("Running model...")
logFileName = "outputs/test-output-references.log"
with open(logFileName,"w") as logFile:
    status = subprocess.run(
        "cd ..; ./Galacticus.exe testSuite/parameters/outputReferences.xml",
        stdout=logFile, stderr=logFile, shell=True
    )
if status.returncode != 0:
    print("FAILED: model run")
    subprocess.run("cat "+logFileName,shell=True)
    sys.exit(0)
errors = subprocess.run(
    "grep -q -i -e fatal -e aborted -e \"task failed\" "+logFileName,
    shell=True
)
if errors.returncode == 0:
    print("FAILED: model run (errors in log)")
    subprocess.run("cat "+logFileName,shell=True)
    sys.exit(0)
print("SUCCESS: model run")

model   = h5py.File("outputs/outputReferences.hdf5","r")
outputs = sorted(model["Outputs"], key=lambda name: int(name.replace("Output","")))
if len(outputs) == 0:
    print("FAILED: the output file contains no outputs")
    sys.exit(0)

referencesChecked = 0
treesChecked      = 0
for outputName in outputs:
    output   = model["Outputs"][outputName]
    nodeData = output["nodeData"]
    treeIndex = output["mergerTreeIndex"     ][:]
    treeStart = output["mergerTreeStartIndex"][:]
    treeCount = output["mergerTreeCount"     ][:]
    for index, start, count in zip(treeIndex,treeStart,treeCount):
        groupName = "mergerTree"+str(index)
        if groupName not in output:
            print("FAILED: "+outputName+" has no group '"+groupName+"' for the tree of that index")
            sys.exit(0)
        treeGroup    = output[groupName]
        treesChecked = treesChecked+1
        # The tree group carries its volume weight, which must match the value recorded for it in
        # the output group.
        if "volumeWeight" not in treeGroup.attrs:
            print("FAILED: '"+groupName+"' in "+outputName+" has no 'volumeWeight' attribute")
            sys.exit(0)
        if treeGroup.attrs["volumeWeight"] != output["mergerTreeWeight"][list(treeIndex).index(index)]:
            print("FAILED: the 'volumeWeight' attribute of '"+groupName+"' in "+outputName+" does not match 'mergerTreeWeight'")
            sys.exit(0)
        # Every reference in the group must resolve to the property dataset of the same name and
        # select exactly this tree's rows of it.
        for propertyName in treeGroup:
            if propertyName not in nodeData:
                print("FAILED: '"+groupName+"' in "+outputName+" references '"+propertyName+"', which is not a property of this output")
                sys.exit(0)
            reference = treeGroup[propertyName][()]
            if not isinstance(reference,h5py.RegionReference):
                print("FAILED: '"+groupName+"/"+propertyName+"' in "+outputName+" is not a region reference")
                sys.exit(0)
            referenced = model[reference]
            if referenced.name != nodeData[propertyName].name:
                print("FAILED: '"+groupName+"/"+propertyName+"' in "+outputName+" resolves to '"+referenced.name+"' rather than '"+nodeData[propertyName].name+"'")
                sys.exit(0)
            selected = referenced[reference]
            expected = nodeData[propertyName][start:start+count]
            if selected.shape != expected.shape:
                print("FAILED: '"+groupName+"/"+propertyName+"' in "+outputName+" selects "+str(selected.shape)+" rows, but the tree spans "+str(expected.shape))
                sys.exit(0)
            if not np.array_equal(selected,expected):
                print("FAILED: '"+groupName+"/"+propertyName+"' in "+outputName+" selects rows which differ from rows "+str(start)+":"+str(start+count)+" of the property")
                sys.exit(0)
            referencesChecked = referencesChecked+1

if referencesChecked == 0:
    print("FAILED: no references were found to check")
    sys.exit(0)
print("SUCCESS: "+str(referencesChecked)+" references across "+str(treesChecked)+" tree groups resolve to the correct rows")
sys.exit(0)
