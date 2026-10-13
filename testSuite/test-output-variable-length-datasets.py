#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check that variable-length output datasets are written with the correct HDF5 datatypes and
# self-consistent shapes. The merger tree outputter writes two kinds of variable-length dataset,
# and the datatype of each is chosen separately in the outputter: a rank-1 variable-length buffer
# becomes a vlen of doubles, while a rank-2 one becomes a vlen *of vlens* of doubles. These are
# easy to confuse - `IO_HDF5` also has a writer taking a rank-2 array of the rank-1 type, which is
# a different thing - and getting it wrong creates the dataset with the wrong datatype.
#
# The metallicity-split star formation history model writes one of each: the mass history is
# split by metallicity and so is rank-2, while the times are rank-1.
# Claude (13-October-2026)

subprocess.run("mkdir -p outputs",shell=True)

# Run the model and check for completion.
print("Running model...")
logFileName = "outputs/test-output-variable-length-datasets.log"
with open(logFileName,"w") as logFile:
    status = subprocess.run(
        "cd ..; ./Galacticus.exe testSuite/parameters/outputVariableLengthDatasets.xml",
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

# Examine the datatypes and shapes of the two variable-length datasets.
model      = h5py.File("outputs/outputVariableLengthDatasets.hdf5","r")
nodeData   = model["Outputs/Output1/nodeData"]
nameMass   = "diskStarFormationHistoryMass"
nameTimes  = "diskStarFormationHistoryTimes"
for name in ( nameMass, nameTimes ):
    if name not in nodeData:
        print("FAILED: dataset '"+name+"' is absent from the output")
        sys.exit(0)

# The mass history must be a variable-length dataset whose elements are themselves
# variable-length arrays of doubles. A dataset of plain variable-length doubles here would mean
# the rank-2 buffer had been written with the datatype belonging to the rank-1 case.
massOuter = h5py.check_vlen_dtype(nodeData[nameMass ].dtype)
if massOuter is None:
    print("FAILED: '"+nameMass+"' is not a variable-length dataset")
    sys.exit(0)
massInner = h5py.check_vlen_dtype(massOuter)
if massInner is None:
    print("FAILED: '"+nameMass+"' is a variable-length dataset of "+str(massOuter)+", but should be one of variable-length arrays")
    sys.exit(0)
if massInner != np.float64:
    print("FAILED: '"+nameMass+"' has innermost type "+str(massInner)+", but should have double precision")
    sys.exit(0)
print("SUCCESS: '"+nameMass+"' is a variable-length dataset of variable-length double arrays")

# The times must be a variable-length dataset of doubles.
timesType = h5py.check_vlen_dtype(nodeData[nameTimes].dtype)
if timesType is None:
    print("FAILED: '"+nameTimes+"' is not a variable-length dataset")
    sys.exit(0)
if h5py.check_vlen_dtype(timesType) is not None:
    print("FAILED: '"+nameTimes+"' is a variable-length dataset of variable-length arrays, but should be one of doubles")
    sys.exit(0)
if timesType != np.float64:
    print("FAILED: '"+nameTimes+"' has element type "+str(timesType)+", but should have double precision")
    sys.exit(0)
print("SUCCESS: '"+nameTimes+"' is a variable-length dataset of doubles")

# Both datasets must have one row per output node.
if nodeData[nameMass].shape != nodeData[nameTimes].shape:
    print("FAILED: '"+nameMass+"' has "+str(nodeData[nameMass].shape[0])+" rows but '"+nameTimes+"' has "+str(nodeData[nameTimes].shape[0]))
    sys.exit(0)

# For each node, the mass history must have one entry per metallicity bin - the bin boundaries
# being recorded as an attribute of the dataset - and each of those must span the same times as
# the times history.
metallicityCount = nodeData[nameMass].attrs["metallicity"].size
for iNode in range(nodeData[nameMass].shape[0]):
    massHistory  = nodeData[nameMass ][iNode]
    timesHistory = nodeData[nameTimes][iNode]
    if len(massHistory) != metallicityCount:
        print("FAILED: node "+str(iNode)+" has "+str(len(massHistory))+" metallicity bins, but the dataset records "+str(metallicityCount))
        sys.exit(0)
    for iMetallicity in range(len(massHistory)):
        massBin = np.asarray(massHistory[iMetallicity])
        if massBin.size != timesHistory.size:
            print("FAILED: node "+str(iNode)+", metallicity bin "+str(iMetallicity)+" spans "+str(massBin.size)+" times, but the time history has "+str(timesHistory.size))
            sys.exit(0)
        if not np.all(np.isfinite(massBin)):
            print("FAILED: node "+str(iNode)+", metallicity bin "+str(iMetallicity)+" contains non-finite masses")
            sys.exit(0)
        if np.any(massBin < 0.0):
            print("FAILED: node "+str(iNode)+", metallicity bin "+str(iMetallicity)+" contains negative masses")
            sys.exit(0)
print("SUCCESS: variable-length histories are self-consistent in shape and finite")
sys.exit(0)
