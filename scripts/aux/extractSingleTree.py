#!/usr/bin/env python3
import os
import sys
import numpy as np
import h5py
import argparse

# Extract a single merger tree from a merger tree file. Useful for creation of test cases.
# Andrew Benson (19-September-2012; ported to Python 23-August-2024)

# Convert command line arguments to integers.
def restricted_integer(x):
    try:
        x = int(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a integer literal" % (x,))
    return x

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='extractSingleTree.py',description='Extract a single forest from a Galacticus merger tree file to its own file')
parser.add_argument('fromFile'                           )
parser.add_argument('toFile'                             )
parser.add_argument('forestIndex',type=restricted_integer)
args = parser.parse_args()

# Remove the output file.
if os.path.exists(args.toFile):
    os.remove(args.toFile)

# Open the files.
fromFile = h5py.File(args.fromFile,'r')
toFile   = h5py.File(args.  toFile,'w')

# Detect file version and set group names appropriately.
fromFormatVersion = fromFile.attrs['formatVersion'] if "formatVersion" in fromFile.attrs else 1
forestName      = "treeIndex" if fromFormatVersion == 1 else "forestIndex"
halosName       = "haloTrees" if fromFormatVersion == 1 else "forestHalos"

# Find the forest to extract.
print("Locating forest...")
forestIndex = fromFile[forestName][forestName     ][:]
firstNode   = fromFile[forestName]['firstNode'    ][:]
nodeCount   = fromFile[forestName]['numberOfNodes'][:]
selected    = forestIndex == args.forestIndex
if len(np.nonzero(selected)) < 1:
    sys.exit("failed to find request forest")
if len(np.nonzero(selected)) > 1:
    sys.exit("found multiple matching forests")
start       = firstNode[np.nonzero(selected)]
count       = nodeCount[np.nonzero(selected)]
end         = start+count
print("...done")

# Read all halo datasets.
print("Extracting datasets...")
toFile.create_group("forestHalos")
for name, h5obj in fromFile[halosName].items():
    # Ignore groups.
    if not isinstance(h5obj,h5py.Dataset):
        continue
    # Skip the particleIndex datasets as we will need to recompute them later.
    if name == "particleIndexStart" or name == "particleIndexCount":
        continue
    print("   "+name)
    dataset    = fromFile[halosName][name][start[0]:end[0],...]
    toFile["forestHalos"].create_dataset(name,data=dataset)
print("...done")

# Create the forestIndex group.
toFile.create_group("forestIndex")
toFile["forestIndex"].create_dataset("forestIndex"  ,data=np.array([args.forestIndex]))
toFile["forestIndex"].create_dataset("firstNode"    ,data=np.array([               0]))
toFile["forestIndex"].create_dataset("numberOfNodes",data=                     count  )

# Copy the particle data if present.
if "particleIndexStart" in fromFile[halosName]:
    print("Extracting particle data...")
    particleIndexStart = fromFile[halosName]['particleIndexStart'][start[0]:end[0]]
    particleIndexCount = fromFile[halosName]['particleIndexCount'][start[0]:end[0]]
    indices = set()
    for i in range(len(particleIndexStart)):
        if particleIndexStart[i] < 0:
            continue
        indices.update(list(range(particleIndexStart[i],particleIndexStart[i]+particleIndexCount[i]+1)))
    indices = np.sort(np.array(list(indices)))
    toFile.create_group("particles")
    for name, h5obj in fromFile['particles'].items():
        # Ignore groups.
        if not isinstance(h5obj,h5py.Dataset):
            continue
        dataset = fromFile['particles'][name][indices]
        toFile["particles"].create_dataset(name,data=dataset)
    particleIndexStartNew = particleIndexStart
    for i in range(len(particleIndexStartNew)):
        if particleIndexStartNew[i] < 0:
            continue
        particleIndexStartNew[i] = np.nonzero(indices == particleIndexStart[i])[0][0]
    toFile[halosName].create_dataset("particleIndexStart",data=particleIndexStartNew)
    toFile[halosName].create_dataset("particleIndexCount",data=particleIndexCount   )
    print("...done")
    
# Set format version.
toFile.attrs["formatVersion"] = 2

# Copy all attributes.
print("Copying attributes....")
for name, h5obj in fromFile.items():
    # Ignore non-groups.
    if not isinstance(h5obj,h5py.Group):
        continue
    # Translate names between versions.
    nameTo = name
    if name == "treeIndex":
        nameTo = "forestIndex"
    if name == "haloTrees":
        nameTo = "forestHalos"
    # Copy the attributes.
    print("   "+name+" -> "+nameTo)
    if not nameTo in toFile:
        toFile.create_group(nameTo)
    toFile[nameTo].attrs.update(h5obj.attrs)
print("...done")
