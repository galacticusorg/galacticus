#!/usr/bin/env python3
import os
import sys
import glob
import subprocess
import numpy as np
import h5py
import argparse

# For each requested node index, find which merger tree file in a given folder contains a node with that
# index (in the "forestHalos/nodeIndex" dataset), then use extractSingleTree.py to extract the containing
# tree into its own file (e.g. "badTree<index>.hdf5").
# Andrew Benson (08-June-2026)

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='findAndExtractTrees.py',description='Locate and extract the merger trees containing given node indices.')
parser.add_argument('folder'                                                                                                      )
parser.add_argument('indices'    ,type=int,nargs='+'                       ,help='one or more node indices to locate and extract' )
parser.add_argument('--outputDir',action='store',default="."               ,help='directory into which extracted trees are placed')
parser.add_argument('--prefix'   ,action='store',default="badTree"         ,help='prefix for the extracted tree file names'       )
parser.add_argument('--glob'     ,action='store',default="*.hdf5"          ,help='glob pattern used to find merger tree files'    )
args = parser.parse_args()

# Locate the extractSingleTree.py script (it lives alongside this script).
scriptDir         = os.path.dirname(os.path.abspath(__file__))
extractScript     = os.path.join(scriptDir,"extractSingleTree.py")
if not os.path.exists(extractScript):
    sys.exit("unable to find extractSingleTree.py at '"+extractScript+"'")

# Find the merger tree files.
files = sorted(glob.glob(os.path.join(args.folder,args.glob)))
if len(files) == 0:
    sys.exit("no files matching '"+args.glob+"' found in '"+args.folder+"'")

# Make the output directory if necessary.
os.makedirs(args.outputDir,exist_ok=True)

# Track which indices we still need to find.
remaining = set(args.indices)
found     = {}

# Scan each file, reading the nodeIndex dataset and checking for any of the remaining indices.
for fileName in files:
    if not remaining:
        break
    print("Scanning '"+fileName+"'...")
    with h5py.File(fileName,'r') as f:
        formatVersion = f.attrs['formatVersion'] if "formatVersion" in f.attrs else 1
        halosName     = "haloTrees" if formatVersion == 1 else "forestHalos"
        if halosName not in f or "nodeIndex" not in f[halosName]:
            print("   (no '"+halosName+"/nodeIndex' dataset; skipping)")
            continue
        nodeIndex = f[halosName]['nodeIndex'][:]
        nodeSet   = set(np.intersect1d(nodeIndex,np.fromiter(remaining,dtype=nodeIndex.dtype)).tolist())
    for index in nodeSet:
        print("   found node "+str(index)+" in '"+fileName+"'")
        found[index] = fileName
        remaining.discard(index)

# Report any indices that were not found.
for index in sorted(remaining):
    print("WARNING: node index "+str(index)+" was not found in any file")

# Extract each located tree.
for index, fileName in found.items():
    toFile = os.path.join(args.outputDir,args.prefix+str(index)+".hdf5")
    print("Extracting tree containing node "+str(index)+" from '"+fileName+"' into '"+toFile+"'...")
    subprocess.run([sys.executable,extractScript,fileName,toFile,str(index),"--indexType","node"],check=True)

print("...done")
