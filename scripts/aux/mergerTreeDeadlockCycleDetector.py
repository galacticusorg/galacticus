#!/usr/bin/env python3
import sys
import re
import os
import subprocess
import argparse
import shutil

# Detect cycles in Galacticus merger tree deadlock files.
# Uses the `gvpr` tool (https://www.graphviz.org/pdf/gvpr.1.pdf) from the GraphViz package (https://www.graphviz.org/).
# Andrew Benson (02-December-2025)

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='mergerTreeDeadlockCycleDetector.py',description='Detect cycles in Galacticus merger tree deadlock files.')
parser.add_argument('deadlockFileName')
args = parser.parse_args()

# Check that we have the `gvpr` tool.
if shutil.which('gvpr') is None:
    print('This script requires that the `gvpr` tool (part of the GraphViz package) be installed.')
    sys.exit(1)

# Run `gvpr` to detect cycles.
status = subprocess.run("gvpr -f "+os.environ['GALACTICUS_EXEC_PATH']+"/scripts/aux/mergerTreeDeadlockCycleDetector.gvpr "+args.deadlockFileName,shell=True)
if status.returncode != 0:
   print("`gvpr` failed")
   sys.exit(1)

# Parse the cycles file.
cycles = []
cycleFile = open(args.deadlockFileName+".cycles",'r')
for line in cycleFile:
    line = line.strip()
    if line == "Cycle:":
        cycles.append([])
    else:
        match = re.match(r'^(\d+)\s+(\d+)\s+\d+:(\d+)\\ntree:\s+(\d+)\\ntime:\s+([\d\.]+)\\n(.+)',line)
        if not match:
            print("Failed to match line")
            sys.exit(1)
        cycles[-1].append({"selfUniqueID": int(match.group(1)), "dependantUniqueID": int(match.group(2)), "treeID": int(match.group(4)), "selfNodeID": int(match.group(3)), "time": float(match.group(5)), "cause": match.group(6)})
cycleFile.close()

# Initialize a dictionary of diagnoses.
diagnoses = {}

# Report on cycles.
for cycle in cycles:
    print(f"Cycle in tree {cycle[0]['treeID']} consisting of {len(cycle)} nodes:")
    for node in cycle:
        print(f"   ({node['selfUniqueID']}) -> ({node['dependantUniqueID']}) [nodeID: {node['selfNodeID']}; time: {node['time']}; reason: {node['cause']}]")
        # Diagnostics.
        ## Mergee assigned, but no mergerTimestepSatellite included.
        match = re.match(r'mergee\s+\(\s*([\d\.]+)\s*\)',node['cause'])
        if match:
            timeMerge = float(match.group(1))
            if abs(timeMerge-node['time']) < 1.0e-3*node['time']:
                # Find the dependant node.
                for nodeDependant in cycle:
                    if nodeDependant['selfUniqueID'] == node['dependantUniqueID']:
                        if nodeDependant['time'] > 1.001*node['time']:
                            if not 'noMergerTimestepSatellite' in diagnoses:
                                diagnoses['noMergerTimestepSatellite'] = [
                                    'Nodes with mergees present, but these exist after their merging time.',
                                    'Check that you have either',
                                    '   <mergerTreeTimestep value="satellite"/>, or',
                                    '   <mergerTreeTimestep value="standard"/>',
                                    ' set in your parameter file.'
                                    ]

# Report diagnoses.
print()
for diagnosis in diagnoses:
    print("Diagnosis:")
    for line in diagnoses[diagnosis]:
        print(f'   {line}')
    print()
