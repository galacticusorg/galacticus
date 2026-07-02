#!/usr/bin/env python3
import subprocess
import sys
import os
import h5py
import numpy as np

# Test that equal-mass progenitor ties in merger trees read from file are resolved deterministically by the primary
# progenitor seniority ordering (descending mass; isolated halos outrank subhalos of equal mass; remaining ties broken
# by descending node index), and validate the structural invariants of the constructed trees (via the
# "validateStructure" merger tree operator, which is enabled in both test parameter files).
# Andrew Benson (02-July-2026)

def treeWrite(fileName, nodes):
    # Write a minimal Galacticus-format merger tree file. `nodes` is a list of (index, descendant, host, mass, redshift).
    tree = h5py.File(fileName, "w")
    tree.attrs['formatVersion'] = 2
    cosmology = tree.create_group('cosmology')
    cosmology.attrs['HubbleParam'       ] = 0.7
    cosmology.attrs['OmegaMatter'       ] = 0.3
    cosmology.attrs['OmegaLambda'       ] = 0.7
    cosmology.attrs['OmegaBaryon'       ] = 0.04
    cosmology.attrs['sigma_8'           ] = 0.93
    cosmology.attrs['powerSpectrumIndex'] = 1.0
    cosmology.attrs['transferFunction'  ] = b'BBKS'
    units = tree.create_group('units')
    units.attrs['massUnitsInSI'          ] = 1.98892e30
    units.attrs['massHubbleExponent'     ] = 0
    units.attrs['massScaleFactorExponent'] = 0
    forestHalos = tree.create_group('forestHalos')
    forestHalos.attrs['forestsAreSelfContained'    ] = 1
    forestHalos.attrs['treesHaveSubhalos'          ] = 1
    forestHalos.attrs['haloMassesIncludeSubhalos'  ] = 0
    forestHalos.attrs['positionsArePeriodic'       ] = 0
    forestHalos.attrs['velocitiesIncludeHubbleFlow'] = 0
    data = np.array(nodes, dtype=[('nodeIndex', 'i8'), ('descendantIndex', 'i8'), ('hostIndex', 'i8'), ('nodeMass', 'f8'), ('redshift', 'f8')])
    for name in data.dtype.names:
        forestHalos.create_dataset(name, data=data[name])
    forestIndex = tree.create_group('forestIndex')
    forestIndex.create_dataset('forestIndex'  , data=np.array([1        ], dtype='i8'))
    forestIndex.create_dataset('firstNode'    , data=np.array([0        ], dtype='i8'))
    forestIndex.create_dataset('numberOfNodes', data=np.array([len(data)], dtype='i8'))
    tree.close()

# Create output path.
try:
    os.mkdir("outputs")
except FileExistsError:
    pass

# Define the test cases. Each case gives the tree structure, the node indices of the satellites expected at the final
# output time, and the node indices which must *not* appear there (the tie loser/winner respectively - see the
# corresponding parameter files for a description of each tree).
cases = {
    "mergerTreeReadTiePromotion": {
        # Subhalos 11 and 12 (equal mass) both descend into isolated node 21 which has no isolated progenitor. Node 12
        # (higher index) must be promoted; node 11 must survive as a satellite.
        "nodes": [
            (10, 20, 10, 1.00e12, 1.0),
            (11, 21, 10, 5.00e10, 1.0),
            (12, 21, 10, 5.00e10, 1.0),
            (20, 30, 20, 1.05e12, 0.5),
            (21, 31, 21, 5.50e10, 0.5),
            (30, -1, 30, 1.10e12, 0.0),
            (31, -1, 31, 6.00e10, 0.0),
        ],
        "satellitesExpected" : {11},
        "satellitesForbidden": {12},
    },
    "mergerTreeReadTieIsolatedPriority": {
        # Subhalo 41 and isolated halo 42 (equal mass) both descend into isolated node 51, with
        # [alwaysPromoteMostMassive]=true. Isolated halos outrank subhalos of equal mass, so 42 remains the primary
        # progenitor; 41 must not be promoted and must survive as a satellite.
        "nodes": [
            (40, 50, 40, 1.00e12, 1.0),
            (41, 51, 40, 5.00e10, 1.0),
            (42, 51, 42, 5.00e10, 1.0),
            (50, 60, 50, 1.05e12, 0.5),
            (51, 61, 51, 6.00e10, 0.5),
            (60, -1, 60, 1.10e12, 0.0),
            (61, -1, 61, 6.50e10, 0.0),
        ],
        "satellitesExpected" : {41},
        "satellitesForbidden": {42},
    },
}

status = "SUCCESS"
for name, case in cases.items():
    print(f"Running case: {name}")
    # Build the input tree file.
    treeWrite(f"outputs/{name}_in.hdf5", case['nodes'])
    # Run the model.
    with open(f"outputs/{name}.log", "w") as log:
        run = subprocess.run(f"cd ..; ./Galacticus.exe testSuite/parameters/{name}.xml", stdout=log, stderr=subprocess.STDOUT, shell=True)
    if run.returncode != 0:
        print(f"FAILED: model run: {name}")
        subprocess.run(f"cat outputs/{name}.log", shell=True)
        status = "FAILED"
        continue
    if subprocess.run(f"grep -q -i -e fatal -e aborted outputs/{name}.log", shell=True).returncode == 0:
        print(f"FAILED: model run (errors): {name}")
        subprocess.run(f"cat outputs/{name}.log", shell=True)
        status = "FAILED"
        continue
    # Check the satellite population at the final output time.
    model          = h5py.File(f"outputs/{name}.hdf5", "r")
    nodeData       = model['Outputs/Output1/nodeData']
    nodeIndex      = nodeData['nodeIndex'     ][...]
    nodeIsIsolated = nodeData['nodeIsIsolated'][...]
    satellites     = set(nodeIndex[nodeIsIsolated == 0])
    model.close()
    print(f"   satellites at final time: {sorted(satellites)}")
    if not case['satellitesExpected'] <= satellites:
        print(f"FAILED: expected satellite(s) {sorted(case['satellitesExpected'] - satellites)} not present: {name}")
        status = "FAILED"
    elif satellites & case['satellitesForbidden']:
        print(f"FAILED: forbidden satellite(s) {sorted(satellites & case['satellitesForbidden'])} present - tie resolved incorrectly: {name}")
        status = "FAILED"
    else:
        print(f"SUCCESS: {name}")
print(status+": merger tree read tie resolution")
