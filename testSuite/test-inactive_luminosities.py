#!/usr/bin/env python3
import subprocess
import sys
import re
import h5py
import numpy as np

# Test inactive luminosity calculations give results identical to active luminosity calculations.
# Andrew Benson (ported to Python)

# Specify types.
types = ("active", "inactive")

# Run the models.
subprocess.run("mkdir -p outputs/inactiveLuminosities", shell=True)
for ltype in types:
    status = subprocess.run(f"cd ..; ./Galacticus.exe testSuite/parameters/{ltype}Luminosities.xml", shell=True)
    if status.returncode != 0:
        print(f"FAILED: model '{ltype}' failed to run")
        sys.exit(0)

# Extract datasets and compare.
models = {}
for ltype in types:
    models[ltype] = h5py.File(f"outputs/inactiveLuminosities/{ltype}Luminosities.hdf5", "r")

success = True
print("Status\tDataset\tFractional error")
try:
    for i in range(1, 5):
        nodes = {ltype: models[ltype][f"Outputs/Output{i}/nodeData"] for ltype in types}
        datasetNames = list(nodes["active"].keys())
        for datasetName in datasetNames:
            if not (re.match(r'^(disk|spheroid)MassStellar$', datasetName) or
                    re.match(r'^(disk|spheroid)LuminositiesStellar', datasetName)):
                continue
            prop = {}
            for ltype in types:
                prop[ltype] = nodes[ltype][datasetName][:]
            errorFractional = np.abs(prop["inactive"] - prop["active"]) / prop["active"]
            status_str = "SUCCESS"
            if errorFractional[0] > 1.5e-3:
                status_str = "FAILED"
                success    = False
            print(f"{status_str}\t{datasetName}\t{errorFractional[0]}")
finally:
    for ltype in types:
        models[ltype].close()

status_str = "SUCCESS" if success else "FAILED: some datasets differed"
print(f"\n{status_str}")
