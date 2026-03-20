#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Test bar instability model.
# Andrew Benson (ported to Python)

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Run the model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/barInstability.xml", shell=True)
if status.returncode != 0:
    print("FAILED: model failed to run")
    sys.exit(0)

# Read model data.
with h5py.File("outputs/barInstability.hdf5", "r") as model:
    spheroidMassStellar = model["Outputs/Output1/nodeData/spheroidMassStellar"][:]

if spheroidMassStellar[0] > 0.0:
    print("SUCCESS: bar instability")
else:
    print("FAILED: bar instability - no stellar mass in spheroid")
    sys.exit(0)
