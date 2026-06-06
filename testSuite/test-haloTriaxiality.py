#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Test the halo triaxiality model of Menker & Benson (2022) by computing the median axis ratios.
# Andrew Benson (ported to Python)

# Run the model.
status = subprocess.run("cd ..; mkdir -p testSuite/outputs; ./Galacticus.exe testSuite/parameters/haloTriaxialityMenkerBenson2022.xml", shell=True)
if status.returncode != 0:
    print("FAIL: Menker & Benson (2022) halo triaxility failed to run")
    sys.exit(0)

# Extract the data.
with h5py.File("outputs/haloTriaxialityMenkerBenson2022.hdf5", "r") as model:
    axisRatio2 = model["Outputs/Output1/nodeData/darkMatterProfileAxisRatiosY"][:]
    axisRatio3 = model["Outputs/Output1/nodeData/darkMatterProfileAxisRatiosZ"][:]

# Get the median axis ratios.
medianAxisRatio2 = np.median(axisRatio2)
medianAxisRatio3 = np.median(axisRatio3)

# Validate the medians.
tolerance              = 0.05
medianAxisRatio2Target = 0.73
medianAxisRatio3Target = 0.53

if abs(medianAxisRatio2 - medianAxisRatio2Target) < tolerance:
    print("success: axis ratio 2 median")
else:
    print(f"FAIL: axis ratio 2 median ({medianAxisRatio2})")

if abs(medianAxisRatio3 - medianAxisRatio3Target) < tolerance:
    print("success: axis ratio 3 median")
else:
    print(f"FAIL: axis ratio 3 median ({medianAxisRatio3})")
