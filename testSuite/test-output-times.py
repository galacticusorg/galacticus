#!/usr/bin/env python3
import subprocess
import sys
import os
import h5py
import numpy as np

# Run a single Galacticus model and ensure that output times are correct.
# Andrew Benson (ported to Python)

# Simply run the models.
subprocess.run("cd ..; mkdir -p testSuite/outputs; ./Galacticus.exe testSuite/parameters/test-output-times.xml", shell=True)

# Check for outputs.
if not os.path.exists("outputs/test-output-times.hdf5"):
    print("test-output-times.py: FAILED to run Galacticus model")
    sys.exit(0)

# Find all output times.
times = []
with h5py.File("outputs/test-output-times.hdf5", "r") as f:
    outputs = f["Outputs"]
    for outputName in outputs.keys():
        output     = outputs[outputName]
        outputTime = output.attrs["outputTime"]
        times.append(float(outputTime))
times = np.array(sorted(times))

# Construct the expected times.
HubbleConstant = 70.0
megaParsec     = 3.08528229e22
kilo           = 1.0e3
gigaYear       = 3.15576e16
ageUniverse    = (2.0/3.0) * megaParsec / kilo / gigaYear / HubbleConstant
# Fixed times.
timesExpected  = [8.0, 9.0]
# Lookback times.
timesExpected.append(ageUniverse - 2.4)
# Redshifts.
redshifts = np.array([0.0, 1.0, 2.0])
timesExpected.extend(ageUniverse / (1.0 + redshifts)**1.5)
timesExpected = np.array(sorted(timesExpected))

# Compare the times.
if len(times) != len(timesExpected):
    print("test-output-times.py: FAILED - number of times does not match")
elif np.any(np.abs(times - timesExpected) / timesExpected > 1.0e-3):
    print("test-output-times.py: FAILED - times do not match")
else:
    print("test-output-times.py: SUCCESS")
