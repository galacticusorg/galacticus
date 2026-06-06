#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Run a test case for bug #1066052: "Crash due to zero-size spheroids".
# Andrew Benson (12-October-2012; ported to Python)

# Run the model and check for successful completion.
status = subprocess.run("./Galacticus.exe testSuite/parameters/bug1066052.xml", shell=True)
if status.returncode != 0:
    print("FAILED: bug1066052.xml model failed to complete")
    sys.exit(1)

# Check that the spheroid size is non-zero whenever the spheroid mass is non-zero.
with h5py.File("testSuite/outputs/bug1066052.hdf5", "r") as model:
    node_data   = model["Outputs/Output1/nodeData"]
    scale_length = node_data["spheroidRadius"     ][:]
    stellar_mass = node_data["spheroidMassStellar" ][:]
    gas_mass     = node_data["spheroidMassGas"     ][:]

mass      = stellar_mass + gas_mass
bad_nodes = np.where((mass > 0.0) & (scale_length <= 0.0))[0]
if len(bad_nodes) > 0:
    print("FAILED: bug1066052.xml model contains zero-sized spheroids")
    sys.exit(1)

print("SUCCESS: bug1066052")
