#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Run a test case which checks that satelliteBoundMass is equal to basicMass in the preset
# satellite class for isolated halos and when no bound mass history is set.
# Andrew Benson (20-September-2016; ported to Python)

# Run the model and check for successful completion.
status = subprocess.run(
    "./Galacticus.exe testSuite/regressions/satellitePresetBoundMassNonZero.xml", shell=True
)
if status.returncode != 0:
    print("FAILED: satellitePresetBoundMassNonZero model failed to complete")
    sys.exit(1)

# Extract required information.
with h5py.File("testSuite/outputs/regressions/satellitePresetBoundMassNonZero.hdf5", "r") as model:
    node_data            = model["Outputs/Output1/nodeData"]
    satellite_bound_mass = node_data["satelliteBoundMass"][:]
    basic_mass           = node_data["basicMass"         ][:]

# Check for consistency.
if not np.all(satellite_bound_mass == basic_mass):
    print("FAILED: satellitePresetBoundMassNonZero: satelliteBoundMass does not always equal basicMass")
    sys.exit(1)

print("SUCCESS: satellitePresetBoundMassNonZero")
