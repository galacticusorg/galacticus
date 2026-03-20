#!/usr/bin/env python3
import subprocess
import sys
import h5py

# Check calculations of mass cooled out of the CGM.
# Andrew Benson (ported to Python)

# This model has a cooling CGM hot halo, but no star formation or feedback. So, the mass cooled out of the CGM should equal
# the initial mass minus the final mass.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/cgmMassCooled.xml", shell=True)
if status.returncode != 0:
    print("FAILED: model run:")
    sys.exit(0)
else:
    print("SUCCESS: model run")

# Read the model data and check for consistency.
with h5py.File("outputs/cgmMassCooled.hdf5", "r") as model:
    massCGM    = model["Outputs/Output1/nodeData/hotHaloMass"][:]
    massCooled = model["Outputs/Output1/nodeData/massCooledCGM"][:]

massCooledTarget = 1.0e12 - massCGM[0]
if abs(massCooledTarget - massCooled[0]) < 1.0e6:
    print("SUCCESS: mass cooled out of CGM")
else:
    print(f"FAILED: mass cooled out of CGM: {massCooled[0]} ≇ {massCooledTarget}")
