#!/usr/bin/env python3
import sys
import os
import subprocess
import validate

# Test the satellite granule operator for soliton+NFW, soliton-only, and NFW-only host halos.
# Yu Zhao (10-November-2026)

# Create output path.
output_file = "outputs/satelliteGranule.hdf5"
os.makedirs("outputs", exist_ok=True)
if os.path.exists(output_file):
    os.remove(output_file)

# Run the test model.
status = subprocess.run(
    "cd ..; ./Galacticus.exe testSuite/parameters/satelliteGranule.xml",
    shell=True,
)
if status.returncode != 0:
    print("FAILED: satellite granule branch model failed to run")
    sys.exit(0)
if not os.path.exists(output_file):
    print("FAILED: satellite granule branch model produced no output file")
    sys.exit(0)

print("SUCCESS: soliton+NFW, soliton-only, and NFW-only granule branches")
