#!/usr/bin/env python3
import subprocess
import sys

# Run a power spectrum task using the CosmicEmu nonlinear power spectrum twice. This tests that a
# pre-computed power spectrum file is correctly re-read on subsequent runs.
# Andrew Benson (13-October-2023; ported to Python)

# Run the model and check for successful completion.
for i in range(2):
    status = subprocess.run("./Galacticus.exe testSuite/regressions/cosmicEmu.xml", shell=True)
    if status.returncode != 0:
        print("FAILED: cosmicEmu model failed to complete")
        sys.exit(1)

print("SUCCESS: cosmicEmu")
