#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Run models that test the `nodePropertyExtractorOutputSelector` class.
# Andrew Benson (ported to Python)

# Make output directory.
subprocess.run("mkdir -p outputs/", shell=True)

# Run the model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/outputSelector.xml", shell=True)
if status.returncode != 0:
    print("FAIL: output selector model failed to run")
    sys.exit(0)

# Check for correctly selected outputs.
success = True
with h5py.File("outputs/outputSelector.hdf5", "r") as model:
    outputs = model["Outputs"]
    for outputName in outputs.keys():
        output          = outputs[outputName]
        nodeData        = output["nodeData"]
        expansionFactor = output.attrs["outputExpansionFactor"]
        redshift        = 1.0 / expansionFactor - 1.0
        found           = "mergerTreeIndex" in nodeData
        if abs(redshift - 0.04) < 1.0e-3 or abs(redshift - 2.34) < 1.0e-3:
            # Output property should be present at this redshift.
            if not found:
                success = False
            print(f"{outputName}\tz={redshift}\t{'found (expected)' if found else 'not found (unexpected)'}")
        else:
            # Output property should not be present at this redshift.
            if found:
                success = False
            print(f"{outputName}\tz={redshift}\t{'found (unexpected)' if found else 'not found (expected)'}")

status_str = "SUCCESS" if success else "FAILED"
print(f"{status_str}: nodePropertyExtractorOutputSelector")
