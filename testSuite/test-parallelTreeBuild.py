#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Run Galacticus models with and without parallelism in tree building.
# Andrew Benson (ported to Python)

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Run models.
status = subprocess.run("export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe testSuite/parameters/parallelTreeBuildSerial.xml", shell=True)
if status.returncode != 0:
    print("FAILED: failed to run serial model")
    sys.exit(0)

status = subprocess.run("export OMP_NUM_THREADS=4; cd ..; ./Galacticus.exe testSuite/parameters/parallelTreeBuildParallel.xml", shell=True)
if status.returncode != 0:
    print("FAILED: failed to run parallel model")
    sys.exit(0)

# Read model data.
data = {}
for modelName in ("parallelTreeBuildSerial", "parallelTreeBuildParallel"):
    with h5py.File(f"outputs/{modelName}.hdf5", "r") as model:
        analysis = model["conditionalMassFunction"]
        data[modelName] = {
            "conditionalMassFunction":      analysis["conditionalMassFunction"][:],
            "conditionalMassFunctionError": analysis["conditionalMassFunctionError"][:],
        }

# Compare results.
serialFlat   = data["parallelTreeBuildSerial"]["conditionalMassFunction"].flat[:]
parallelFlat = data["parallelTreeBuildParallel"]["conditionalMassFunction"].flat[:]
errorFlat    = data["parallelTreeBuildSerial"]["conditionalMassFunctionError"].flat[:]
selection    = np.where(serialFlat > 0.0)[0]
errorNormalized = np.abs(parallelFlat[selection] - serialFlat[selection]) / errorFlat[selection]
status_str   = "FAIL" if np.any(errorNormalized > 4.0) else "SUCCESS"
print(f"{status_str}: parallel tree build")
if status_str == "FAIL":
    print("Conditional mass function:")
    print(f"\t  Serial build: {serialFlat[selection]}")
    print(f"\tParallel build: {parallelFlat[selection]}")
    print(f"\t   Uncertainty: {errorFlat[selection]}")
