#!/usr/bin/env python3
import subprocess
import sys
import h5py

# Test that the deterministic spin model matches the constraint.
# Andrew Benson (ported to Python)

# Run model.
status = subprocess.run("cd ..; mkdir -p testSuite/outputs; ./Galacticus.exe testSuite/parameters/constrainDeterministicSpins.xml", shell=True)
if status.returncode != 0:
    print("FAILED: Galacticus model failed to run")
    sys.exit(0)

# Extract the log-likelihood and check it is sufficiently high.
with h5py.File("outputs/constrainDeterministicSpins.hdf5", "r") as model:
    logLikelihood = model["analyses/spinDistributionBett2007"].attrs["logLikelihood"]

if logLikelihood > -305.50:
    print(f"SUCCESS: deterministic spin model matches constraint [log \u2112 = {logLikelihood:10.2f}]")
else:
    print(f"FAILED: deterministic spin model does not match constraint [log \u2112 = {logLikelihood:10.2f}]")
