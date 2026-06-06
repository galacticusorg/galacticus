#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Test that the dark matter halo mass function matches the constraint.
# Andrew Benson (ported to Python)

# Run model.
status = subprocess.run("cd ..; mkdir -p testSuite/outputs; ./Galacticus.exe testSuite/parameters/constrainHaloMassFunction.xml", shell=True)
if status.returncode != 0:
    print("FAILED: Galacticus model failed to run")
    sys.exit(0)

# Read model results.
with h5py.File("outputs/constrainHaloMassFunction.hdf5", "r") as model:
    output            = model["Outputs/Output1"]
    massHaloModel     = output["haloMass"][:]
    massFunctionModel = output["haloMassFunctionLnMBinAveraged"][:]

# Read target data.
with h5py.File("data/darkMatterHaloMassFunctionMDPL2.hdf5", "r") as target:
    massHaloTarget          = target["massHalo"][:]
    massFunctionTarget      = target["massFunction"][:]
    massFunctionErrorTarget = target["massFunctionError"][:]

emptyBins = np.where(massFunctionTarget <= 0.0)[0]
haloCount = 1.0 / (massFunctionErrorTarget / massFunctionTarget)**2
haloCount[emptyBins] = 0.0

# Check that halo masses agree.
massDifferenceFractional = np.abs(massHaloModel - massHaloTarget) / massHaloTarget
if not np.all(massDifferenceFractional < 1.0e-6):
    print("FAILED: halo masses differ from target dataset")
    sys.exit(0)

# Specify ranges.
massRangeMinimum = 2.5e12
haloCountMinimum = 30

# Select bins within range.
inRange = np.where((massHaloTarget > massRangeMinimum) & (haloCount >= haloCountMinimum))[0]

# Compute normalized offsets.
offsets = np.abs(massFunctionModel[inRange] - massFunctionTarget[inRange]) / massFunctionErrorTarget[inRange]

# Check offsets are sufficiently small.
if np.all(offsets < 4.0):
    print("SUCCESS: halo mass function matches target")
else:
    print("FAILED: halo mass function does not match target")
