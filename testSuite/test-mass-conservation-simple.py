#!/usr/bin/env python3
import subprocess
import sys
import glob
import h5py
import numpy as np

# Run a simple model to test mass conservation.
# Andrew Benson (ported to Python)

# Run the model.
subprocess.run(
    "cd ..; scripts/aux/launch.pl testSuite/parameters/test-mass-conservation-simple.xml --launchMethod local --threadMaximum 1 --ompThreads 4",
    shell=True
)

# Check for failed models.
logFiles = glob.glob("outputs/test-mass-conservation-simple/galacticus_*/galacticus.log")
failures = []
for logFile in logFiles:
    result = subprocess.run(f"grep -q -i -e fatal -e aborted {logFile}", shell=True)
    if result.returncode == 0:
        failures.append(logFile)

if failures:
    for failure in failures:
        print(f"FAILED: log from {failure}:")
        with open(failure) as f:
            print(f.read())
else:
    print("SUCCESS: model run")

# Extract masses and check they add up to the expected value.
with h5py.File("outputs/test-mass-conservation-simple/galacticus_0:1/galacticus.hdf5", "r") as galacticus:
    # Read parameters.
    cosmologyParams = galacticus["Parameters/cosmologyParameters"]
    OmegaBaryon     = cosmologyParams.attrs["OmegaBaryon"]
    OmegaMatter     = cosmologyParams.attrs["OmegaMatter"]
    baryonFraction  = OmegaBaryon / OmegaMatter

    # Read datasets.
    nodeData = galacticus["Outputs/Output1/nodeData"]
    props    = {}
    for key in ("mergerTreeWeight", "diskMassStellar", "diskMassGas", "hotHaloMass",
                "hotHaloOutflowedMass", "nodeIsIsolated", "massBertschinger",
                "hotHaloUnaccretedMass", "mergerTreeIndex"):
        props[key] = nodeData[key][:]

# Find centrals and satellites.
centrals   = np.where(props["nodeIsIsolated"] == 1)[0]
satellites = np.where(props["nodeIsIsolated"] == 0)[0]

# Find masses in satellites.
massSatellites = np.zeros(len(centrals))
for i, ci in enumerate(centrals):
    treeIdx = props["mergerTreeIndex"][ci]
    inTree  = np.where(props["mergerTreeIndex"][satellites] == treeIdx)[0]
    satMass = (
        props["diskMassStellar"][satellites][inTree].sum()
        + props["diskMassGas"][satellites][inTree].sum()
        + props["hotHaloMass"][satellites][inTree].sum()
        + props["hotHaloOutflowedMass"][satellites][inTree].sum()
        + props["hotHaloUnaccretedMass"][satellites][inTree].sum()
    )
    massSatellites[i] = satMass / baryonFraction / props["massBertschinger"][ci]

# Normalize central masses.
normFactor = baryonFraction * props["massBertschinger"]
massTotal = (
      props["diskMassStellar"      ][centrals] / normFactor[centrals]
    + props["diskMassGas"          ][centrals] / normFactor[centrals]
    + props["hotHaloMass"          ][centrals] / normFactor[centrals]
    + props["hotHaloOutflowedMass" ][centrals] / normFactor[centrals]
    + props["hotHaloUnaccretedMass"][centrals] / normFactor[centrals]
    + massSatellites
)

# Check that all masses are unity.
if np.any(np.abs(massTotal - 1.0) > 1.0e-6):
    print(f"FAILED: mass conservation [simple] failure -> {massTotal - 1.0}")
else:
    print("SUCCESS: mass conservation")
