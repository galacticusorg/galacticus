#!/usr/bin/env python3
import subprocess
import sys
import os
import glob
import h5py
import numpy as np

# Run a cold mode model with a fixed, non-zero IGM metallicity and zero metal yield to test metal conservation. All baryons
# enter the halo from the IGM with the same metallicity, and no metals are produced, so every tree, and every gas and stellar
# component of every node, should end up with exactly the IGM metallicity (to rounding error).
# Andrew Benson, Claude

# Run the model.
subprocess.run(
    "cd ..; ./scripts/aux/launch.py testSuite/parameters/test-metal-conservation-coldMode.xml --launchMethod local --threadMaximum 1 --ompThreads 4",
    shell=True
)

# Check for failed models.
logFiles = glob.glob("outputs/test-metal-conservation-coldMode/galacticus_*/galacticus.log")
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

fileName = "outputs/test-metal-conservation-coldMode/galacticus_0:1/galacticus.hdf5"
if not os.path.exists(fileName):
    print(f"FAILED: model output '{fileName}' was not found")
    sys.exit(0)

# Pairs of (mass, metal mass) properties for all gas and stellar components. Black holes carry no metals, so are excluded.
components = (
    ("diskMassGas"          , "diskAbundancesGasMetals"          ),
    ("diskMassStellar"      , "diskAbundancesStellarMetals"      ),
    ("spheroidMassGas"      , "spheroidAbundancesGasMetals"      ),
    ("spheroidMassStellar"  , "spheroidAbundancesStellarMetals"  ),
    ("hotHaloMass"          , "hotHaloAbundancesMetals"          ),
    ("hotHaloMassCold"      , "hotHaloAbundancesColdMetals"      ),
    ("hotHaloOutflowedMass" , "hotHaloOutflowedAbundancesMetals" ),
    ("hotHaloUnaccretedMass", "hotHaloUnaccretedAbundancesMetals"),
)
with h5py.File(fileName, "r") as galacticus:
    nodeData        = galacticus["Outputs/Output1/nodeData"]
    mergerTreeIndex = nodeData["mergerTreeIndex"][:]
    masses          = {mass  : nodeData[mass ][:] for mass, metals in components}
    massesMetals    = {metals: nodeData[metals][:] for mass, metals in components}

# Find the metallicity of each tree.
massTotal       = sum(masses      [mass  ] for mass, metals in components)
massMetalsTotal = sum(massesMetals[metals] for mass, metals in components)
trees           = np.unique(mergerTreeIndex)
massTree        = np.array([massTotal      [mergerTreeIndex == tree].sum() for tree in trees])
massMetalsTree  = np.array([massMetalsTotal[mergerTreeIndex == tree].sum() for tree in trees])
hasMass         = massTree > 0.0
if not np.any(hasMass):
    print("FAILED: metal conservation [coldMode] - no tree contains any baryonic mass")
    sys.exit(0)
metallicityTree = massMetalsTree[hasMass]/massTree[hasMass]
metallicityIGM  = np.median(metallicityTree)
# Comparisons are written such that a NaN metallicity counts as a failure.
if not metallicityIGM > 0.0:
    print(f"FAILED: metal conservation [coldMode] - median tree metallicity is {metallicityIGM}, expected a positive value")
    sys.exit(0)

# Check that each tree has the same metallicity.
toleranceRelative = 1.0e-6
errorTree         = np.abs(metallicityTree/metallicityIGM-1.0)
if not np.all(errorTree <= toleranceRelative):
    print(f"FAILED: metal conservation [coldMode] - tree metallicities differ: {metallicityTree[~(errorTree <= toleranceRelative)]} vs. {metallicityIGM}")
else:
    print("SUCCESS: metal conservation [coldMode] - tree metallicities")

# Check that each component of each node has the same metallicity.
status = "SUCCESS"
for mass, metals in components:
    hasMass = masses[mass] > 0.0
    if not np.any(hasMass):
        continue
    metallicityNode = massesMetals[metals][hasMass]/masses[mass][hasMass]
    errorNode       = np.abs(metallicityNode/metallicityIGM-1.0)
    if not np.all(errorNode <= toleranceRelative):
        print(f"FAILED: metal conservation [coldMode] - component '{mass}' metallicity differs from {metallicityIGM} in {np.count_nonzero(~(errorNode <= toleranceRelative))} of {len(errorNode)} nodes (range {metallicityNode.min()} to {metallicityNode.max()})")
        status = "FAILED"
if status == "SUCCESS":
    print("SUCCESS: metal conservation [coldMode] - component metallicities")
