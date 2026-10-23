#!/usr/bin/env python3
import subprocess
import sys
import os
import glob
import h5py
import numpy as np
from lxml import etree

# Run cold mode models with a non-zero IGM metallicity and zero metal yield to test metal conservation.
#  1. With a fixed IGM metallicity, all baryons enter the halo from the IGM with the same metallicity, and no metals are produced,
#     so every tree, and every gas and stellar component of every node, should end up with exactly the IGM metallicity (to
#     rounding error).
#  2. With an IGM metallicity that increases with time (the polynomial used by the reference model) exact metallicities are not
#     predictable, but no component can have a negative metal mass, or a metallicity exceeding that of the IGM at the output
#     time. This catches metals being moved out of the unaccreted gas reservoir at the current IGM metallicity, rather than at the
#     reservoir's own metallicity.
# Andrew Benson, Claude

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

# Fixed IGM metallicity (in Solar units) of the base model, and polynomial coefficients for the time-varying variant.
metallicityFixed        = 0.1
metallicityCoefficients = (-1.3, -1.9)

def runModel(parameterFile, name):
    # Run a model, check its logs for failures, and return its node data (or None if the model failed).
    subprocess.run(
        f"cd ..; ./scripts/aux/launch.py {parameterFile} --launchMethod local --threadMaximum 1 --ompThreads 4",
        shell=True
    )
    logFiles = glob.glob(f"outputs/{name}/galacticus_*/galacticus.log")
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
        print(f"SUCCESS: model run [{name}]")
    fileName = f"outputs/{name}/galacticus_0:1/galacticus.hdf5"
    if not os.path.exists(fileName):
        print(f"FAILED: model output '{fileName}' was not found")
        return None
    with h5py.File(fileName, "r") as galacticus:
        nodeData = galacticus["Outputs/Output1/nodeData"]
        data     = {property: nodeData[property][:] for pair in components for property in pair}
        data["mergerTreeIndex"] = nodeData["mergerTreeIndex"][:]
    return data

# Comparisons below are written such that a NaN metallicity counts as a failure.

## Fixed IGM metallicity.
data = runModel("testSuite/parameters/test-metal-conservation-coldMode.xml", "test-metal-conservation-coldMode")
if data is None:
    sys.exit(0)

# Find the metallicity of each tree.
massTotal       = sum(data[mass  ] for mass, metals in components)
massMetalsTotal = sum(data[metals] for mass, metals in components)
mergerTreeIndex = data["mergerTreeIndex"]
trees           = np.unique(mergerTreeIndex)
massTree        = np.array([massTotal      [mergerTreeIndex == tree].sum() for tree in trees])
massMetalsTree  = np.array([massMetalsTotal[mergerTreeIndex == tree].sum() for tree in trees])
hasMass         = massTree > 0.0
if not np.any(hasMass):
    print("FAILED: metal conservation [coldMode] - no tree contains any baryonic mass")
    sys.exit(0)
metallicityTree = massMetalsTree[hasMass]/massTree[hasMass]
metallicityIGM  = np.median(metallicityTree)
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
    hasMass = data[mass] > 0.0
    if not np.any(hasMass):
        continue
    metallicityNode = data[metals][hasMass]/data[mass][hasMass]
    errorNode       = np.abs(metallicityNode/metallicityIGM-1.0)
    if not np.all(errorNode <= toleranceRelative):
        print(f"FAILED: metal conservation [coldMode] - component '{mass}' metallicity differs from {metallicityIGM} in {np.count_nonzero(~(errorNode <= toleranceRelative))} of {len(errorNode)} nodes (range {metallicityNode.min()} to {metallicityNode.max()})")
        status = "FAILED"
if status == "SUCCESS":
    print("SUCCESS: metal conservation [coldMode] - component metallicities")

## Time-varying IGM metallicity.
# The Solar metallicity used by the model follows from the fixed-metallicity run.
metallicitySolar = metallicityIGM/metallicityFixed
# Build the parameter file for this variant from that of the fixed-metallicity model.
name       = "test-metal-conservation-coldMode-variable"
parameters = etree.parse("parameters/test-metal-conservation-coldMode.xml")
parameters.find(".//modelRootDirectory").text = f"testSuite/outputs/{name}"
state = parameters.find(".//intergalacticMediumState[@value='metallicityFixed']")
state.set("value", "metallicityPolynomial")
state.remove(state.find("metallicity"))
coefficients = etree.Element("coefficients", value=" ".join(str(coefficient) for coefficient in metallicityCoefficients))
state.insert(0, coefficients)
parameterFile = f"outputs/{name}.xml"
parameters.write(parameterFile)
data = runModel(f"testSuite/{parameterFile}", name)
if data is None:
    sys.exit(0)

# The output is at z=0, and the IGM metallicity increases monotonically with time, so no gas or stars can be more metal-rich than
# the IGM at z=0.
metallicityIGMPresent = metallicitySolar*10.0**metallicityCoefficients[0]
status = "SUCCESS"
for mass, metals in components:
    if not np.all(data[metals] >= 0.0):
        print(f"FAILED: metal conservation [coldMode, variable] - component '{mass}' has negative metal mass in {np.count_nonzero(~(data[metals] >= 0.0))} of {len(data[metals])} nodes (minimum {data[metals].min()})")
        status = "FAILED"
    hasMass = data[mass] > 0.0
    if not np.any(hasMass):
        continue
    metallicityNode = data[metals][hasMass]/data[mass][hasMass]
    if not np.all(metallicityNode <= metallicityIGMPresent*(1.0+toleranceRelative)):
        print(f"FAILED: metal conservation [coldMode, variable] - component '{mass}' metallicity exceeds that of the IGM at z=0 ({metallicityIGMPresent}) in {np.count_nonzero(~(metallicityNode <= metallicityIGMPresent*(1.0+toleranceRelative)))} of {len(metallicityNode)} nodes (maximum {metallicityNode.max()})")
        status = "FAILED"
# Every unaccreted reservoir received gas from the IGM at a non-zero metallicity, so must be enriched.
hasMass = data["hotHaloUnaccretedMass"] > 0.0
if not np.all(data["hotHaloUnaccretedAbundancesMetals"][hasMass] > 0.0):
    print(f"FAILED: metal conservation [coldMode, variable] - {np.count_nonzero(~(data['hotHaloUnaccretedAbundancesMetals'][hasMass] > 0.0))} unaccreted gas reservoirs contain no metals")
    status = "FAILED"
if status == "SUCCESS":
    print("SUCCESS: metal conservation [coldMode, variable] - component metallicities")
