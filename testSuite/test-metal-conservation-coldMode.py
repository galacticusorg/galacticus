#!/usr/bin/env python3
import subprocess
import sys
import os
import glob
import h5py
import numpy as np
from lxml import etree

# Run models with a cold mode CGM, a non-zero IGM metallicity, and zero metal yield to test metal conservation. Each is run with
# the coldMode and naozBarkana2007 halo accretion classes.
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

# Relative tolerance on metallicities.
toleranceRelative = 1.0e-6
# Components holding less than this fraction of their node's halo mass (i.e. about 1e-6 of the halo's baryons) are excluded from
# metallicity checks: their metallicities are set by the ODE solver's absolute tolerances (or by masses which have been depleted
# to subnormal values), not by the physics.
fractionMassMinimum = 1.0e-7

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
        data["basicMass"      ] = nodeData["basicMass"      ][:]
    return data

def isSignificant(data, mass):
    # Return a mask of nodes in which the given component holds a significant mass compared to the node's halo mass.
    return data[mass] > fractionMassMinimum*data["basicMass"]

def makeVariant(name, accretionHalo, metallicityVariable):
    # Build the parameter file for a variant of the base model, and return its path (relative to the Galacticus root).
    parameters = etree.parse("parameters/test-metal-conservation-coldMode.xml")
    parameters.find(".//modelRootDirectory").text = f"testSuite/outputs/{name}"
    if accretionHalo == "naozBarkana2007":
        # Use the naozBarkana2007 accretion class, with the filtering mass and reionization of the reference model. Reionization
        # raises the filtering mass, which drives a corrective flow of gas from the hot halo to the unaccreted reservoir.
        accretion = parameters.find(".//accretionHalo[@value='coldMode']")
        parent    = accretion.getparent()
        index     = parent.index(accretion)
        parent.remove(accretion)
        parent.insert(index  , etree.Element("accretionHalo"                   , value="naozBarkana2007"))
        parent.insert(index+1, etree.Element("intergalacticMediumFilteringMass", value="gnedin2000"     ))
        state        = parameters.find(".//intergalacticMediumState[@value='metallicityFixed']")
        stateFile    = state.find("intergalacticMediumState")
        state.remove(stateFile)
        reionization = etree.SubElement(state, "intergalacticMediumState", value="instantReionization")
        for parameter, value in (("electronScatteringOpticalDepth", "0.0633"), ("reionizationTemperature", "2.0e4"), ("presentDayTemperature", "1.0e3")):
            etree.SubElement(reionization, parameter, value=value)
        reionization.append(stateFile)
    if metallicityVariable:
        state = parameters.find(".//intergalacticMediumState[@value='metallicityFixed']")
        state.set("value", "metallicityPolynomial")
        state.remove(state.find("metallicity"))
        state.insert(0, etree.Element("coefficients", value=" ".join(str(coefficient) for coefficient in metallicityCoefficients)))
    parameterFile = f"outputs/{name}.xml"
    parameters.write(parameterFile)
    return f"testSuite/{parameterFile}"

# Comparisons below are written such that a NaN metallicity counts as a failure.

def checkFixed(data, label):
    # Check that every tree, and every component of every node, has the same, non-zero, metallicity. Return that metallicity (or
    # None if it could not be determined).
    massTotal       = sum(data[mass  ] for mass, metals in components)
    massMetalsTotal = sum(data[metals] for mass, metals in components)
    mergerTreeIndex = data["mergerTreeIndex"]
    trees           = np.unique(mergerTreeIndex)
    massTree        = np.array([massTotal      [mergerTreeIndex == tree].sum() for tree in trees])
    massMetalsTree  = np.array([massMetalsTotal[mergerTreeIndex == tree].sum() for tree in trees])
    hasMass         = massTree > 0.0
    if not np.any(hasMass):
        print(f"FAILED: metal conservation [{label}] - no tree contains any baryonic mass")
        return None
    metallicityTree = massMetalsTree[hasMass]/massTree[hasMass]
    metallicityIGM  = np.median(metallicityTree)
    if not metallicityIGM > 0.0:
        print(f"FAILED: metal conservation [{label}] - median tree metallicity is {metallicityIGM}, expected a positive value")
        return None
    # Check that each tree has the same metallicity.
    errorTree = np.abs(metallicityTree/metallicityIGM-1.0)
    if not np.all(errorTree <= toleranceRelative):
        print(f"FAILED: metal conservation [{label}] - tree metallicities differ: {metallicityTree[~(errorTree <= toleranceRelative)]} vs. {metallicityIGM}")
    else:
        print(f"SUCCESS: metal conservation [{label}] - tree metallicities")
    # Check that each component of each node has the same metallicity.
    status = "SUCCESS"
    for mass, metals in components:
        hasMass = isSignificant(data, mass)
        if not np.any(hasMass):
            continue
        metallicityNode = data[metals][hasMass]/data[mass][hasMass]
        errorNode       = np.abs(metallicityNode/metallicityIGM-1.0)
        if not np.all(errorNode <= toleranceRelative):
            print(f"FAILED: metal conservation [{label}] - component '{mass}' metallicity differs from {metallicityIGM} in {np.count_nonzero(~(errorNode <= toleranceRelative))} of {len(errorNode)} nodes (range {metallicityNode.min()} to {metallicityNode.max()})")
            status = "FAILED"
    if status == "SUCCESS":
        print(f"SUCCESS: metal conservation [{label}] - component metallicities")
    return metallicityIGM

def checkVariable(data, label, metallicityIGMPresent):
    # The output is at z=0, and the IGM metallicity increases monotonically with time, so no gas or stars can be more metal-rich
    # than the IGM at z=0.
    status = "SUCCESS"
    for mass, metals in components:
        if not np.all(data[metals] >= 0.0):
            print(f"FAILED: metal conservation [{label}] - component '{mass}' has negative metal mass in {np.count_nonzero(~(data[metals] >= 0.0))} of {len(data[metals])} nodes (minimum {data[metals].min()})")
            status = "FAILED"
        hasMass = isSignificant(data, mass)
        if not np.any(hasMass):
            continue
        metallicityNode = data[metals][hasMass]/data[mass][hasMass]
        if not np.all(metallicityNode <= metallicityIGMPresent*(1.0+toleranceRelative)):
            print(f"FAILED: metal conservation [{label}] - component '{mass}' metallicity exceeds that of the IGM at z=0 ({metallicityIGMPresent}) in {np.count_nonzero(~(metallicityNode <= metallicityIGMPresent*(1.0+toleranceRelative)))} of {len(metallicityNode)} nodes (maximum {metallicityNode.max()})")
            status = "FAILED"
    # Every unaccreted reservoir received gas from the IGM at a non-zero metallicity, so must be enriched.
    hasMass = isSignificant(data, "hotHaloUnaccretedMass")
    if not np.all(data["hotHaloUnaccretedAbundancesMetals"][hasMass] > 0.0):
        print(f"FAILED: metal conservation [{label}] - {np.count_nonzero(~(data['hotHaloUnaccretedAbundancesMetals'][hasMass] > 0.0))} unaccreted gas reservoirs contain no metals")
        status = "FAILED"
    if status == "SUCCESS":
        print(f"SUCCESS: metal conservation [{label}] - component metallicities")

for accretionHalo in ( "coldMode", "naozBarkana2007" ):
    ## Fixed IGM metallicity.
    nameBase = "test-metal-conservation-coldMode"+("" if accretionHalo == "coldMode" else f"-{accretionHalo}")
    name     = nameBase
    if accretionHalo == "coldMode":
        parameterFile = "testSuite/parameters/test-metal-conservation-coldMode.xml"
    else:
        parameterFile = makeVariant(name, accretionHalo, metallicityVariable=False)
    data = runModel(parameterFile, name)
    if data is None:
        continue
    metallicityIGM = checkFixed(data, accretionHalo)
    if metallicityIGM is None:
        continue
    ## Time-varying IGM metallicity.
    # The Solar metallicity used by the model follows from the fixed-metallicity run.
    metallicitySolar      = metallicityIGM/metallicityFixed
    metallicityIGMPresent = metallicitySolar*10.0**metallicityCoefficients[0]
    name                  = f"{nameBase}-variable"
    data                  = runModel(makeVariant(name, accretionHalo, metallicityVariable=True), name)
    if data is None:
        continue
    checkVariable(data, f"{accretionHalo}, variable", metallicityIGMPresent)
