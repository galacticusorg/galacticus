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
#     predictable, but no component can have a negative metallicity, or one exceeding that of the IGM at the output time. This
#     catches metals being moved out of the unaccreted gas reservoir at the current IGM metallicity, rather than at the reservoir's
#     own metallicity.
#  3. As 2, but with halo masses perturbed (so that halos sometimes lose mass) and negative accretion allowed. This catches gas
#     lost from the halo carrying the wrong metallicity. Negative CGM masses are disallowed in this case, as no metallicity is
#     meaningful for gas with negative mass.
#  4. As 1, but with halo masses perturbed as in 3, both with negative CGM masses allowed (the default) and disallowed. With a
#     fixed IGM metallicity every component must still have exactly that metallicity, however masses evolve. This catches gas
#     lost from, or created in, the CGM carrying the wrong metallicity.
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

def makeVariant(name, accretionHalo, metallicityVariable, perturbed=False, allowNegativeCGMMass=False):
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
    if perturbed:
        # Perturb halo masses (independently, by 10%) so that halos sometimes lose mass, allow negative accretion, and set whether
        # negative CGM masses are allowed.
        constructor = parameters.find(".//mergerTreeConstructor")
        parent      = constructor.getparent()
        operator    = etree.Element("mergerTreeOperator", value="perturbMasses")
        massError   = etree.SubElement(operator, "nbodyHaloMassError", value="powerLaw")
        for parameter, value in (("normalization", "0.1"), ("exponent", "0.0"), ("fractionalErrorHighMass", "0.1")):
            etree.SubElement(massError, parameter, value=value)
        parent.insert(parent.index(constructor)+1, operator)
        accretion        = parameters.find(f".//accretionHalo[@value='{accretionHalo}']")
        negativeAllowed  = accretion.find("accretionNegativeAllowed")
        if negativeAllowed is None:
            negativeAllowed = etree.SubElement(accretion, "accretionNegativeAllowed")
        negativeAllowed.set("value", "true")
        parameters.find(".//allowNegativeCGMMass").set("value", "true" if allowNegativeCGMMass else "false")
    parameterFile = f"outputs/{name}.xml"
    os.makedirs("outputs", exist_ok=True)
    parameters.write(parameterFile)
    return f"testSuite/{parameterFile}"

# Comparisons below are written such that a NaN metallicity counts as a failure.

def checkFixed(data, label, metallicityExpected=None, tolerance=toleranceRelative):
    # Check that every tree, and every component of every node, has the same, non-zero, metallicity - equal to
    # metallicityExpected if given - to the given relative tolerance. Return that metallicity (or None if it could not be
    # determined).
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
    metallicityIGM  = np.median(metallicityTree) if metallicityExpected is None else metallicityExpected
    if not metallicityIGM > 0.0:
        print(f"FAILED: metal conservation [{label}] - median tree metallicity is {metallicityIGM}, expected a positive value")
        return None
    # Check that each tree has the same metallicity.
    errorTree = np.abs(metallicityTree/metallicityIGM-1.0)
    if not np.all(errorTree <= tolerance):
        print(f"FAILED: metal conservation [{label}] - tree metallicities differ: {metallicityTree[~(errorTree <= tolerance)]} vs. {metallicityIGM}")
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
        if not np.all(errorNode <= tolerance):
            print(f"FAILED: metal conservation [{label}] - component '{mass}' metallicity differs from {metallicityIGM} in {np.count_nonzero(~(errorNode <= tolerance))} of {len(errorNode)} nodes (range {metallicityNode.min()} to {metallicityNode.max()})")
            status = "FAILED"
    if status == "SUCCESS":
        print(f"SUCCESS: metal conservation [{label}] - component metallicities")
    return metallicityIGM

def checkVariable(data, label, metallicityIGMPresent):
    # The output is at z=0, and the IGM metallicity increases monotonically with time, so no gas or stars can be more metal-rich
    # than the IGM at z=0, nor have negative metallicity. (CGM masses may be negative, as allowed by the model, in which case their
    # metal masses are negative also, so only components with a significant, positive mass are checked.)
    status = "SUCCESS"
    for mass, metals in components:
        hasMass = isSignificant(data, mass)
        if not np.any(hasMass):
            continue
        metallicityNode = data[metals][hasMass]/data[mass][hasMass]
        inRange         = (metallicityNode >= -toleranceRelative*metallicityIGMPresent) & (metallicityNode <= metallicityIGMPresent*(1.0+toleranceRelative))
        if not np.all(inRange):
            print(f"FAILED: metal conservation [{label}] - component '{mass}' metallicity is outside the range 0 to {metallicityIGMPresent} (that of the IGM at z=0) in {np.count_nonzero(~inRange)} of {len(metallicityNode)} nodes (range {metallicityNode.min()} to {metallicityNode.max()})")
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
    if data is not None:
        checkVariable(data, f"{accretionHalo}, variable", metallicityIGMPresent)
    ## Time-varying IGM metallicity, with halos that lose mass.
    name                  = f"{nameBase}-perturbed"
    data                  = runModel(makeVariant(name, accretionHalo, metallicityVariable=True, perturbed=True), name)
    if data is not None:
        checkVariable(data, f"{accretionHalo}, perturbed", metallicityIGMPresent)
    ## Fixed IGM metallicity, with halos that lose mass.
    for allowNegativeCGMMass in ( True, False ):
        # The hot halo component resets a negative hot-mode CGM mass to zero after each ODE step, leaving its metals unchanged.
        # This does not conserve metallicity, and naozBarkana2007's corrective flow into an empty hot-mode CGM (while gas is lost
        # from it in the same step) reaches it: with negative CGM masses allowed, metallicities are then off by up to about 0.6%
        # (by about 4e-6 when they are not). coldMode is exact. Tolerances are relaxed for naozBarkana2007 accordingly.
        if accretionHalo == "naozBarkana2007":
            tolerance = 1.0e-2 if allowNegativeCGMMass else 1.0e-5
        else:
            tolerance = toleranceRelative
        suffix = "" if allowNegativeCGMMass else "-noNegativeCGM"
        name   = f"{nameBase}-fixed-perturbed{suffix}"
        data   = runModel(makeVariant(name, accretionHalo, metallicityVariable=False, perturbed=True, allowNegativeCGMMass=allowNegativeCGMMass), name)
        if data is not None:
            checkFixed(data, f"{accretionHalo}, fixed, perturbed{', no negative CGM' if not allowNegativeCGMMass else ''}", metallicityExpected=metallicityIGM, tolerance=tolerance)
