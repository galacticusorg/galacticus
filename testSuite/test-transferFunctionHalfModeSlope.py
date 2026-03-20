#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np
import xml.etree.ElementTree as ET

# Test that the "half-mode slope" transfer function agrees with results originated from Daniel Gilman.
# Andrew Benson (ported to Python)

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Read the data from Daniel Gilman's calculation.
dataFile = "data/transferFunctionHalfModeSlopeGilman.txt"
rawData  = np.loadtxt(dataFile)
wavenumberGilman          = rawData[:, 0]
transferFunctionGilman = {
    "-1.0": rawData[:, 1],
    "-2.5": rawData[:, 2],
    "-4.0": rawData[:, 3],
}

# Parse the base parameter file.
tree       = ET.parse("parameters/transferFunctionHalfModeSlope.xml")
root       = tree.getroot()

# Iterate over slopes.
for slope in sorted(transferFunctionGilman.keys()):
    # Set Hubble constant - values in Daniel's file are in "h" units.
    h            = 0.674
    # Set the half-mode mass.
    massHalfMode = 10.0**7.7 / h

    nPoints = len(wavenumberGilman)
    wavenumberMin = wavenumberGilman[0]  * h
    wavenumberMax = wavenumberGilman[-1] * h
    pointsPerDecade = (nPoints - 1) / np.log10(wavenumberGilman[-1] / wavenumberGilman[0])

    # Update XML parameters.
    for elem in root.iter():
        if elem.tag == "task":
            for child in elem:
                if child.tag == "wavenumberMinimum":
                    child.set("value", str(wavenumberMin))
                elif child.tag == "wavenumberMaximum":
                    child.set("value", str(wavenumberMax))
                elif child.tag == "pointsPerDecade":
                    child.set("value", str(pointsPerDecade))
        if elem.tag == "transferFunction":
            for child in elem:
                if child.tag == "massHalfMode":
                    child.set("value", str(massHalfMode))
                elif child.tag == "slopeHalfMode":
                    child.set("value", slope)
        if elem.tag == "outputFileName":
            elem.set("value", "testSuite/outputs/transferFunctionHalfModeSlope.hdf5")

    # Write updated parameter file.
    with open("outputs/transferFunctionHalfModeSlope.xml", "wb") as outFile:
        tree.write(outFile)

    # Run the model.
    status = subprocess.run("cd ..; ./Galacticus.exe testSuite/outputs/transferFunctionHalfModeSlope.xml", shell=True)
    if status.returncode != 0:
        print(f"FAILED: [d\u33c1T/d\u33c1k={slope}] model did not run")
        continue

    # Read the results.
    with h5py.File("outputs/transferFunctionHalfModeSlope.hdf5", "r") as model:
        output           = model["Outputs/Output1"]
        wavenumber       = output["wavenumber"][:]
        transferFunction = output["transferFunction"][:]
    wavenumber /= h

    if len(transferFunction) != len(wavenumberGilman):
        print(f"FAILED: [d\u33c1T/d\u33c1k={slope}] count of wavenumbers does not match")
    else:
        print(f"SUCCESS: [d\u33c1T/d\u33c1k={slope}] count wavenumbers matches")

    if np.any(np.abs(wavenumber - wavenumberGilman) > 1.0e-3 * wavenumberGilman):
        print(f"FAILED: [d\u33c1T/d\u33c1k={slope}] wavenumbers do not match")
    else:
        print(f"SUCCESS: [d\u33c1T/d\u33c1k={slope}] wavenumbers match")

    if np.any(np.abs(transferFunction - transferFunctionGilman[slope]) > 1.0e-3 * transferFunctionGilman[slope]):
        print(f"FAILED: [d\u33c1T/d\u33c1k={slope}] transfer functions do not match")
    else:
        print(f"SUCCESS: [d\u33c1T/d\u33c1k={slope}] transfer functions match")
