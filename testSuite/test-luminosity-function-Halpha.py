#!/usr/bin/env python3
"""Check that the Halpha luminosity function analysis reaches the dust attenuation framework.

The analysis no longer carries a dust normalization of its own. Instead it wraps its emission line extractor in
`nodePropertyExtractorDustAttenuation` and scalarizes the summed, attenuated luminosity:

    scalarizer -> dustAttenuation(outputSumOnly) -> Panuzzo2003

so the dust model given to the analysis is what determines the attenuation. This runs the analysis twice, once with
no dust and once with a screen whose optical depth follows the surface density of metals, and checks that the dust
model actually reaches the luminosities: attenuation must move galaxies to lower luminosity, thinning the bright end
of the luminosity function.

Andrew Benson, with assistance from Claude.
"""
import os
import subprocess
import sys

import h5py
import lxml.etree as etree
import numpy as np

parameterFile = "testSuite/parameters/luminosityFunctionHalpha.xml"
analysisName  = "luminosityFunctionHalphaGunawardhana2013SDSS"

os.makedirs("outputs", exist_ok=True)


def run(label, dustAttenuation):
    """Run the analysis with the given dust model, returning its luminosity function."""
    parameters = etree.parse(os.path.join("..", parameterFile))
    analysis   = parameters.find(".//outputAnalysis[@value='luminosityFunctionGunawardhana2013SDSS']")
    if analysis is None:
        print("FAILED: the Halpha analysis is missing from the parameter file")
        return None
    existing = analysis.find("dustAttenuation")
    if existing is not None:
        analysis.remove(existing)
    analysis.append(dustAttenuation)
    outputFile = f"testSuite/outputs/luminosityFunctionHalpha_{label}.hdf5"
    parameters.find("outputFileName").set("value", outputFile)
    parameterPath = os.path.join("outputs", f"luminosityFunctionHalpha_{label}.xml")
    parameters.write(parameterPath, pretty_print=True)

    outputPath = os.path.join("..", outputFile)
    if os.path.exists(outputPath):
        os.remove(outputPath)
    status = subprocess.run(f"cd ..; ./Galacticus.exe testSuite/{parameterPath}", shell=True, capture_output=True)
    log = status.stdout.decode() + status.stderr.decode()
    if status.returncode != 0 or not os.path.exists(outputPath):
        print(f"FAILED: the Halpha luminosity function model did not run with {label} dust")
        print(log[-2000:])
        return None

    with h5py.File(outputPath, "r") as f:
        if "analyses" not in f or analysisName not in f["analyses"]:
            print(f"FAILED: the {label} run produced no Halpha luminosity function analysis")
            return None
        analysisGroup = f["analyses"][analysisName]
        return analysisGroup["luminosity"][:], analysisGroup["luminosityFunction"][:]


noDust = etree.Element("dustAttenuation", value="zero")

screen = etree.Element("dustAttenuation", value="screenSurfaceDensityMetals")
etree.SubElement(screen, "coefficient", value="1.0")
curve  = etree.SubElement(screen, "dustExtinctionCurve", value="powerLaw")
etree.SubElement(curve, "exponent"           , value="0.7"   )
etree.SubElement(curve, "wavelengthReference", value="5500.0")

unattenuated = run("unattenuated", noDust)
if unattenuated is None:
    sys.exit(0)
attenuated   = run("attenuated"  , screen)
if attenuated is None:
    sys.exit(0)

luminosity, functionUnattenuated = unattenuated
_         , functionAttenuated   = attenuated

if not bool(np.any(functionUnattenuated > 0.0)):
    print("FAILED: the unattenuated luminosity function is zero in every bin, so nothing is being tested")
    sys.exit(0)
print(f"SUCCESS: the analysis is populated in {int((functionUnattenuated > 0.0).sum())} of"
      f" {functionUnattenuated.size} luminosity bins")

# The dust model must actually reach the luminosities.
if bool(np.array_equal(functionAttenuated, functionUnattenuated)):
    print("FAILED: the attenuated and unattenuated luminosity functions are identical, so the dust model is"
          " not reaching the analysis")
    sys.exit(0)

# Dust can only move galaxies to lower luminosity, so the distribution must shift down. Comparing the weighted mean
# of log10(L) rather than individual bins avoids resting the test on a single sparsely populated bin.
def meanLogLuminosity(function):
    return float(np.sum(function * np.log10(luminosity)) / np.sum(function))

meanUnattenuated = meanLogLuminosity(functionUnattenuated)
meanAttenuated   = meanLogLuminosity(functionAttenuated  )
if not meanAttenuated < meanUnattenuated:
    print(f"FAILED: attenuation did not shift the luminosity function to lower luminosity: mean log10(L) went from"
          f" {meanUnattenuated:.4f} to {meanAttenuated:.4f}")
    sys.exit(0)
print(f"SUCCESS: dust shifts the luminosity function to lower luminosity, mean log10(L/erg/s) from"
      f" {meanUnattenuated:.4f} to {meanAttenuated:.4f}")

sys.exit(0)
