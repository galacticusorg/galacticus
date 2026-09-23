#!/usr/bin/env python3
"""Test the options of the Sobral et al. (2013) HiZELS H-alpha luminosity function analysis used for emulator calibration.

Runs the model of testSuite/parameters/luminosityFunctionHalpha.xml, with its outputs moved to span the first HiZELS
redshift interval (z=0.40), through the `luminosityFunctionSobral2013HiZELS` analysis with the
`stellarMassRedshift` dust attenuation, three times:

* with the mean attenuation only;
* with emission lines from active galactic nuclei added (`nodePropertyExtractorAGN`);
* with scatter in the attenuation (`rootVarianceAttenuation`).

It checks that each run populates the analysis, that adding AGN emission changes the luminosity function and can only
increase the number of galaxies brighter than any luminosity (adding luminosity to every galaxy, with the same
attenuation and observational error, can move none of them fainter---whereas the mean luminosity within the range of
the analysis need not increase, as galaxies are brought into range at its faint end), and that scatter in the attenuation changes and broadens the luminosity function. (The number of
galaxies within the range of the analysis is not conserved by the scatter---the truncated attenuation can only make a
galaxy fainter than it is unattenuated, so galaxies are scattered out of the faint end---so conservation is tested
instead, exactly, by the tests.output_analyses.attenuation_scatter unit test.)

Black holes in the base model accrete too little for their H-alpha emission (at most ~1e38 erg/s) to move any galaxy
between luminosity bins, so Bondi-Hoyle accretion is made stronger (and not restricted to the hot mode), which makes it
significant. This is applied to all three runs, which therefore differ only in the options of the analysis.

Andrew Benson, Claude (23-September-2026).
"""

import os
import subprocess
import sys

import h5py
import lxml.etree as etree
import numpy as np

parameterFile = "testSuite/parameters/luminosityFunctionHalpha.xml"
analysisName  = "luminosityFunctionHalphaSobral2013HiZELSZ1"

os.makedirs("outputs", exist_ok=True)


def reportFailure(message, log):
    """Print a failure together with the part of the model log which explains it."""
    print(f"FAILED: {message}")
    for marker in ("Fatal error", "FATAL", "Error occurred", "experienced a segfault"):
        index = log.find(marker)
        if index >= 0:
            print(log[max(0, index - 500):index + 1500])
            return
    print(log[-2000:])


def analysis(rootVarianceAttenuation=0.0, includeAGN=False):
    """Return the Sobral (2013) HiZELS analysis element, with the given options."""
    element = etree.Element("outputAnalysis", value="luminosityFunctionSobral2013HiZELS")
    etree.SubElement(element, "redshiftInterval"                    , value="1"     )
    etree.SubElement(element, "randomErrorMinimum"                  , value="+0.1"  )
    etree.SubElement(element, "randomErrorMaximum"                  , value="+0.1"  )
    etree.SubElement(element, "randomErrorPolynomialCoefficient"    , value="+0.1"  )
    etree.SubElement(element, "systematicErrorPolynomialCoefficient", value="0.0"   )
    etree.SubElement(element, "sizeSourceLensing"                   , value="2.0e-3")
    etree.SubElement(element, "rootVarianceAttenuation"             , value=str(rootVarianceAttenuation))
    dust = etree.SubElement(element, "dustAttenuation", value="stellarMassRedshift")
    etree.SubElement(dust, "dustExtinctionCurve", value="calzetti2000")
    if includeAGN:
        agn = etree.SubElement(element, "nodePropertyExtractorAGN", value="luminosityEmissionLineAGN")
        etree.SubElement(agn, "lineNames", value="balmerAlpha6565")
    return element


def run(label, analysisElement):
    """Run the model with the given analysis, returning its luminosities and luminosity function."""
    parameters = etree.parse(os.path.join("..", parameterFile))
    # Replace the analyses with the one under test.
    analyses = parameters.find("outputAnalysis")
    for child in list(analyses):
        analyses.remove(child)
    analyses.append(analysisElement)
    # Outputs spanning the first HiZELS redshift interval (z=0.401+/-0.010).
    parameters.find("outputTimes/redshifts").set("value", "0.38 0.39 0.40 0.41 0.42")
    # More strongly accreting black holes, so that AGN emission is significant.
    accretion = parameters.find("blackHoleAccretionRate")
    accretion.find("bondiHoyleAccretionEnhancementSpheroid").set("value", "500.0")
    accretion.find("bondiHoyleAccretionEnhancementHotHalo" ).set("value", "600.0")
    accretion.find("bondiHoyleAccretionHotModeOnly"        ).set("value", "false")
    outputFile = f"testSuite/outputs/luminosityFunctionHalphaSobral_{label}.hdf5"
    parameters.find("outputFileName").set("value", outputFile)
    parameterPath = os.path.join("outputs", f"luminosityFunctionHalphaSobral_{label}.xml")
    parameters.write(parameterPath, pretty_print=True)
    outputPath = os.path.join("..", outputFile)
    if os.path.exists(outputPath):
        os.remove(outputPath)
    status = subprocess.run(f"cd ..; ./Galacticus.exe testSuite/{parameterPath}", shell=True, capture_output=True)
    log = status.stdout.decode() + status.stderr.decode()
    if status.returncode != 0 or not os.path.exists(outputPath):
        reportFailure(f"the {label} model did not run", log)
        return None
    with h5py.File(outputPath, "r") as f:
        if "analyses" not in f or analysisName not in f["analyses"]:
            print(f"FAILED: the {label} run produced no Sobral (2013) H-alpha luminosity function analysis")
            return None
        group = f["analyses"][analysisName]
        return group["luminosity"][:], group["luminosityFunction"][:]


def moments(luminosity, function):
    """Return the total, and the weighted mean and standard deviation of log10(L), of a luminosity function."""
    logL  = np.log10(luminosity)
    total = float(np.sum(function))
    mean  = float(np.sum(function * logL) / total)
    sigma = float(np.sqrt(np.sum(function * (logL - mean) ** 2) / total))
    return total, mean, sigma


def main():
    results = {}
    for label, element in (("mean"   , analysis()                            ),
                           ("AGN"    , analysis(includeAGN=True)             ),
                           ("scatter", analysis(rootVarianceAttenuation=0.5))):
        result = run(label, element)
        if result is None:
            return
        results[label] = result
    luminosity = results["mean"][0]
    for label, (_, function) in results.items():
        if not bool(np.any(function > 0.0)):
            print(f"FAILED: the {label} luminosity function is zero in every bin, so nothing is being tested")
            return
    print("SUCCESS: all three runs populate the analysis")
    _, _, sigma = moments(luminosity, results["mean"][1])
    # AGN emission adds luminosity to every galaxy, so must change the luminosity function, and can only increase the number
    # of galaxies brighter than each luminosity.
    brighter    = np.cumsum(results["mean"][1][::-1])[::-1]
    brighterAGN = np.cumsum(results["AGN" ][1][::-1])[::-1]
    if bool(np.array_equal(results["AGN"][1], results["mean"][1])):
        print("FAILED: adding AGN emission did not change the luminosity function")
    elif bool(np.any(brighterAGN < brighter * (1.0 - 1.0e-9))):
        print("FAILED: adding AGN emission reduced the number of galaxies brighter than some luminosity")
    else:
        print(f"SUCCESS: AGN emission increases the number of galaxies brighter than each luminosity (by up to a factor {np.max(brighterAGN[brighter > 0.0] / brighter[brighter > 0.0]):.2f})")
    # Scatter in the attenuation must change, and broaden, the luminosity function.
    _, _, sigmaScatter = moments(luminosity, results["scatter"][1])
    if bool(np.array_equal(results["scatter"][1], results["mean"][1])):
        print("FAILED: scatter in the attenuation did not change the luminosity function")
    elif not sigmaScatter > sigma:
        print(f"FAILED: scatter in the attenuation did not broaden the luminosity function: width went from {sigma:.4f} to {sigmaScatter:.4f} dex")
    else:
        print(f"SUCCESS: scatter in the attenuation broadens the luminosity function, width from {sigma:.4f} to {sigmaScatter:.4f} dex")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # Report any unexpected error as a failure, rather than a bare traceback.
        print(f"FAILED: unexpected error: {type(error).__name__}: {error}")
    sys.exit(0)
