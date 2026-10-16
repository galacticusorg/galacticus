#!/usr/bin/env python3
"""Check the SDSS analyses which reach the dust attenuation framework.

`outputAnalysisColorDistributionSDSS` and `outputAnalysisLuminosityFunctionMonteroDorta2009SDSS` both used to build
`nodePropertyExtractorLmnstyStllrCF2000`, which applied Charlot & Fall (2000) dust internally. They now build

    scalarizer -> dustAttenuation(outputSumOnly) -> luminosityStellar(disk), luminosityStellar(spheroid)

Note this changed their results. The retired extractor evaluated its extinction curve at the *observed* filter
wavelength; the framework evaluates it in the rest frame of the emitting galaxy, which is where the dust is. Both
analyses use observed-frame SDSS filters at a band redshift of 0.1, so their optical depths increased by a factor
(1.1)^0.7 = 1.0693.

What is checked here is that both analyses construct and run, and that the luminosity function they produce is
populated and finite. The color distribution is checked only for finiteness: the small model run here populates no
galaxies in the stellar mass and redshift bin that distribution 1 selects, so its function is legitimately zero.

Andrew Benson, with assistance from Claude.
"""
import os
import subprocess
import sys

import h5py
import numpy as np


def reportFailure(message, log):
    """Print a failure together with the part of the model log which explains it.

    The tail of a Galacticus log is the bottom of a backtrace, which says where the code was rather than what went
    wrong. The message itself sits above that, so print from it where one can be found.
    """
    print(f"FAILED: {message}")
    for marker in ("Fatal error", "FATAL", "Error occurred", "experienced a segfault"):
        index = log.find(marker)
        if index >= 0:
            print(log[max(0, index - 500):index + 1500])
            return
    print(log[-2000:])

parameterFile = "testSuite/parameters/dustAnalysesSDSS.xml"
outputFile    = "testSuite/outputs/dustAnalysesSDSS.hdf5"
outputPath    = os.path.join("..", outputFile)

os.makedirs("outputs", exist_ok=True)
if os.path.exists(outputPath):
    os.remove(outputPath)

status = subprocess.run(f"cd ..; ./Galacticus.exe {parameterFile}", shell=True, capture_output=True)
log = status.stdout.decode() + status.stderr.decode()
if status.returncode != 0 or not os.path.exists(outputPath):
    reportFailure("the SDSS dust analyses model did not run", log)
    sys.exit(0)

with h5py.File(outputPath, "r") as f:
    if "analyses" not in f:
        print("FAILED: the run produced no analyses")
        sys.exit(0)
    analyses = f["analyses"]
    for name in ("colorDistributionSDSS01", "luminosityFunctionMonteroDorta2009SDSSr"):
        if name not in analyses:
            print(f"FAILED: analysis '{name}' is missing from the output: {list(analyses.keys())}")
            sys.exit(0)
    color             = analyses["colorDistributionSDSS01"]["colorDistributionSDSSFunction"][:]
    luminosityFunction = analyses["luminosityFunctionMonteroDorta2009SDSSr"]["luminosityFunction"][:]

print("SUCCESS: both SDSS analyses construct and run through the dust attenuation framework")

for name, values in (("color distribution", color), ("luminosity function", luminosityFunction)):
    if not bool(np.all(np.isfinite(values))):
        print(f"FAILED: the {name} contains non-finite values")
        sys.exit(0)
    if bool(np.any(values < 0.0)):
        print(f"FAILED: the {name} contains negative values")
        sys.exit(0)
print("SUCCESS: both analyses produce finite, non-negative distributions")

# The luminosity function must actually be populated -- this is what exercises the attenuated luminosity end to end.
countPopulated = int((luminosityFunction > 0.0).sum())
if countPopulated == 0:
    print("FAILED: the luminosity function is zero in every bin, so the attenuated luminosity is not being tested")
    sys.exit(0)
print(f"SUCCESS: the luminosity function is populated in {countPopulated} of {luminosityFunction.size} bins")

sys.exit(0)
