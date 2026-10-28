#!/usr/bin/env python3
"""Check broadband luminosities computed directly from the star formation history, for lightcone output.

The same luminosities are available by two routes, and this test compares them. The reference route integrates a filter response
over an SED (`nodePropertyExtractorLuminosityStellarFromSED`): it needs no tabulation in redshift and performs the filter
integral itself, so it is exact up to the resolution of the SED. The route under test contracts the star formation history
against a template of filter-convolved simple stellar population luminosities
(`nodePropertyExtractorLuminosityStellarFromStarFormationHistory`), which for observed-frame luminosities must be tabulated on a
grid of redshifts and interpolated to the redshift of each galaxy.

Three things are checked:

 * that the two routes agree, in both frames, for every filter;
 * that a deliberately coarse grid of redshifts agrees *less* well, so that the comparison above is known to be measuring the
   interpolation in redshift rather than reporting that two routes happen to agree;
 * that the two conditions the direct route requires are reported as errors rather than silently violated: a star formation
   history which does not use fixed ages, and a galaxy which lies outside the tabulated range of redshift.

Note that the agreement of the two routes is the substantive check here, and is stronger than a comparison against an
independent integration of the output SED would be: such an integration shares the metallicity and age averaging of the SED
route, and so cannot detect an error in that averaging, whereas the two routes perform it independently.
"""

import os
import re
import subprocess
import sys

import h5py
import numpy as np

# Tolerance on the agreement between the two routes. A residual difference remains however fine the grid of redshifts, since the
# reference route integrates the filter over an SED of finite resolution while the route under test uses the exact filter
# integral tabulated for each filter. Measured at 3.7e-3 over all filters and both frames.
TOLERANCE_ROUTES      = 5.0e-3
# A coarse grid of redshifts must disagree by at least this factor more than the default grid, and by at least this absolute
# amount, for the comparison above to be demonstrably sensitive to the interpolation. Measured: the coarse grid of the companion
# model is in error by 2e-2, against 2e-4 for the default grid.
FACTOR_COARSE_MINIMUM = 5.0
COARSE_MINIMUM        = 3.0e-3

# The filters, and the frames, for which both routes are computed.
FILTERS = ( "Galex_FUV", "Galex_NUV", "SDSS_r", "adaptiveResolutionTopHat_1300_200" )
FRAMES  = ( "observed", "rest" )
# The filter used for the coarse-grid comparison. At the redshifts of this lightcone it samples the spectrum across the Lyman
# break, which makes it the filter most sensitive to the grid.
FILTER_COARSE = "adaptiveResolutionTopHat_1300_200"

# This script runs with the working directory set to testSuite/, while the model is run from the repository root, so the two
# need different paths to the same file.
outputDirectory = "outputs/test-luminosities-from-SFH-lightcone"
outputPath      = os.path.join(outputDirectory,"galacticus.hdf5"      )
coarsePath      = os.path.join(outputDirectory,"galacticusCoarse.hdf5")
agesPath        = os.path.join(outputDirectory,"galacticusAges.hdf5"  )
rangePath       = os.path.join(outputDirectory,"galacticusRange.hdf5" )
os.makedirs(outputDirectory,exist_ok=True)
for path in ( outputPath, coarsePath, agesPath, rangePath ):
    if os.path.exists(path):
        os.remove(path)

def runModel(arguments,logName):
    """Run a model, returning its exit status and log."""
    status = subprocess.run("cd ..; ./Galacticus.exe "+arguments,shell=True,capture_output=True,stdin=subprocess.DEVNULL)
    log    = status.stdout.decode()+status.stderr.decode()
    with open(os.path.join(outputDirectory,logName),"w") as logFile:
        logFile.write(log)
    return status.returncode, log

def ratios(nodes,filter,frame):
    """Return the ratio of the two routes' luminosities, for the galaxies which both find luminous."""
    direct    = nodes[f"diskStellarSFH:identity:{filter}:{frame}"][:]
    reference = nodes[f"diskStellarSED:identity:{filter}:{frame}"][:]
    luminous  = ( direct > 0.0 ) & ( reference > 0.0 )
    return luminous, direct[luminous]/reference[luminous]

print("Running model...")
status, log = runModel("testSuite/parameters/test-luminosities-from-SFH-lightcone.xml","galacticus.log")
if status != 0 or re.search(r"fatal|aborted|Galacticus experienced an error in the GSL library",log,re.IGNORECASE) or not os.path.exists(outputPath):
    print("FAILED: model run:")
    print(log[-4000:])
    sys.exit(0)
print("SUCCESS: model run")

# Compare the two routes.
with h5py.File(outputPath,"r") as model:
    try:
        nodes    = model["Lightcone/Output1/nodeData"]
        redshift = nodes["redshift"][:]
        # Check that the lightcone holds galaxies which both routes find luminous, as otherwise nothing below is tested.
        luminous, _ = ratios(nodes,"SDSS_r","rest")
        countLuminous = int(luminous.sum())
        if countLuminous == 0:
            print("FAILED: no galaxy in the lightcone is luminous by both routes, so nothing is being tested")
            sys.exit(0)
        print(f"SUCCESS: {countLuminous} of {len(redshift)} galaxies in the lightcone are luminous, spanning redshifts {redshift[luminous].min():.3f} to {redshift[luminous].max():.3f}")
        differenceFine = {}
        for frame in FRAMES:
            for filter in FILTERS:
                _, ratio   = ratios(nodes,filter,frame)
                difference = np.abs(ratio-1.0)
                if filter == FILTER_COARSE and frame == "observed":
                    differenceFine[filter] = np.median(difference)
                if np.any(~np.isfinite(difference)) or np.any(difference > TOLERANCE_ROUTES):
                    indexWorst = np.nanargmax(np.where(np.isfinite(difference),difference,np.inf))
                    print(f"FAILED: {filter}:{frame} computed directly from the star formation history against the from-SED route: maximum relative difference is {difference[indexWorst]:.3e} (tolerance {TOLERANCE_ROUTES:.1e})")
                else:
                    print(f"SUCCESS: {filter}:{frame} computed directly from the star formation history agrees with the from-SED route to {difference.max():.3e} (tolerance {TOLERANCE_ROUTES:.1e})")
    except KeyError as e:
        print(f"FAILED: expected dataset missing from the output: {e}")
        sys.exit(0)

# Confirm that the comparison is sensitive to the interpolation in redshift, by repeating it with a deliberately coarse grid.
print("Running model with a coarse grid of redshifts...")
status, log = runModel("testSuite/parameters/test-luminosities-from-SFH-lightcone.xml testSuite/parameters/test-luminosities-from-SFH-lightcone-coarse.xml","galacticusCoarse.log")
if status != 0 or re.search(r"fatal|aborted",log,re.IGNORECASE) or not os.path.exists(coarsePath):
    print("FAILED: model run with a coarse grid of redshifts:")
    print(log[-4000:])
    sys.exit(0)
with h5py.File(coarsePath,"r") as model:
    try:
        _, ratio         = ratios(model["Lightcone/Output1/nodeData"],FILTER_COARSE,"observed")
        differenceCoarse = np.median(np.abs(ratio-1.0))
    except KeyError as e:
        print(f"FAILED: expected dataset missing from the coarse-grid output: {e}")
        sys.exit(0)
if differenceCoarse < COARSE_MINIMUM or differenceCoarse < FACTOR_COARSE_MINIMUM*differenceFine[FILTER_COARSE]:
    print(f"FAILED: a coarse grid of redshifts is in error by only {differenceCoarse:.3e}, against {differenceFine[FILTER_COARSE]:.3e} for the default grid, so the comparison of the two routes is not demonstrably sensitive to the interpolation in redshift")
else:
    print(f"SUCCESS: a coarse grid of redshifts is in error by {differenceCoarse:.3e}, against {differenceFine[FILTER_COARSE]:.3e} for the default grid, so the comparison of the two routes is sensitive to the interpolation in redshift")

# The direct route requires a star formation history with fixed ages, since its template is shared by every galaxy.
print("Running model with a star formation history which does not use fixed ages...")
status, log = runModel("testSuite/parameters/test-luminosities-from-SFH-lightcone.xml testSuite/parameters/test-luminosities-from-SFH-lightcone-ages.xml","galacticusAges.log")
if status != 0 and re.search(r"must use fixed ages",log):
    print("SUCCESS: a star formation history which does not use fixed ages is reported as an error")
else:
    print(f"FAILED: a star formation history which does not use fixed ages was not reported as an error (exit status {status})")

# A galaxy outside the tabulated range of redshift must be reported, not extrapolated.
print("Running model with galaxies outside the tabulated range of redshift...")
status, log = runModel("testSuite/parameters/test-luminosities-from-SFH-lightcone.xml testSuite/parameters/test-luminosities-from-SFH-lightcone-range.xml","galacticusRange.log")
if status != 0 and re.search(r"lies outside the tabulated range",log):
    print("SUCCESS: a galaxy outside the tabulated range of redshift is reported as an error")
else:
    print(f"FAILED: a galaxy outside the tabulated range of redshift was not reported as an error (exit status {status})")

sys.exit(0)
