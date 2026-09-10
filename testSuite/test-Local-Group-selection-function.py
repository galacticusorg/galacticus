#!/usr/bin/env python3
"""Test the observational selection function applied to Local Group satellite galaxies.

The `localGroupDetection` weight operator weights each model satellite by the probability that it
would be detected in the DELVE Milky Way Census I, interpolating a tabulation of that probability as
a function of galactocentric radius, absolute V-band magnitude, and projected half-light radius.

This test runs a small model which writes both the `localGroupLuminosityFunction` analysis and a
catalog of the node properties which that weight operator uses, then recomputes the weighted
satellite counts directly from the same tabulation and checks that the two agree. It therefore
exercises the whole path -- the units of the tabulated axes, the properties extracted from each
node, and the trilinear interpolation -- rather than any one piece of it.

The random error convolution of the analysis is disabled in the parameter file so that the model
distribution is a plain weighted histogram, directly comparable with the recomputation.

Andrew Benson
"""

import os
import subprocess
import sys

import h5py
import numpy as np

# Detection probabilities are interpolated, so agreement should be at round-off level. This is set
# well above that, but far below any plausible error in the interpolation.
toleranceFractional = 1.0e-6

# The ratio of projected half-light radius to three-dimensional stellar half-mass radius. Must match
# the default of the `localGroupDetection` weight operator.
radiusHalfLightRatio = 0.75

# The outer radius which defines the sample. Must match that used by the analysis.
radiusOuter = 0.3

fileNameModel = "outputs/test-Local-Group-selection-function.hdf5"
fileNameTable = os.path.join(
    os.environ.get("GALACTICUS_DATA_PATH", ""),
    "static/observations/localGroup/localGroupSelectionFunctionDELVE.hdf5",
)


def interpolationFactors(values, axis):
    """Return the lower index and upper weight for linear interpolation onto a monotonic axis."""
    indexUpper = np.clip(np.searchsorted(axis, values, side="left"), 1, axis.size - 1)
    weightUpper = np.clip((values - axis[indexUpper - 1]) / (axis[indexUpper] - axis[indexUpper - 1]), 0.0, 1.0)
    return indexUpper - 1, weightUpper


def detectionProbability(table, radius, magnitude, radiusHalfLight):
    """Recompute the detection probability, following the weight operator."""
    axisRadius = table["radiusGalactocentric"][:]
    axisMagnitude = table["magnitudeAbsoluteV"][:]
    axisRadiusHalf = table["radiusHalfLight"][:]
    # Transposed on read: the file is written such that Galacticus sees
    # (radiusGalactocentric,magnitudeAbsoluteV,radiusHalfLight).
    probabilityTable = table["detectionProbability"][:].T
    inRange = (radius >= axisRadius[0]) & (radius <= axisRadius[-1]) & np.isfinite(magnitude)
    magnitude = np.clip(magnitude, axisMagnitude[0], axisMagnitude[-1])
    radiusHalfLight = np.clip(radiusHalfLight, axisRadiusHalf[0], axisRadiusHalf[-1])
    indexRadius, weightRadius = interpolationFactors(
        np.log10(np.maximum(radius, axisRadius[0])), np.log10(axisRadius)
    )
    indexMagnitude, weightMagnitude = interpolationFactors(magnitude, axisMagnitude)
    indexRadiusHalf, weightRadiusHalf = interpolationFactors(np.log10(radiusHalfLight), np.log10(axisRadiusHalf))
    probability = np.zeros(radius.size)
    for offsetRadius, factorRadius in ((0, 1.0 - weightRadius), (1, weightRadius)):
        for offsetMagnitude, factorMagnitude in ((0, 1.0 - weightMagnitude), (1, weightMagnitude)):
            for offsetRadiusHalf, factorRadiusHalf in ((0, 1.0 - weightRadiusHalf), (1, weightRadiusHalf)):
                probability += (
                    factorRadius
                    * factorMagnitude
                    * factorRadiusHalf
                    * probabilityTable[
                        indexRadius + offsetRadius,
                        indexMagnitude + offsetMagnitude,
                        indexRadiusHalf + offsetRadiusHalf,
                    ]
                )
    return np.where(inRange, np.clip(probability, 0.0, 1.0), 0.0)


# Create output path.
try:
    os.mkdir("outputs")
except FileExistsError:
    pass

# Run the model.
status = subprocess.run(
    "cd ..; ./Galacticus.exe testSuite/parameters/test-Local-Group-selection-function.xml", shell=True
)
if status.returncode != 0:
    print("FAILED: Local Group selection function model failed to run")
    sys.exit(0)
print("SUCCESS: Local Group selection function model run")

if not os.path.exists(fileNameTable):
    print(f"FAILED: tabulated selection function '{fileNameTable}' not found")
    sys.exit(0)

with h5py.File(fileNameModel, "r") as model, h5py.File(fileNameTable, "r") as table:
    nodes = model["Outputs/Output1/nodeData"]
    isIsolated = nodes["nodeIsIsolated"][:]
    radius = nodes["radiusOrbitalCurrent"][:]
    radiusHalfMass = nodes["radiusHalfMassStellar"][:]
    weightTree = nodes["mergerTreeWeight"][:]
    # The total stellar luminosity -- not the per-component luminosities, which are also present.
    luminosity = nodes["luminosityStellar:Buser_V:rest"][:]
    analysis = model["analyses/localGroupLuminosityFunction"]
    magnitudeBin = analysis["magnitudeAbsoluteV"][:]
    luminosityFunction = analysis["luminosityFunction"][:]

    with np.errstate(divide="ignore"):
        magnitude = np.where(
            luminosity > 0.0, -2.5 * np.log10(np.maximum(luminosity, np.finfo(float).tiny)), np.inf
        )
    probability = detectionProbability(table, radius, magnitude, radiusHalfLightRatio * radiusHalfMass)
    probability = np.where(radius <= radiusOuter, probability, 0.0)

    isSatellite = isIsolated == 0
    isCentral = isIsolated == 1
    if not np.any(probability[isSatellite] > 0.0):
        print("FAILED: no model satellite has a non-zero detection probability - the test is not meaningful")
        sys.exit(0)

    # Bin the satellites as the analysis does, normalizing to the number of host galaxies.
    binWidth = magnitudeBin[1] - magnitudeBin[0]
    edges = np.concatenate([magnitudeBin - 0.5 * binWidth, [magnitudeBin[-1] + 0.5 * binWidth]])
    counts, _ = np.histogram(
        magnitude[isSatellite], bins=edges, weights=(weightTree * probability)[isSatellite]
    )
    counts /= np.sum(weightTree[isCentral])

    scale = max(np.max(luminosityFunction), np.max(counts))
    deviation = np.max(np.abs(luminosityFunction - counts)) / scale
    if deviation < toleranceFractional:
        print(f"SUCCESS: Local Group selection function weighting (maximum deviation {deviation:.3e})")
    else:
        print(f"FAILED: Local Group selection function weighting (maximum deviation {deviation:.3e})")
        for magnitude_, fortran, python in zip(magnitudeBin, luminosityFunction, counts):
            print(f"   M_V = {magnitude_:6.2f}: analysis = {fortran:12.6g}, recomputed = {python:12.6g}")
