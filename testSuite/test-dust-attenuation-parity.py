#!/usr/bin/env python3
"""Check that the dust attenuation framework reproduces the class it replaces.

`nodePropertyExtractorDustAttenuation`, applying `charlotFall2000` to per-component stellar luminosities decomposed
by stellar population age, must reproduce `nodePropertyExtractorLmnstyStllrCF2000` -- the class it is intended to
retire -- when configured equivalently.

Two configuration details make the comparison exact rather than approximate:

* The old class normalizes its power-law extinction curve at a round 5500 A, whereas the framework normalizes at the
  effective wavelength of the Buser V filter, as the stellar spectra classes do. The parameter file sets
  `wavelengthReference` to 5500 A so that this known difference is removed; leaving it at the default shifts
  luminosities by a few parts in ten thousand.
* Rest-frame filters only. The old class evaluates its extinction curve at the observed filter wavelength without
  correcting to the rest frame, so for observed-frame filters the two differ by construction -- the framework being
  the correct one.

Also checked: `outputSum` must equal the sum of the per-component attenuated luminosities exactly, since dust has to
be applied component by component and the total recovered afterwards.

Andrew Benson, with assistance from Claude.
"""
import os
import subprocess
import sys

import h5py
import numpy as np

# Tolerance on the comparison. The two codes are algebraically identical once the reference wavelength is matched,
# but derive the optical-depth-to-magnitudes conversion through different expressions whose values differ in the last
# place, so agreement is to a few ulp rather than exactly zero.
TOLERANCE = 1.0e-12

# This script runs with the working directory set to testSuite/, while the model is run from the repository root, so
# the two need different paths to the same file. Getting this wrong would leave a stale output in place and let a
# failed model run be scored against the previous run's results.
parameterFile = "testSuite/parameters/dustAttenuationParity.xml"
outputFile    = "testSuite/outputs/dustAttenuationParity.hdf5"
outputPath    = os.path.join("..", outputFile)

os.makedirs("outputs", exist_ok=True)
if os.path.exists(outputPath):
    os.remove(outputPath)

status = subprocess.run(f"cd ..; ./Galacticus.exe {parameterFile}", shell=True, capture_output=True)
log = status.stdout.decode() + status.stderr.decode()
if status.returncode != 0 or not os.path.exists(outputPath):
    print("FAILED: the dust attenuation parity model did not run")
    print(log[-2000:])
    sys.exit(0)

with h5py.File(outputPath, "r") as f:
    nodes = f["Outputs/Output1/nodeData"]
    try:
        old      = nodes["luminosityStellar:SDSS_r:rest:dustCharlotFall2000"][:]
        disk     = nodes["luminosityStellar:SDSS_r:rest:disk:dustAttenuated:charlotFall2000"][:]
        spheroid = nodes["luminosityStellar:SDSS_r:rest:spheroid:dustAttenuated:charlotFall2000"][:]
        total    = nodes["luminosityTotal:dustAttenuated:charlotFall2000"][:]
    except KeyError as e:
        print(f"FAILED: expected dataset missing from the output: {e}")
        sys.exit(0)

new       = disk + spheroid
comparable = (old > 0.0) & (new > 0.0)
countGalaxies = int(comparable.sum())

if countGalaxies == 0:
    print("FAILED: no galaxy has a non-zero luminosity in both the old and new extractors")
    sys.exit(0)

# The summed output must equal the sum of its parts exactly -- it is formed from the same numbers.
differenceSum = float(np.nanmax(np.abs(total - new)))
if differenceSum > 0.0:
    print(f"FAILED: outputSum differs from the sum of the per-component luminosities by {differenceSum:.3e}")
    sys.exit(0)
print(f"SUCCESS: outputSum equals the sum of the per-component attenuated luminosities exactly")

difference = np.abs(new[comparable] - old[comparable]) / old[comparable]
worst      = float(difference.max())
if worst > TOLERANCE:
    print(f"FAILED: dust attenuation framework differs from lmnstyStllrCF2000 by up to {worst:.3e}"
          f" (tolerance {TOLERANCE:.1e}) over {countGalaxies} galaxies")
    sys.exit(0)

print(f"SUCCESS: dust attenuation framework reproduces lmnstyStllrCF2000 to {worst:.3e}"
      f" over {countGalaxies} galaxies")

# ---------------------------------------------------------------------------------------------------------------------
# Spectral energy distributions. These check the decomposition machinery itself rather than agreement with a previous
# implementation, since no previous implementation attenuated a spectrum.
# ---------------------------------------------------------------------------------------------------------------------
with h5py.File(outputPath, "r") as f:
    nodes = f["Outputs/Output1/nodeData"]
    try:
        sedRaw      = nodes["diskStellarSED:inoue2014"][:]
        sedZero     = nodes["diskStellarSED:inoue2014:dustAttenuated:zero"][:]
        sedScreen   = nodes["diskStellarSED:inoue2014:dustAttenuated:screenSurfaceDensityMetals"][:]
        sedSequence = nodes["diskStellarSED:inoue2014:dustAttenuated:sequence"][:]
        sedDisk     = nodes["diskStellarSED:inoue2014:dustAttenuated:charlotFall2000"][:]
        sedSpheroid = nodes["spheroidStellarSED:inoue2014:dustAttenuated:charlotFall2000"][:]
        sedTotal    = nodes["sedTotal:dustAttenuated:charlotFall2000"][:]
        wavelengths = nodes["diskStellarSED:inoue2014ColumnValues"][:]
    except KeyError as e:
        print(f"FAILED: expected spectral energy distribution dataset missing from the output: {e}")
        sys.exit(0)

emitting = sedRaw > 0.0
if not emitting.any():
    print("FAILED: every spectral energy distribution is zero, so nothing is being tested")
    sys.exit(0)
print(f"SUCCESS: spectral energy distributions are non-zero at {int(emitting.sum())} wavelength points")

# Decomposing a spectrum and recombining it, with an attenuator that attenuates nothing, must return it unchanged.
# This tests the decomposition independently of any dust physics.
roundTrip = float(np.nanmax(np.abs(sedZero[emitting] - sedRaw[emitting]) / sedRaw[emitting]))
if roundTrip > TOLERANCE:
    print(f"FAILED: the zero attenuator does not reproduce the unattenuated spectrum; differs by {roundTrip:.3e}")
    sys.exit(0)
print(f"SUCCESS: decomposing and recombining a spectrum is lossless to {roundTrip:.3e}")

# The same attenuation reached two ways. `screenSurfaceDensityMetals` is age independent, so the spectrum is
# decomposed with one parcel per wavelength. The sequence contains a birth cloud member, which requests a young/old
# split and so forces a decomposition per star formation history bin, but that member has zero optical depth. Light
# dropped or double counted by the age resolved path would show up as a difference here.
branch = float(np.nanmax(np.abs(sedSequence[emitting] - sedScreen[emitting]) / sedScreen[emitting]))
if branch > TOLERANCE:
    print(f"FAILED: age-resolved and age-independent decompositions disagree by {branch:.3e}")
    sys.exit(0)
print(f"SUCCESS: age-resolved and age-independent decompositions agree to {branch:.3e}")

# Summing over children must be exact, as for the broad-band case.
differenceSED = float(np.nanmax(np.abs(sedTotal - (sedDisk + sedSpheroid))))
if differenceSED > 0.0:
    print(f"FAILED: sedTotal differs from the sum of the per-component spectra by {differenceSED:.3e}")
    sys.exit(0)
print("SUCCESS: sedTotal equals the sum of the per-component attenuated spectra exactly")

# Attenuation can only remove light, and a power-law extinction curve must remove more of it at shorter wavelengths.
if not bool(np.all(sedDisk[emitting] <= sedRaw[emitting] * (1.0 + TOLERANCE))):
    print("FAILED: an attenuated spectrum exceeds the unattenuated one at some wavelength")
    sys.exit(0)
brightest    = int(np.argmax(sedRaw.sum(axis=1)))
emittingHere = sedRaw[brightest] > 0.0
transmission = sedDisk[brightest][emittingHere] / sedRaw[brightest][emittingHere]
if not bool(np.all(np.diff(transmission) >= -TOLERANCE)):
    print("FAILED: transmission is not monotonically increasing with wavelength for a power-law extinction curve")
    sys.exit(0)
print(f"SUCCESS: transmission increases with wavelength, from {transmission[0]:.4f} at"
      f" {wavelengths[emittingHere][0]:.0f} A to {transmission[-1]:.4f} at {wavelengths[emittingHere][-1]:.0f} A")

# ---------------------------------------------------------------------------------------------------------------------
# Emission lines. The extractor sums over disk and spheroid internally, so a correct decomposition has to split by
# component and recombine after attenuation; a multi-line extractor is also a child emitting several elements.
# ---------------------------------------------------------------------------------------------------------------------
lineWavelengths = {"oxygenII3727": 3727.0, "balmerBeta4863": 4863.0, "balmerAlpha6565": 6565.0}

with h5py.File(outputPath, "r") as f:
    nodes = f["Outputs/Output1/nodeData"]
    lineRaw        = {}
    lineAttenuated = {}
    try:
        for line in lineWavelengths:
            lineRaw       [line] = nodes[f"luminosityEmissionLineTotal:{line}"][:]
            lineAttenuated[line] = nodes[f"luminosityEmissionLineTotal:{line}:dustAttenuated:charlotFall2000"][:]
    except KeyError as e:
        print(f"FAILED: expected emission line dataset missing from the output: {e}")
        sys.exit(0)

emittingLines = lineRaw["balmerAlpha6565"] > 0.0
if not emittingLines.any():
    print("FAILED: every emission line luminosity is zero, so nothing is being tested")
    sys.exit(0)
print(f"SUCCESS: emission lines are non-zero for {int(emittingLines.sum())} galaxies")

# Attenuation can only remove light.
for line in lineWavelengths:
    if not bool(np.all(lineAttenuated[line][emittingLines] <= lineRaw[line][emittingLines] * (1.0 + TOLERANCE))):
        print(f"FAILED: the attenuated {line} luminosity exceeds the unattenuated one")
        sys.exit(0)

# A power-law extinction curve must transmit more of a redder line than of a bluer one.
ordered = sorted(lineWavelengths, key=lambda line: lineWavelengths[line])
lineTransmission = {
    line: lineAttenuated[line][emittingLines] / lineRaw[line][emittingLines] for line in ordered
}
for bluer, redder in zip(ordered[:-1], ordered[1:]):
    if not bool(np.all(lineTransmission[bluer] <= lineTransmission[redder] * (1.0 + TOLERANCE))):
        print(f"FAILED: {bluer} is transmitted more than {redder}, but it is the bluer line")
        sys.exit(0)
print(f"SUCCESS: emission line transmission increases with wavelength, from"
      f" {lineTransmission[ordered[0]].mean():.4f} at {lineWavelengths[ordered[0]]:.0f} A to"
      f" {lineTransmission[ordered[-1]].mean():.4f} at {lineWavelengths[ordered[-1]]:.0f} A")

# Reddening must raise the Balmer decrement: Halpha is absorbed less than Hbeta, so their observed ratio exceeds the
# intrinsic one. This is the standard observational signature of dust, and gets the sign of the effect right rather
# than merely its monotonicity.
decrementIntrinsic  = lineRaw       ["balmerAlpha6565"][emittingLines] / lineRaw       ["balmerBeta4863"][emittingLines]
decrementAttenuated = lineAttenuated["balmerAlpha6565"][emittingLines] / lineAttenuated["balmerBeta4863"][emittingLines]
if not bool(np.all(decrementAttenuated > decrementIntrinsic * (1.0 - TOLERANCE))):
    print("FAILED: dust does not redden the Balmer decrement")
    sys.exit(0)
print(f"SUCCESS: dust reddens the Balmer decrement, from {decrementIntrinsic.mean():.4f} to"
      f" {decrementAttenuated.mean():.4f}")

# ---------------------------------------------------------------------------------------------------------------------
# Emitting only the sum. A consumer which expects a single value per galaxy -- an output analysis, for instance --
# cannot use an extractor which emits one property per component, so the wrapper can be asked to emit just the total.
# ---------------------------------------------------------------------------------------------------------------------
with h5py.File(outputPath, "r") as f:
    nodes    = f["Outputs/Output1/nodeData"]
    emitted  = [name for name in nodes.keys() if name.startswith("lineTotalSumOnly")]
    if len(emitted) != 1:
        print(f"FAILED: outputSumOnly emitted {len(emitted)} properties, expected exactly one: {emitted}")
        sys.exit(0)
    sumOnly = nodes[emitted[0]][:]

if sumOnly.ndim != 1:
    print(f"FAILED: outputSumOnly emitted a property of rank {sumOnly.ndim - 1}, expected a scalar")
    sys.exit(0)

# Check the total against the per-line luminosities emitted separately elsewhere in the same run, rather than against
# anything this wrapper produced itself.
expected  = lineAttenuated["balmerAlpha6565"] + lineAttenuated["balmerBeta4863"]
worstSum  = float(np.nanmax(np.abs(sumOnly[emittingLines] - expected[emittingLines]) / expected[emittingLines]))
if worstSum > TOLERANCE:
    print(f"FAILED: the outputSumOnly total differs from the sum of the separately emitted lines by {worstSum:.3e}")
    sys.exit(0)
print(f"SUCCESS: outputSumOnly emits a single scalar equal to the sum of its children, to {worstSum:.3e}")

# ---------------------------------------------------------------------------------------------------------------------
# Scalarizing that sum. This is the chain by which an output analysis reaches the framework: the wrapper emits a single
# summed property, and a `scalarizer` presents it as the scalar the analysis requires.
# ---------------------------------------------------------------------------------------------------------------------
with h5py.File(outputPath, "r") as f:
    nodes      = f["Outputs/Output1/nodeData"]
    emitted    = [name for name in nodes.keys() if name.startswith("lineTotalScalarized")]
    if len(emitted) != 1:
        print(f"FAILED: the scalarizer emitted {len(emitted)} properties, expected exactly one: {emitted}")
        sys.exit(0)
    scalarized = nodes[emitted[0]][:]

if scalarized.ndim != 1:
    print(f"FAILED: the scalarizer emitted a property of rank {scalarized.ndim - 1}, expected a scalar")
    sys.exit(0)

worstScalarized = float(np.nanmax(np.abs(scalarized[emittingLines] - expected[emittingLines]) / expected[emittingLines]))
if worstScalarized > TOLERANCE:
    print(f"FAILED: the scalarized total differs from the sum of the separately emitted lines by {worstScalarized:.3e}")
    sys.exit(0)
print(f"SUCCESS: the scalarizer reproduces the summed, attenuated luminosity, to {worstScalarized:.3e}")

sys.exit(0)
