#!/usr/bin/env python3
"""Check AGN tori whose axes are misaligned from those of their galaxies.

Each galaxy is given a random orientation by `nodeOperatorInclinationRandom`, and the axis of its torus a random
misalignment from that of the galaxy by `nodeOperatorMisalignmentRandom`, used by `galacticInclinationMisaligned`. The
model outputs the inclinations of each galaxy and of its torus (whose axis is independent of that of the galaxy), and
the spectrum of the AGN seen through a torus aligned with its galaxy, through one misaligned within a cone of
half-angle zero, and through one whose axis is independent of that of its galaxy. It is checked that:

* inclinations lie in [0°,90°];
* a torus misaligned within a cone of half-angle zero transmits exactly as one aligned with its galaxy;
* the axes of tori independent of their galaxies are isotropic (cos θ uniform, by a Kolmogorov-Smirnov test), and
  uncorrelated with the axes of their galaxies;
* for each galaxy, the AGN is seen through more of the misaligned torus than of the aligned one when the misaligned
  torus is the closer of the two to edge-on, and through less when it is the further.

The statistical checks use fixed seeds, so the outcome is reproducible; their thresholds correspond to a probability
of about 10⁻³ of a correct model failing them.

Andrew Benson, Claude
"""
import os
import subprocess
import sys

import h5py
import numpy as np

# Critical value of the Kolmogorov-Smirnov statistic, times √N, for a significance of 10⁻³.
KS_CRITICAL            = 1.95
# Critical value of the correlation coefficient, times √N, for a two-sided significance of 10⁻³.
CORRELATION_CRITICAL   = 3.29
# Relative tolerance on the comparison of transmissions which should be identical, and the relative difference in
# transmission beyond which a galaxy is counted as seeing its AGN through a different amount of torus.
TOLERANCE_EXACT        = 1.0e-12
DIFFERENCE_SIGNIFICANT = 1.0e-3

# This script runs with the working directory set to testSuite/, while the model is run from the repository root.
parameterFile = "testSuite/parameters/agnTorusMisalignment.xml"
outputFile    = "testSuite/outputs/agnTorusMisalignment.hdf5"
outputPath    = os.path.join("..", outputFile)

os.makedirs("outputs", exist_ok=True)
if os.path.exists(outputPath):
    os.remove(outputPath)

status = subprocess.run(f"cd ..; ./Galacticus.exe {parameterFile}", shell=True, capture_output=True)
log = status.stdout.decode() + status.stderr.decode()
if status.returncode != 0 or not os.path.exists(outputPath):
    print("FAILED: the AGN torus misalignment model did not run")
    print(log[-2000:])
    sys.exit(0)

with h5py.File(outputPath, "r") as f:
    try:
        nodes             = f["Outputs"]["Output1"]["nodeData"]
        inclinationGalaxy = nodes["inclination"      ][:]
        inclinationTorus  = nodes["inclination:torus"][:]
        aligned           = np.sum(nodes["agnSED:dustAttenuated:AGNTorus:aligned"    ][:], axis=1)
        misalignedNone    = np.sum(nodes["agnSED:dustAttenuated:AGNTorus:misaligned0"][:], axis=1)
        misaligned        = np.sum(nodes["agnSED:dustAttenuated:AGNTorus:misaligned" ][:], axis=1)
    except KeyError as e:
        print(f"FAILED: expected dataset missing from the output: {e}")
        sys.exit(0)

countGalaxies = len(inclinationGalaxy)
if countGalaxies < 100:
    print(f"FAILED: only {countGalaxies} galaxies were output, too few for the statistical checks")
    sys.exit(0)

# Inclinations lie in [0°,90°].
for name, inclination in (("galaxies", inclinationGalaxy), ("tori", inclinationTorus)):
    if np.any(inclination < 0.0) or np.any(inclination > 90.0):
        print(f"FAILED: inclinations of {name} lie outside [0°,90°]: {inclination.min()}° to {inclination.max()}°")
        sys.exit(0)
print(f"SUCCESS: inclinations of galaxies and tori lie in [0°,90°], for {countGalaxies} galaxies")

# A torus misaligned within a cone of half-angle zero is aligned with its galaxy.
worst = float(np.max(np.abs(misalignedNone - aligned) / np.maximum(aligned, np.finfo(float).tiny)))
if worst > TOLERANCE_EXACT:
    print(f"FAILED: a torus misaligned within a cone of half-angle zero differs from an aligned torus by up to {worst:.3e}")
    sys.exit(0)
print(f"SUCCESS: a torus misaligned within a cone of half-angle zero transmits as an aligned torus, to {worst:.3e}")

# The axes of independent tori are isotropic, so cos θ is uniform on [0,1].
cosineTorus  = np.sort(np.cos(np.radians(inclinationTorus)))
empirical    = np.arange(1, countGalaxies + 1) / countGalaxies
statisticKS  = float(max(np.max(empirical - cosineTorus), np.max(cosineTorus - (empirical - 1.0 / countGalaxies))))
if statisticKS * np.sqrt(countGalaxies) > KS_CRITICAL:
    print(f"FAILED: cos(θ) of tori independent of their galaxies is not uniform: Kolmogorov-Smirnov statistic {statisticKS:.3f}")
    sys.exit(0)
print(f"SUCCESS: cos(θ) of tori independent of their galaxies is uniform: Kolmogorov-Smirnov statistic {statisticKS:.3f}")

# ...and uncorrelated with the axes of their galaxies.
correlation = float(np.corrcoef(np.cos(np.radians(inclinationGalaxy)), np.cos(np.radians(inclinationTorus)))[0, 1])
if abs(correlation) * np.sqrt(countGalaxies) > CORRELATION_CRITICAL:
    print(f"FAILED: the orientations of tori independent of their galaxies are correlated with them: r={correlation:.3f}")
    sys.exit(0)
print(f"SUCCESS: the orientations of tori independent of their galaxies are uncorrelated with them: r={correlation:.3f}")

# A torus closer to edge-on transmits less. Compare, for each galaxy with an AGN, the misaligned torus with the aligned
# one, which is seen at the inclination of the galaxy.
active     = aligned > 0.0
ratio      = misaligned[active] / aligned[active]
closerFace = inclinationTorus[active] < inclinationGalaxy[active]
violations = int(np.sum(( closerFace & (ratio < 1.0 - TOLERANCE_EXACT)) | (~closerFace & (ratio > 1.0 + TOLERANCE_EXACT))))
different  = int(np.sum(np.abs(ratio - 1.0) > DIFFERENCE_SIGNIFICANT))
if different < 10:
    print(f"FAILED: only {different} galaxies see their AGN through significantly different amounts of torus, so the comparison is not being tested")
    sys.exit(0)
if violations > 0:
    print(f"FAILED: for {violations} of {int(active.sum())} galaxies the misaligned torus transmits more when closer to edge-on, or less when closer to face-on")
    sys.exit(0)
print(f"SUCCESS: misaligned tori transmit more when closer to face-on than their galaxies, and less when closer to edge-on, for {int(active.sum())} galaxies ({different} significantly different)")
sys.exit(0)
