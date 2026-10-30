#!/usr/bin/env python3
"""Check that black hole mass growth accounts for the energy carried off by jets.

A black hole accreting rest mass at a rate Ṁ grows at a rate Ṁ(1-ε_rad) - P_jet/c², where ε_rad is the radiative
efficiency of its accretion disk and P_jet its jet power (`nodeOperatorBlackHolesAccretion`). The model evolves a single isolated halo
(so that the black hole never merges) whose spheroid gas feeds, by Bondi-Hoyle accretion, a rapidly spinning black hole
through an ADAF, which launches a strong jet. Star formation and feedback are switched off. The accretion rate, radiative
efficiency and jet power are output at closely spaced times, and the growth of the black hole mass is checked against
the integral of the growth rate above. Omitting the jet term changes that integral by about 20%, so the check would
catch the jet power being dropped, or mis-scaled (as it was, by a factor 10¹², in the conversion from jet power to
efficiency fixed in October 2026).

Andrew Benson, Claude
"""
import os
import subprocess
import sys

import h5py
import numpy as np

# Speed of light (km/s): jet power is in units of M☉ (km/s)² Gyr⁻¹, and accretion rates are in M☉ Gyr⁻¹.
SPEED_LIGHT              = 2.99792458e5
# Relative tolerance on the integrated black hole mass growth. The trapezoidal integration over the output times is
# accurate to about 3×10⁻⁴, while omitting the jet term gives a difference of about 0.2.
TOLERANCE_GROWTH         = 1.0e-2
# Minimum fractional growth of the black hole, and minimum mean jet efficiency, for the check to be sensitive.
GROWTH_FRACTION_MINIMUM  = 1.0
EFFICIENCY_JET_MINIMUM   = 0.05

# This script runs with the working directory set to testSuite/, while the model is run from the repository root.
parameterFile = "testSuite/parameters/blackHoleJetEfficiency.xml"
outputFile    = "testSuite/outputs/blackHoleJetEfficiency.hdf5"
outputPath    = os.path.join("..", outputFile)

os.makedirs("outputs", exist_ok=True)
if os.path.exists(outputPath):
    os.remove(outputPath)

status = subprocess.run(f"cd ..; ./Galacticus.exe {parameterFile}", shell=True, capture_output=True)
log = status.stdout.decode() + status.stderr.decode()
if status.returncode != 0 or not os.path.exists(outputPath):
    print("FAILED: the black hole jet efficiency model did not run")
    print(log[-2000:])
    sys.exit(0)
print("SUCCESS: model run")

# Read the black hole properties at each output, in order of time.
time          = []
massBlackHole = []
rateAccretion = []
powerJet      = []
efficiencyRad = []
with h5py.File(outputPath, "r") as f:
    try:
        for name in sorted(f["Outputs"].keys(), key=lambda s: int(s.replace("Output", ""))):
            nodes = f["Outputs"][name]["nodeData"]
            if nodes["nodeIndex"].shape[0] != 1:
                print(f"FAILED: expected one node in {name}, found {nodes['nodeIndex'].shape[0]}")
                sys.exit(0)
            # Black hole properties are output as one list per node, with one entry per black hole.
            rates = np.atleast_1d(nodes["massAccretionRateBlackHoles"  ][0])
            if rates.size != 1:
                print(f"FAILED: expected one black hole in {name}, found {rates.size}")
                sys.exit(0)
            time         .append(nodes["time"         ][0])
            massBlackHole.append(nodes["blackHoleMass"][0])
            rateAccretion.append(rates[0])
            powerJet     .append(np.atleast_1d(nodes["powerJetBlackHoles"           ][0])[0])
            efficiencyRad.append(np.atleast_1d(nodes["radiativeEfficiencyBlackHoles"][0])[0])
    except KeyError as e:
        print(f"FAILED: expected dataset missing from the output: {e}")
        sys.exit(0)
time, massBlackHole, rateAccretion, powerJet, efficiencyRad = map(np.array, (time, massBlackHole, rateAccretion, powerJet, efficiencyRad))

# At the first output (the initial time of the tree) the structure of the galaxy has not yet been solved for, so no
# accretion rate is reported. Begin the comparison at the first output with accretion.
accreting = np.nonzero(rateAccretion > 0.0)[0]
if accreting.size < 10:
    print(f"FAILED: the black hole is accreting at only {accreting.size} of {time.size} outputs")
    sys.exit(0)
start = accreting[0]
if np.any(rateAccretion[start:] <= 0.0):
    print("FAILED: black hole accretion stops part way through the model")
    sys.exit(0)
time, massBlackHole, rateAccretion, powerJet, efficiencyRad = (x[start:] for x in (time, massBlackHole, rateAccretion, powerJet, efficiencyRad))

# The check is sensitive to the jet term only if the black hole grows substantially and the jet is strong.
efficiencyJet = powerJet / rateAccretion / SPEED_LIGHT**2
growth        = massBlackHole - massBlackHole[0]
if growth[-1] < GROWTH_FRACTION_MINIMUM * massBlackHole[0]:
    print(f"FAILED: the black hole grew by only a factor {massBlackHole[-1]/massBlackHole[0]:.3f}, too little for the check to be sensitive")
    sys.exit(0)
if np.mean(efficiencyJet) < EFFICIENCY_JET_MINIMUM:
    print(f"FAILED: the mean jet efficiency is only {np.mean(efficiencyJet):.3e}, too small for the check to be sensitive")
    sys.exit(0)
print(f"SUCCESS: the black hole grows by a factor {massBlackHole[-1]/massBlackHole[0]:.2f}, with mean jet efficiency {np.mean(efficiencyJet):.3f}")

# Integrate the growth rate over the output times (trapezoidal rule), and compare with the growth of the black hole mass.
def integrate(rate):
    return np.concatenate(([0.0], np.cumsum(0.5 * (rate[1:] + rate[:-1]) * np.diff(time))))
growthExpected  = integrate(rateAccretion * (1.0 - efficiencyRad) - powerJet / SPEED_LIGHT**2)
growthNoJet     = integrate(rateAccretion * (1.0 - efficiencyRad)                            )
errorRelative   = float(np.max(np.abs(growth - growthExpected)) / growth[-1])
differenceNoJet = float(abs(growth[-1] - growthNoJet[-1]) / growth[-1])
if errorRelative > TOLERANCE_GROWTH:
    print(f"FAILED: black hole mass growth differs from the accreted mass less radiated and jet energy by up to {errorRelative:.3e} (relative to the total growth); with the jet energy omitted the difference would be {differenceNoJet:.3e}")
else:
    print(f"SUCCESS: black hole mass growth equals the accreted mass less radiated and jet energy, to {errorRelative:.3e} (omitting the jet energy would give {differenceNoJet:.3e})")
sys.exit(0)
