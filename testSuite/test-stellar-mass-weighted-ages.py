#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check calculations of stellar mass-weighted ages.
# Andrew Benson (ported to Python)

# These models used fixed timescales for star formation, and have no inflow or outflow, resulting in analytically solvable
# (exponential) star formation histories. The mean age can therefore be computed directly.

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Simple model - a disk of gas evolves in isolation forming stars.
# Compute expectations.
timescaleDisk = 1.0
timeStart     = 1.0
timeEnd       = 13.8
massDiskStart = 1.0e11
massIntegral     = massDiskStart * (1.0 - np.exp(-(timeEnd - timeStart) / timescaleDisk))
massTimeIntegral = massDiskStart * (timeStart + timescaleDisk - (timeEnd + timescaleDisk) * np.exp(-(timeEnd - timeStart) / timescaleDisk))
ageDiskTarget    = timeEnd - massTimeIntegral / massIntegral

# Run the simple model.
with open("outputs/stellarMassWeightedAgesSimple.log", "w") as logFile:
    status = subprocess.run(
        "cd ..; ./Galacticus.exe testSuite/parameters/stellarMassWeightedAgesSimple.xml",
        shell=True, stdout=logFile, stderr=subprocess.STDOUT
    )
if status.returncode != 0:
    print("FAILED: simple model run:")
    with open("outputs/stellarMassWeightedAgesSimple.log") as f:
        print(f.read())
else:
    print("SUCCESS: simple model run")

# Read the model data and check for consistency.
with h5py.File("outputs/stellarMassWeightedAgesSimple.hdf5", "r") as model:
    ageDisk = model["Outputs/Output1/nodeData/diskAgeStellarMassWeighted"][:]

if abs(ageDisk[0] - ageDiskTarget) < 1.0e-3:
    print("SUCCESS: simple model age")
else:
    print(f"FAILED: simple model age: {ageDisk[0]} \u2247 {ageDiskTarget}")


# Merging model - two disks of gas evolve in isolation up to t=6 Gyr, then merge, forming a spheroid.
# Compute expectations.
timescaleDisk     = 1.00
timescaleSpheroid = 0.75
timeStartCentral  = 1.00
timeStartSatellite= 3.00
timeMerge         = 6.00
timeEnd           = 13.80
massDiskStart     = 1.00e11

massIntegralCentralPreMerge       = massDiskStart * (1.0 - np.exp(-(timeMerge - timeStartCentral) / timescaleDisk))
massTimeIntegralCentralPreMerge   = massDiskStart * (timeStartCentral + timescaleDisk - (timeMerge + timescaleDisk) * np.exp(-(timeMerge - timeStartCentral) / timescaleDisk))
massGasCentralMerge               = massDiskStart * np.exp(-(timeMerge - timeStartCentral) / timescaleDisk)
massIntegralSatellitePreMerge     = massDiskStart * (1.0 - np.exp(-(timeMerge - timeStartSatellite) / timescaleDisk))
massTimeIntegralSatellitePreMerge = massDiskStart * (timeStartSatellite + timescaleDisk - (timeMerge + timescaleDisk) * np.exp(-(timeMerge - timeStartSatellite) / timescaleDisk))
massGasSatelliteMerge             = massDiskStart * np.exp(-(timeMerge - timeStartSatellite) / timescaleDisk)
massSpheroidStart                 = massGasCentralMerge + massGasSatelliteMerge
massIntegralCentralPostMerge      = massSpheroidStart * (1.0 - np.exp(-(timeEnd - timeMerge) / timescaleSpheroid))
massTimeIntegralCentralPostMerge  = massSpheroidStart * (timeMerge + timescaleSpheroid - (timeEnd + timescaleSpheroid) * np.exp(-(timeEnd - timeMerge) / timescaleSpheroid))
massIntegral                      = massIntegralCentralPreMerge + massIntegralSatellitePreMerge + massIntegralCentralPostMerge
massTimeIntegral                  = massTimeIntegralCentralPreMerge + massTimeIntegralSatellitePreMerge + massTimeIntegralCentralPostMerge
ageSpheroidTarget                 = timeEnd - massTimeIntegral / massIntegral

# Run the merging model.
with open("outputs/stellarMassWeightedAgesMerging.log", "w") as logFile:
    status = subprocess.run(
        "cd ..; ./Galacticus.exe testSuite/parameters/stellarMassWeightedAgesMerging.xml",
        shell=True, stdout=logFile, stderr=subprocess.STDOUT
    )
if status.returncode != 0:
    print("FAILED: merging model run:")
    with open("outputs/stellarMassWeightedAgesMerging.log") as f:
        print(f.read())
else:
    print("SUCCESS: merging model run")

# Read the model data and check for consistency.
with h5py.File("outputs/stellarMassWeightedAgesMerging.hdf5", "r") as model:
    ageSpheroid = model["Outputs/Output1/nodeData/spheroidAgeStellarMassWeighted"][:]

if abs(ageSpheroid[0] - ageSpheroidTarget) < 1.0e-3:
    print("SUCCESS: merging model age")
else:
    print(f"FAILED: merging model age: {ageSpheroid[0]} \u2247 {ageSpheroidTarget}")
