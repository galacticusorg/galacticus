#!/usr/bin/env python3
import subprocess
import sys
import argparse
import h5py
import numpy as np

# Test the radiative transfer code by attempting to reproduce a Strömgren sphere solution.
# Andrew Benson (ported to Python)

PI = np.pi

# Parse command line options.
parser = argparse.ArgumentParser()
parser.add_argument("--processesPerNode", type=int, default=1    )
parser.add_argument("--allowRunAsRoot"  , type=str, default="no" )
parser.add_argument("--oversubscribe"   , type=str, default="yes")
args, _ = parser.parse_known_args()

# We need at least 2 processes to run this test.
if args.processesPerNode < 2:
    print("SKIPPED: at least 2 processes per node are required for this test")
    sys.exit(0)

allowRunAsRoot = " --allow-run-as-root" if args.allowRunAsRoot == "yes" else ""

oversubscribe  = " --oversubscribe"     if args.oversubscribe  == "yes" else ""

# Run the calculation.
status = subprocess.run(
    f"cd ..; mpirun -np {args.processesPerNode}{allowRunAsRoot}{oversubscribe} Galacticus.exe testSuite/parameters/test-radiativeTransfer-StromgrenSphere.xml",
    shell=True
)
if status.returncode != 0:
    print("FAILED: failed to run calculation")
    sys.exit(0)

# Read model output and parameters.
with h5py.File("outputs/radiativeTransfer-StromgrenSphere:MPI0000.hdf5", "r") as outputFile:
    parameters              = outputFile["Parameters"]
    recombination           = parameters["atomicRecombinationRateRadiative"]
    computationalDomain     = parameters["computationalDomain"]
    model                   = outputFile["radiativeTransferModel"]

    densityNumber          = model["densityNumberH"][:]
    fractionHydrogenII     = model["fractionHII"][:]
    recombinationCoefficient  = recombination.attrs["rateCoefficient"]
    xBoundaries            = computationalDomain.attrs["xBoundaries"]
    yBoundaries            = computationalDomain.attrs["yBoundaries"]
    zBoundaries            = computationalDomain.attrs["zBoundaries"]
    countCells             = computationalDomain.attrs["countCells"]
    rateLymanContinuumEmitted = model.attrs["rateLymanContinuumEmitted"]

# Assert that the grid is the same in each direction.
if (countCells[1] != countCells[0] or countCells[2] != countCells[0] or
    yBoundaries[0] != xBoundaries[0] or zBoundaries[0] != xBoundaries[0] or
    yBoundaries[1] != xBoundaries[1] or zBoundaries[1] != xBoundaries[1]):
    print("FAILED: grid is not the same along each axis")
    sys.exit(0)

# Compute computational domain cell volumes.
megaParsec = 3.086e+22
centi      = 1.000e-02
cellVolume = (
    (xBoundaries[1] - xBoundaries[0])
    * (yBoundaries[1] - yBoundaries[0])
    * (zBoundaries[1] - zBoundaries[0])
    / np.prod(countCells)
    * (megaParsec / centi)**3
)
cellSize = xBoundaries[1] - xBoundaries[0]

# Compute the total recombination rate.
recombinationRate = np.sum((densityNumber * fractionHydrogenII)**2) * recombinationCoefficient * cellVolume

# Compute the Strömgren radius.
radiusStromgren = (3.0 * rateLymanContinuumEmitted / 4.0 / PI / densityNumber.max()**2 / recombinationCoefficient)**(1.0/3.0) * centi / megaParsec

# Compute boundary layer volume relative to Strömgren sphere volume.
boundaryVolumeFraction = 3.0 * cellSize / radiusStromgren

# Test for success in the recombination rate.
successRate = "success" if abs(recombinationRate - rateLymanContinuumEmitted) < 2.0 * boundaryVolumeFraction * rateLymanContinuumEmitted else "FAIL"
print(f"{successRate}: ionization balance: recombination / ionization rate = {recombinationRate} / {rateLymanContinuumEmitted} s\u207b\u00b9")

# Check each cell ionization state.
nx, ny, nz      = countCells[0], countCells[1], countCells[2]
insideIsIonized  = True
outsideIsNeutral = True

for i in range(nx):
    xLower   = xBoundaries[0] + (xBoundaries[1] - xBoundaries[0]) / nx * i
    xUpper   = xBoundaries[0] + (xBoundaries[1] - xBoundaries[0]) / nx * (i + 1)
    xMinimum = min(abs(xLower), abs(xUpper))
    xMaximum = max(abs(xLower), abs(xUpper))
    for j in range(ny):
        yLower   = yBoundaries[0] + (yBoundaries[1] - yBoundaries[0]) / ny * j
        yUpper   = yBoundaries[0] + (yBoundaries[1] - yBoundaries[0]) / ny * (j + 1)
        yMinimum = min(abs(yLower), abs(yUpper))
        yMaximum = max(abs(yLower), abs(yUpper))
        for k in range(nz):
            zLower   = zBoundaries[0] + (zBoundaries[1] - zBoundaries[0]) / nz * k
            zUpper   = zBoundaries[0] + (zBoundaries[1] - zBoundaries[0]) / nz * (k + 1)
            zMinimum = min(abs(zLower), abs(zUpper))
            zMaximum = max(abs(zLower), abs(zUpper))
            radiusMinimum = np.sqrt(xMinimum**2 + yMinimum**2 + zMinimum**2)
            radiusMaximum = np.sqrt(xMaximum**2 + yMaximum**2 + zMaximum**2)
            if radiusMaximum < radiusStromgren - cellSize:
                if fractionHydrogenII[i, j, k] <= 0.99:
                    insideIsIonized = False
            elif radiusMinimum > radiusStromgren + cellSize:
                if fractionHydrogenII[i, j, k] >= 0.02:
                    outsideIsNeutral = False

print(f"{'success' if insideIsIonized  else 'FAIL'}: inside Str\u00f6mgren radius is ionized")
print(f"{'success' if outsideIsNeutral else 'FAIL'}: outside Str\u00f6mgren radius is neutral")
