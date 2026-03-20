#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check calculations of mean interoutput star formation rates.
# Andrew Benson (ported to Python)

# Run the model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/interoutputStarFormationRate.xml", shell=True)
if status.returncode != 0:
    print("FAILED: model run")
else:
    print("SUCCESS: model run")

# Read the model data and check for consistency.
with h5py.File("outputs/interoutputStarFormationRate.hdf5", "r") as model:
    recycledFraction = model["Parameters/stellarPopulation"].attrs["recycledFraction"]
    outputs          = model["Outputs"]
    timePrevious     = 0.0
    massPrevious     = np.zeros(6)
    allGood          = True
    for i in range(1, 4):
        output                       = outputs[f"Output{i}"]
        nodeData                     = output["nodeData"]
        time                         = output.attrs["outputTime"]
        rateStarFormationInterOutput = nodeData["diskStarFormationRateInterOutputMean"][:]
        massStellar                  = nodeData["diskMassStellar"][:]
        treeIndex                    = nodeData["mergerTreeIndex"][:]
        order                        = np.argsort(treeIndex)
        # Compute increase in mass directly and from star formation rate.
        massIncrease         = massStellar[order] - massPrevious
        timeIncrease         = time - timePrevious
        massIncreaseFromRate = rateStarFormationInterOutput[order] * timeIncrease * (1.0 - recycledFraction)
        if np.any((np.abs(massIncreaseFromRate - massIncrease) > 1.0e-2 * massIncrease) & (massIncrease > 1.0)):
            allGood = False
        timePrevious = time
        massPrevious = massStellar[order]

success = "success" if allGood else "FAIL"
print(f"{success}: mean interoutput star formation rate")
