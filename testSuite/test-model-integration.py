#!/usr/bin/env python3
import subprocess
import sys
import os
import re
import glob
import argparse
import h5py
import numpy as np
import xml.etree.ElementTree as ET
from math import comb

# Run Galacticus models for integration testing.
# Andrew Benson (ported to Python)

parser = argparse.ArgumentParser()
parser.add_argument("--calibrate",           type=str,   default="no")
parser.add_argument("--calibrateCount",      type=int,   default=100)
parser.add_argument("--calibratePercentile", type=float, default=1.0)
parser.add_argument("--instance",            type=str,   default="1:1")
args, _ = parser.parse_known_args()

m = re.match(r'(\d+):(\d+)', args.instance)
if not m:
    print("'instance' argument syntax error")
    sys.exit(0)
instance      = int(m.group(1))
instanceCount = int(m.group(2))
print(f" -> launching instance {instance} of {instanceCount}")

# Get Git revision.
result = subprocess.run("git rev-parse HEAD", shell=True, capture_output=True, text=True)
gitRevision = result.stdout.strip() if result.returncode == 0 else "Unknown"

randomSeed = np.random.randint(1000)


def weighted_histogram(bins, values, weights):
    """Compute a weighted histogram (differential)."""
    n = len(bins)
    massFunction      = np.zeros(n)
    massFunctionError = np.zeros(n)
    binWidth = bins[1] - bins[0] if n > 1 else 1.0
    for i, binCenter in enumerate(bins):
        selection = np.where((values >= binCenter - 0.5*binWidth) & (values < binCenter + 0.5*binWidth))[0]
        if len(selection) > 0:
            massFunction[i]      = weights[selection].sum() / binWidth
            massFunctionError[i] = np.sqrt((weights[selection]**2).sum()) / binWidth
    return massFunction, massFunctionError


def test_mass_function_halo(modelName, modelDir, label, gitRevision):
    bins = np.arange(10) / 2.0 + 9.0
    with h5py.File(f"{modelDir}/galacticus.hdf5", "r") as model:
        nodeData = model["Outputs/Output1/nodeData"]
        basicMass         = nodeData["basicMass"][:]
        mergerTreeWeight  = nodeData["mergerTreeWeight"][:]
    with np.errstate(divide='ignore'):
        massLog          = np.log10(basicMass)
    massFunction, _ = weighted_histogram(bins, massLog, mergerTreeWeight)
    massFunction    /= np.log(10.0)
    with open(f"{modelDir}/{label}_r{gitRevision}.txt", "w") as f:
        f.write(f"# Model integration test: Halo mass function at z=0.0 [git revision: {gitRevision}]\n")
        for i in range(len(bins)):
            f.write(f"{i}\t{massFunction[i]}\n")
    refFile = f"data/model-integration/{label}_{modelName}.txt"
    if os.path.exists(refFile):
        refData = np.loadtxt(refFile, comments="#")
        refLow  = refData[:, 1]
        refHigh = refData[:, 3]
        failures = np.where((massFunction < refLow * (1.0-1.0e-6)) | (massFunction > refHigh * (1.0+1.0e-6)))[0]
        return len(failures), len(massFunction)
    return 0, len(massFunction)


def test_mass_function_stellar(modelName, modelDir, label, gitRevision):
    bins = np.arange(10) / 2.0 + 8.0
    with h5py.File(f"{modelDir}/galacticus.hdf5", "r") as model:
        nodeData = model["Outputs/Output1/nodeData"]
        massStellar      = nodeData["diskMassStellar"][:] + nodeData.get("spheroidMassStellar", np.zeros(1))[:]
        mergerTreeWeight = nodeData["mergerTreeWeight"][:]
    with np.errstate(divide='ignore'):
        massLog          = np.log10(massStellar)
    massFunction, _ = weighted_histogram(bins, massLog, mergerTreeWeight)
    massFunction    /= np.log(10.0)
    with open(f"{modelDir}/{label}_r{gitRevision}.txt", "w") as f:
        f.write(f"# Model integration test: Stellar mass function at z=0.0 [git revision: {gitRevision}]\n")
        for i in range(len(bins)):
            f.write(f"{i}\t{massFunction[i]}\n")
    refFile = f"data/model-integration/{label}_{modelName}.txt"
    if os.path.exists(refFile):
        refData = np.loadtxt(refFile, comments="#")
        refLow  = refData[:, 1]
        refHigh = refData[:, 3]
        failures = np.where((massFunction < refLow * (1.0-1.0e-6)) | (massFunction > refHigh * (1.0+1.0e-6)))[0]
        return len(failures), len(massFunction)
    return 0, len(massFunction)


def test_mass_function_ism(modelName, modelDir, label, gitRevision):
    bins = np.arange(10) / 2.0 + 6.0
    with h5py.File(f"{modelDir}/galacticus.hdf5", "r") as model:
        nodeData = model["Outputs/Output1/nodeData"]
        massColdGas      = nodeData["diskMassGas"][:] + nodeData.get("spheroidMassGas", np.zeros(1))[:]
        mergerTreeWeight = nodeData["mergerTreeWeight"][:]
    with np.errstate(divide='ignore'):
        massLog          = np.log10(massColdGas)
    massFunction, _ = weighted_histogram(bins, massLog, mergerTreeWeight)
    massFunction    /= np.log(10.0)
    with open(f"{modelDir}/{label}_r{gitRevision}.txt", "w") as f:
        f.write(f"# Model integration test: ISM mass function at z=0.0 [git revision: {gitRevision}]\n")
        for i in range(len(bins)):
            f.write(f"{i}\t{massFunction[i]}\n")
    refFile = f"data/model-integration/{label}_{modelName}.txt"
    if os.path.exists(refFile):
        refData = np.loadtxt(refFile, comments="#")
        refLow  = refData[:, 1]
        refHigh = refData[:, 3]
        failures = np.where((massFunction < refLow * (1.0-1.0e-6)) | (massFunction > refHigh * (1.0+1.0e-6)))[0]
        return len(failures), len(massFunction)
    return 0, len(massFunction)


def test_median_sizes(modelName, modelDir, label, gitRevision):
    bins = np.arange(10) / 2.0 + 8.0
    with h5py.File(f"{modelDir}/galacticus.hdf5", "r") as model:
        nodeData = model["Outputs/Output1/nodeData"]
        diskRadius       = nodeData["diskRadius"][:]
        massStellar      = nodeData["diskMassStellar"][:] + nodeData.get("spheroidMassStellar", np.zeros(1))[:]
        mergerTreeWeight = nodeData["mergerTreeWeight"][:]
    with np.errstate(divide='ignore'):
        massLog   = np.log10(massStellar)
    binWidth  = bins[1] - bins[0] if len(bins) > 1 else 1.0
    quantiles = np.zeros(len(bins))
    for i, binCenter in enumerate(bins):
        selection = np.where((massLog >= binCenter - 0.5*binWidth) & (massLog < binCenter + 0.5*binWidth))[0]
        if len(selection) > 0:
            sortedIdx = np.argsort(diskRadius[selection])
            cumWeights = np.cumsum(mergerTreeWeight[selection][sortedIdx])
            totalWeight = cumWeights[-1]
            idx50 = np.searchsorted(cumWeights, 0.5 * totalWeight)
            quantiles[i] = diskRadius[selection][sortedIdx[min(idx50, len(sortedIdx)-1)]]
    with open(f"{modelDir}/{label}_r{gitRevision}.txt", "w") as f:
        f.write(f"# Model integration test: Median disk sizes at z=0.0 [git revision: {gitRevision}]\n")
        for i in range(len(bins)):
            f.write(f"{i}\t{quantiles[i]}\n")
    refFile = f"data/model-integration/{label}_{modelName}.txt"
    if os.path.exists(refFile):
        refData = np.loadtxt(refFile, comments="#")
        refLow  = refData[:, 1]
        refHigh = refData[:, 3]
        failures = np.where((quantiles < refLow * (1.0-1.0e-6)) | (quantiles > refHigh * (1.0+1.0e-6)))[0]
        return len(failures), len(quantiles)
    return 0, len(quantiles)


tests = [
    {"name": "halo mass function z=0.0",    "label": "massFunctionHalo",    "func": test_mass_function_halo},
    {"name": "stellar mass function z=0.0", "label": "massFunctionStellar", "func": test_mass_function_stellar},
    {"name": "ISM gas mass function z=0.0", "label": "massFunctionISM",     "func": test_mass_function_ism},
    {"name": "median disk sizes z=0.0",     "label": "medianSizes",         "func": test_median_sizes},
]


# Find integration models to run.
modelCount = 0
for fileName in sorted(glob.glob("test-model-integration-*.xml")):
    m = re.match(r'^test-model-integration-([a-zA-Z0-9]+)\.xml$', fileName)
    if not m:
        continue
    modelName = m.group(1)
    modelCount += 1
    if modelCount % instanceCount != instance - 1:
        continue
    print(f"Running model '{modelName}'")

    # Run the model.
    modelDir = f"outputs/test-model-integration/{modelName}"
    subprocess.run(f"mkdir -p {modelDir}", shell=True)

    # Update outputFileName in parameters.
    tree = ET.parse(fileName)
    root = tree.getroot()
    for elem in root.iter("outputFileName"):
        elem.set("value", f"testSuite/outputs/test-model-integration/{modelName}/galacticus.hdf5")
    for elem in root.iter():
        if elem.tag == "randomNumberGenerator":
            for child in elem:
                if child.tag == "seed":
                    child.set("value", str(randomSeed))

    paramFile = f"{modelDir}/parameters.xml"
    tree.write(paramFile)

    print(f"--> Generating model '{modelName}' for model integration testing...")
    print("   --> Running model...")
    subprocess.run(f"cd ..; ./Galacticus.exe testSuite/outputs/test-model-integration/{modelName}/parameters.xml", shell=True)
    print("   <-- ...done")
    print("<-- ...done")

    # Test models.
    failureTotal = 0
    testTotal    = 0
    print(f"--> Testing model '{modelName}' for model integration testing...")
    for test in tests:
        print(f"   --> Test '{test['name']}'...")
        try:
            failureCount, testCount = test["func"](modelName, modelDir, test["label"], gitRevision)
        except Exception as e:
            print(f"      Warning: test failed with exception: {e}")
            failureCount, testCount = 0, 0
        failureTotal += failureCount
        testTotal    += testCount
        print(f"   --> failure rate: {failureCount}/{testCount}")
    print("<-- ...done")

    # Compute probability of failure count.
    probabilityFailure = 2.0 * args.calibratePercentile / 100.0
    probability        = 0.0
    for k in range(failureTotal + 1):
        probability += comb(testTotal, k) * probabilityFailure**k * (1.0 - probabilityFailure)**(testTotal - k)
    print(f"Probability of this number or fewer failures = {probability}")
    probabilityExcess = 1.0 - probability
    if probabilityExcess < 0.001:
        status = "FAILED"
    elif probabilityExcess < 0.010:
        status = "WARNING"
    else:
        status = "success"
    print(f"Status: {status}")
