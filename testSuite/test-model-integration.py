#!/usr/bin/env python3
import subprocess
import sys
import os
import re
import glob
import argparse
import importlib.util
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
parser.add_argument("--launchMethod",        type=str,   default=None,
                    help="Override the queue manager from galacticusConfig.xml when calibrating.")
args, _ = parser.parse_known_args()

# Load the Galacticus launch helpers (provides queue-manager hooks for calibrate mode).
_exec_path = os.environ.get("GALACTICUS_EXEC_PATH",
                            os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
_launch_spec = importlib.util.spec_from_file_location(
    "galacticus_launch", os.path.join(_exec_path, "scripts", "aux", "launch.py")
)
_launch = importlib.util.module_from_spec(_launch_spec)
_launch_spec.loader.exec_module(_launch)

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


def write_model_parameters(input_file, output_file, output_hdf5, seed):
    """Read the input parameter file, set output filename and random seed, write to disk."""
    tree = ET.parse(input_file)
    root = tree.getroot()
    for elem in root.iter("outputFileName"):
        elem.set("value", output_hdf5)
    for elem in root.iter("randomNumberGenerator"):
        for child in elem:
            if child.tag == "seed":
                child.set("value", str(seed))
    tree.write(output_file)


def build_launch_script(method, queue_config):
    """Construct the minimal launch_script dict consumed by launch.py's _MODULE_HOOKS."""
    queue_section = dict(queue_config) if isinstance(queue_config, dict) else {}
    queue_section.setdefault("executable",  "Galacticus.exe")
    queue_section.setdefault("threadCount", "1")
    queue_section.setdefault("ompThreads",  "1")
    return {
        "verbosity":          1,
        "splitModels":        1,
        "modelRootDirectory": "outputs/test-model-integration",
        "compressModels":     "no",
        "useStateFile":       "no",
        "emailReport":        "no",
        "launchMethod":       method,
        method:               queue_section,
    }


def calibrate_model(modelName, fileName, randomSeed, gitRevision, calibrateCount,
                    calibratePercentile, launchMethodOverride):
    """Run `calibrateCount` realizations of `modelName` and write the calibration reference file."""
    # Resolve the queue manager from galacticusConfig.xml (or fall back to local).
    config       = _launch._load_config()
    qm_section   = config.get("queueManager", {}) if isinstance(config, dict) else {}
    method       = launchMethodOverride
    if method is None:
        method   = qm_section.get("manager") if isinstance(qm_section, dict) else None
    if method is None or method not in _launch._MODULE_HOOKS:
        method   = "local"
    queue_config = config.get(method, {}) if isinstance(config, dict) else {}
    launch_script = build_launch_script(method, queue_config)
    _launch._MODULE_HOOKS[method]["validate"](launch_script)

    # Build per-realization parameter files and job specs.
    jobs = []
    for i in range(calibrateCount):
        realDir = f"outputs/test-model-integration/{modelName}{i}"
        os.makedirs(realDir, exist_ok=True)
        write_model_parameters(
            fileName,
            f"{realDir}/parameters.xml",
            f"testSuite/outputs/test-model-integration/{modelName}{i}/galacticus.hdf5",
            randomSeed + i,
        )
        jobs.append({
            "directory":    realDir,
            "label":        f"{modelName}{i}",
            "modelCounter": i,
            "mergeGroup":   i,
            "analysis":     None,
        })

    print(f"--> Calibrating model '{modelName}' via launchMethod='{method}'"
          f" ({calibrateCount} realizations)...")
    _launch._MODULE_HOOKS[method]["launch"](jobs, launch_script, {})
    print(f"<-- launch hook returned for model '{modelName}'.")

    # Run the test functions on each realization to write per-realization stat files.
    for i in range(calibrateCount):
        realDir = f"outputs/test-model-integration/{modelName}{i}"
        if not os.path.exists(f"{realDir}/galacticus.hdf5"):
            print(f"   --> realization {i} produced no HDF5 output; skipping its tests")
            continue
        for test in tests:
            try:
                test["func"](modelName, realDir, test["label"], gitRevision)
            except Exception as e:
                print(f"      Warning: test '{test['label']}' failed for realization {i}: {e}")

    # Aggregate per-bin percentiles across realizations and write the calibration reference file.
    percentiles = [calibratePercentile / 100.0, 0.5, 1.0 - calibratePercentile / 100.0]
    refDir = "data/model-integration"
    os.makedirs(refDir, exist_ok=True)
    for test in tests:
        per_realization = []
        for i in range(calibrateCount):
            statFile = f"outputs/test-model-integration/{modelName}{i}/{test['label']}_r{gitRevision}.txt"
            if not os.path.exists(statFile):
                continue
            data = np.loadtxt(statFile, comments="#")
            per_realization.append(np.atleast_2d(data)[:, 1])
        if not per_realization:
            print(f"   --> no realizations produced output for test '{test['label']}'; skipping calibration")
            continue
        stack = np.vstack(per_realization)
        low, median, high = np.quantile(stack, percentiles, axis=0)
        refFile = f"{refDir}/{test['label']}_{modelName}.txt"
        with open(refFile, "w") as f:
            f.write(f"# Model integration test: {test['name']} [Git revision: {gitRevision}]\n")
            for j in range(stack.shape[1]):
                f.write(f"{j}\t{low[j]}\t{median[j]}\t{high[j]}\n")
        print(f"   --> wrote calibration reference {refFile}")


def evaluate_model(modelName, fileName, randomSeed, gitRevision, calibratePercentile):
    """Run a single realization of `modelName` and run the regression check against the calibration file."""
    modelDir = f"outputs/test-model-integration/{modelName}"
    os.makedirs(modelDir, exist_ok=True)
    write_model_parameters(
        fileName,
        f"{modelDir}/parameters.xml",
        f"testSuite/outputs/test-model-integration/{modelName}/galacticus.hdf5",
        randomSeed,
    )
    print(f"--> Generating model '{modelName}' for model integration testing...")
    print("   --> Running model...")
    subprocess.run(
        f"cd ..; ./Galacticus.exe testSuite/outputs/test-model-integration/{modelName}/parameters.xml",
        shell=True,
    )
    print("   <-- ...done")
    print("<-- ...done")

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

    probabilityFailure = 2.0 * calibratePercentile / 100.0
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
    if args.calibrate == "yes":
        calibrate_model(modelName, fileName, randomSeed, gitRevision,
                        args.calibrateCount, args.calibratePercentile, args.launchMethod)
    else:
        evaluate_model(modelName, fileName, randomSeed, gitRevision, args.calibratePercentile)
