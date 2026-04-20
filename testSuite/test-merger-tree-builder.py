#!/usr/bin/env python3
import subprocess
import sys
import os
import glob
import argparse
import shutil
import h5py
import numpy as np

# Compare conditional mass functions with a reference set.
# Andrew Benson (ported to Python)

# Parse command line options.
parser = argparse.ArgumentParser()
parser.add_argument("--launchMethod",  type=str, default="local")
parser.add_argument("--threadMaximum", type=int, default=1)
parser.add_argument("--ompThreads",    type=int, default=4)
parser.add_argument("--instance",      type=str, default=None)
args, _ = parser.parse_known_args()

# Specify models.
models = [
    {
        "type":      "reference",
        "fileName":  "data/mergerTreeBuilding/test-merger-tree-builder-reference.hdf5",
    },
    {
        "type":      "referenceGeneric",
        "fileName":  "data/mergerTreeBuilding/test-merger-tree-builder-reference-generic.hdf5",
    },
    {
        "type":            "Cole et al. (2000); intervalStep=true",
        "fileName":        "outputs/mergerTreeBuilderCole2000_intervalStepTrue.hdf5",
        "parameterFile":   "parameters/mergerTreeBuilderCole2000_intervalStepTrue.xml",
        "reference":       "reference",
        "toleranceSigma":       6.0,
        "toleranceFractional":  1.0e-3,
    },
    {
        "type":            "Cole et al. (2000); intervalStep=false",
        "fileName":        "outputs/mergerTreeBuilderCole2000_intervalStepFalse.hdf5",
        "parameterFile":   "parameters/mergerTreeBuilderCole2000_intervalStepFalse.xml",
        "reference":       "reference",
        "toleranceSigma":       9.0,
        "toleranceFractional":  1.0e-3,
    },
    {
        "type":            "genericLinearBarrier",
        "fileName":        "outputs/mergerTreeBuilderGenericLinearBarrier.hdf5",
        "parameterFile":   "parameters/mergerTreeBuilderGenericLinearBarrier.xml",
        "reference":       "referenceGeneric",
        "toleranceSigma":       6.0,
        "toleranceFractional":  1.0e-3,
    },
    {
        "type":            "genericSolver",
        "fileName":        "outputs/mergerTreeBuilderGenericSolver.hdf5",
        "parameterFile":   "parameters/mergerTreeBuilderGenericSolver.xml",
        "reference":       "referenceGeneric",
        "toleranceSigma":       9.0,
        "toleranceFractional":  1.0e-3,
    },
]

# Build the parameterGrid document.
subprocess.run("mkdir -p outputs/test-merger-tree-builder", shell=True)

paramFiles = []
for model in models:
    if "parameterFile" in model:
        paramFiles.append(model["parameterFile"])

# Run each model directly (in local mode).
launchOptions = f"--launchMethod {args.launchMethod} --threadMaximum {args.threadMaximum} --ompThreads {args.ompThreads}"
if args.instance:
    launchOptions += f" --instance {args.instance}"

# Build parameterGrid XML and run via launch.py.
with open("outputs/test-merger-tree-builder.xml", "w") as f:
    f.write("<parameterGrid>\n")
    f.write("  <emailReport>no</emailReport>\n")
    f.write("  <doAnalysis>no</doAnalysis>\n")
    f.write("  <modelRootDirectory>testSuite/outputs/test-merger-tree-builder</modelRootDirectory>\n")
    f.write("  <baseParameters>testSuite/parameters/mergerTreeBuilderCole2000_intervalStepTrue.xml</baseParameters>\n")
    for paramFile in paramFiles:
        with open(paramFile) as pf:
            f.write(pf.read())
    f.write("</parameterGrid>\n")

subprocess.run(
    f"cd ..; mkdir -p testSuite/outputs/test-merger-tree-builder; ./scripts/aux/launch.py testSuite/outputs/test-merger-tree-builder.xml {launchOptions}",
    shell=True
)

# Check for failed models.
logFiles = glob.glob("outputs/test-merger-tree-builder/galacticus_*/galacticus.log")
failures = []
for logFile in logFiles:
    result = subprocess.run(
        f'grep -q -i -e fatal -e aborted -e "Galacticus experienced an error in the GSL library" {logFile}',
        shell=True
    )
    if result.returncode == 0:
        failures.append(logFile)

if failures:
    for failure in failures:
        print(f"FAILED: log from {failure}:")
        with open(failure) as f:
            print(f.read())
else:
    print("SUCCESS: model run")

# Read test and reference data.
modelData = {}
i = -1
for model in models:
    if "parameterFile" in model:
        i += 1
        tmpName = f"outputs/test-merger-tree-builder/galacticus_{i}:1/galacticus.hdf5"
        if os.path.exists(tmpName):
            shutil.copy(tmpName, model["fileName"])
            model["exists"] = True
        else:
            model["exists"] = False
            continue
    else:
        model["exists"] = True

    if not model["exists"]:
        continue

    with h5py.File(model["fileName"], "r") as f:
        group = f["conditionalMassFunction"]
        modelData[model["type"]] = {}
        for datasetName in ("massParent", "massRatio", "redshiftParent", "redshiftProgenitor",
                            "conditionalMassFunction", "conditionalMassFunctionError"):
            modelData[model["type"]][datasetName] = group[datasetName][:]
        modelData[model["type"]]["matchMeasureMaximum"]        = 0.0
        modelData[model["type"]]["deviationFractionalMaximum"] = 0.0

# Iterate over redshifts and parent masses.
subprocess.run("mkdir -p outputs/mergerTreeBuilder", shell=True)
if "reference" in modelData:
    nRedshift = len(modelData["reference"]["redshiftParent"])
    nMass     = len(modelData["reference"]["massParent"])
    for iRedshift in range(nRedshift):
        for iMass in range(nMass):
            for model in models:
                if not model.get("exists") or "reference" not in model:
                    continue
                refType  = model["reference"]
                modType  = model["type"]
                if refType not in modelData or modType not in modelData:
                    continue
                refCMF   = modelData[refType]["conditionalMassFunction"][iRedshift, iMass, :]
                modCMF   = modelData[modType]["conditionalMassFunction"][iRedshift, iMass, :]
                refErr   = modelData[refType]["conditionalMassFunctionError"][iRedshift, iMass, :]
                modErr   = modelData[modType]["conditionalMassFunctionError"][iRedshift, iMass, :]
                select   = np.where((refCMF > 0.0) & (modCMF > 0.0))[0]
                if len(select) == 0:
                    continue
                matchMeasure        = (modCMF[select] - refCMF[select])**2 / (modErr[select]**2 + refErr[select]**2)
                deviationFractional = np.abs(modCMF[select] - refCMF[select]) / np.maximum(modCMF[select], refCMF[select])
                if matchMeasure.max() > modelData[modType]["matchMeasureMaximum"]:
                    modelData[modType]["matchMeasureMaximum"] = matchMeasure.max()
                if deviationFractional.max() > modelData[modType]["deviationFractionalMaximum"]:
                    modelData[modType]["deviationFractionalMaximum"] = deviationFractional.max()

# Report on the maximum deviations.
typeLength = max(len(m["type"]) for m in models)
for model in models:
    if not model.get("exists") or "reference" not in model:
        continue
    modType = model["type"]
    if modType not in modelData:
        continue
    matchMeasMax  = modelData[modType]["matchMeasureMaximum"]
    devFracMax    = modelData[modType]["deviationFractionalMaximum"]
    status        = "FAILED" if (matchMeasMax > model["toleranceSigma"] and devFracMax > model["toleranceFractional"]) else "success"
    padding       = " " * (typeLength - len(modType))
    print(f"{status}: {modType}{padding}: \u0394\u03c3\u2098\u2090\u2093 = {matchMeasMax:3.1f} {{limit: {model['toleranceSigma']:3.1f}}}; \u0394(log f)\u2098\u2090\u2093 = {devFracMax:7.1e} {{limit: {model['toleranceFractional']:7.1e}}}")
