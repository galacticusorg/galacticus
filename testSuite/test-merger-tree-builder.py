#!/usr/bin/env python3
import subprocess
import sys
import os
import argparse
import shutil
import xml.etree.ElementTree as ET
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

# Build the parameterGrid XML and run via launch.py. The parameter files are parsed and their root elements grafted into the
# document, rather than their text being concatenated into it: each is a document in its own right and so begins with an XML
# declaration, which is well formed only at the very start of a document. Concatenating them produced a document which
# launch.py could not parse, so no model ran at all.
grid = ET.Element("parameterGrid")
for name, text in (
        ("emailReport",        "no"                                                                  ),
        ("doAnalysis",         "no"                                                                  ),
        ("modelRootDirectory", "testSuite/outputs/test-merger-tree-builder"                          ),
        ("baseParameters",     "testSuite/parameters/mergerTreeBuilderCole2000_intervalStepTrue.xml" ),
):
    ET.SubElement(grid, name).text = text
for paramFile in paramFiles:
    grid.append(ET.parse(paramFile).getroot())
ET.indent(grid)
ET.ElementTree(grid).write("outputs/test-merger-tree-builder.xml", encoding="UTF-8", xml_declaration=True)

subprocess.run(
    f"cd ..; mkdir -p testSuite/outputs/test-merger-tree-builder; ./scripts/aux/launch.py testSuite/outputs/test-merger-tree-builder.xml {launchOptions}",
    shell=True
)

# Identify which models this invocation actually ran. `launch.py` distributes models over instances - CI runs this script four
# times with `--instance N:4`, so each run launches only one of the four models - and it creates an output directory only for
# those models it runs. The presence of that directory, not of a log, is therefore what separates a model which ran from one
# belonging to another instance; a model which ran but died leaves the directory behind, with no log or no output in it.
modelIndex = -1
for model in models:
    if "parameterFile" not in model:
        continue
    modelIndex              += 1
    model["outputDirectory"] = f"outputs/test-merger-tree-builder/galacticus_{modelIndex}:1"
    model["launched"       ] = os.path.isdir(model["outputDirectory"])

launched = [model for model in models if model.get("launched")]

# Check for failed models. A model which never ran leaves no log to grep, so an empty set of logs would otherwise yield an empty
# set of failures and report as a success - hence the check that this invocation ran something at all, and that everything it did
# run left a log.
failures = []
if not launched:
    print(f"FAILED: none of the {len(paramFiles)} models was launched - no output directory was created")
else:
    for model in launched:
        logFile = f"{model['outputDirectory']}/galacticus.log"
        if not os.path.exists(logFile):
            print(f"FAILED: model '{model['type']}' was launched but produced no log at {logFile}")
            failures.append(logFile)
            continue
        result = subprocess.run(
            f'grep -q -i -e fatal -e aborted -e "Galacticus experienced an error in the GSL library" {logFile}',
            shell=True
        )
        if result.returncode == 0:
            print(f"FAILED: log from {logFile}:")
            with open(logFile) as f:
                print(f.read())
            failures.append(logFile)
    if not failures:
        print(f"SUCCESS: model run ({len(launched)} of {len(paramFiles)} models in this instance)")

# Read test and reference data.
modelData = {}
for model in models:
    if "parameterFile" in model:
        if not model["launched"]:
            # This model belongs to another instance - it is not expected to have run here, and is simply not compared.
            model["exists"] = False
            continue
        tmpName = f"{model['outputDirectory']}/galacticus.hdf5"
        if os.path.exists(tmpName):
            shutil.copy(tmpName, model["fileName"])
            model["exists"] = True
        else:
            # Report this: a model which ran but produced no output would otherwise be silently dropped from the comparisons
            # below, leaving the run to report as a success.
            print(f"FAILED: model '{model['type']}' was launched but produced no output at {tmpName}")
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
