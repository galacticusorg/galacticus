#!/usr/bin/env python3
import subprocess
import sys
import os
import re
import glob
import argparse
import xml.etree.ElementTree as ET

# Check for memory leaks in various models.
# Andrew Benson (ported to Python)

# Extract options.
parser = argparse.ArgumentParser()
parser.add_argument("--mpi", type=int, default=0)
args, _ = parser.parse_known_args()

# Models to check.
models = [
    {
        "label":      "formationHalos",
        "parameters": "testSuite/parameters/memoryLeakFormationHalos.xml",
    },
    {
        "label":      "MCMC",
        "parameters": "testSuite/parameters/mcmcConfig.xml",
        "mpi": {
            "processes": 4,
            "threads":   1,
        },
    },
]

overallStatus = "SUCCESS"

for model in models:
    hasMPI = "mpi" in model
    # Skip models that do not match MPI option.
    if args.mpi and not hasMPI:
        continue
    if not args.mpi and hasMPI:
        continue

    print(f"Running model '{model['label']}'...")
    subprocess.run(f"mkdir -p outputs/memoryLeaks/{model['label']}", shell=True)

    if hasMPI:
        mpiPrefix = f"export OMP_NUM_THREADS={model['mpi']['threads']}; mpirun --oversubscribe --allow-run-as-root --n {model['mpi']['processes']} "
    else:
        mpiPrefix = ""

    cmd = (
        f"cd ..; {mpiPrefix}valgrind --leak-check=full --xml=yes "
        f"--xml-file=testSuite/outputs/memoryLeaks/{model['label']}/memory-leaks-%p.xml "
        f"./Galacticus.exe {model['parameters']}"
    )
    status = subprocess.run(cmd, shell=True)
    if status.returncode != 0:
        print("\tFAILED:  model run:")
    else:
        print("\tSUCCESS: model run")

    # Find the Valgrind outputs.
    leakDir = f"outputs/memoryLeaks/{model['label']}"
    valgrindOutputs = glob.glob(f"{leakDir}/memory-leaks-*.xml")

    for valgrindOutput in valgrindOutputs:
        print(f"\tAnalyzing output file '{os.path.basename(valgrindOutput)}'")
        try:
            tree = ET.parse(valgrindOutput)
            root = tree.getroot()
        except ET.ParseError:
            print("   malformed XML")
            continue

        for error in root.findall("error"):
            kind = error.findtext("kind", "")
            if kind not in ("Leak_PossiblyLost", "Leak_DefinitelyLost"):
                continue

            ignore    = False
            frameFinal = None
            stack = error.find("stack")
            if stack is not None:
                for frame in stack.findall("frame"):
                    obj = frame.findtext("obj", "")
                    if frameFinal is None and "Galacticus.exe" in obj:
                        frameFinal = frame
                    if kind == "Leak_PossiblyLost" and "libgomp" in obj:
                        ignore = True
                    if "libmpi" in obj:
                        ignore = True

            if ignore or frameFinal is None:
                continue

            fFile = frameFinal.findtext("file", "(UNKNOWN)")
            fLine = frameFinal.findtext("line", "(UNKNOWN)")
            print(f"\t\tMemory leak ({kind}) for model '{model['label']}' in: '{fFile}' line {fLine}")
            overallStatus = "FAILED"

print(f"{overallStatus}: memory leaks")
