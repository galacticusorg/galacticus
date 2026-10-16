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
    {
        "label":      "MCMCHMF",
        "parameters": "testSuite/parameters/memoryLeakMCMCHMFConfig.xml",
        "mpi": {
            "processes": 4,
            "threads":   1,
        },
    }
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
            # Consider only the two kinds which indicate a leak of memory that can no longer be
            # reached. `Leak_StillReachable` is excluded because a pointer to the allocation still
            # exists at exit, which is the normal state of anything allocated once and held for the
            # lifetime of the run; `Leak_IndirectlyLost` is excluded because it is a consequence of
            # some `Leak_DefinitelyLost` block which is itself reported.
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
                    # Ignore possibly-lost allocations made beneath libgomp.
                    #
                    # Read this filter as broader than it looks. Under gfortran the body of an
                    # `!$omp parallel` region is *outlined* into a separate function (`..._omp_fn.N`)
                    # which is invoked through `GOMP_parallel`, and worker threads run beneath
                    # `gomp_thread_start` — both of which live in libgomp. Essentially every
                    # allocation made inside any parallel region therefore carries a libgomp frame in
                    # its allocation stack, not merely the threadprivate allocations this filter was
                    # presumably introduced for. In effect, then, *no* possibly-lost allocation made
                    # by parallel code is reported.
                    #
                    # That is not necessarily wrong: "possibly lost" means valgrind found only an
                    # interior pointer to the block, which is routine for the descriptors gfortran
                    # builds for threadprivate and allocatable data, so the reports it suppresses are
                    # largely false positives. What matters is knowing the extent of it, and that the
                    # check is gated on `kind`, so `Leak_DefinitelyLost` inside parallel regions is
                    # still reported and genuine leaks are not hidden wholesale.
                    if kind == "Leak_PossiblyLost" and "libgomp" in obj:
                        ignore = True
                    # Ignore anything allocated beneath the MPI library. Note that, unlike the
                    # libgomp filter above, this one is *not* gated on `kind`, so it suppresses
                    # definitely-lost reports as well as possibly-lost ones.
                    if "libmpi" in obj:
                        ignore = True

            # `frameFinal` is the innermost frame within Galacticus itself, used to report where the
            # allocation was made. A leak whose stack contains no Galacticus frame at all — one made
            # entirely within a library — is therefore skipped rather than reported against an
            # unknown location.
            if ignore or frameFinal is None:
                continue

            fFile = frameFinal.findtext("file", "(UNKNOWN)")
            fLine = frameFinal.findtext("line", "(UNKNOWN)")
            print(f"\t\tMemory leak ({kind}) for model '{model['label']}' in: '{fFile}' line {fLine}")
            print(error)
            overallStatus = "FAILED"

print(f"{overallStatus}: memory leaks")
