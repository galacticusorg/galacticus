#!/usr/bin/env python3
import subprocess
import sys
import os
import glob
import argparse

# Run a set of short Galacticus models spanning a full range of method options.
# Andrew Benson (ported to Python)

# Parse command line options.
parser = argparse.ArgumentParser()
parser.add_argument("--launchMethod",  type=str, default="local")
parser.add_argument("--threadMaximum", type=int, default=1)
parser.add_argument("--ompThreads",    type=int, default=4)
args, _ = parser.parse_known_args()

launchOptions = f"--launchMethod {args.launchMethod} --threadMaximum {args.threadMaximum} --ompThreads {args.ompThreads}"

# Determine dynamic data path.
dataDynamicPath = os.environ.get("GALACTICUS_DYNAMIC_DATA_PATH", os.environ.get("GALACTICUS_DATA_PATH", ""))

# Remove automatically generated files to force them to be regenerated.
# Core files (older than 7 days).
for pattern in ("core.*", "vgcore.*"):
    subprocess.run(f"find ../ -name '{pattern}' -ctime +7 -exec rm {{}} \\;", shell=True)
# Others are checked only if we have the dynamic datasets path.
if dataDynamicPath != "":
    # FSPS stellar population synthesis code and associated file.
    subprocess.run(f"rm -f {dataDynamicPath}/dynamic/stellarPopulations/SSP_Spectra_Conroy-et-al_v2.5_imfSalpeter.hdf5", shell=True)
    subprocess.run(f"rm -rf {dataDynamicPath}/dynamic/FSPS_v2.5", shell=True)
    # Noninstantaneous recycling files (older than 14 days).
    for pattern in ("yield*.hdf5", "recycledFraction*.hdf5", "energyOutput*.hdf5"):
        subprocess.run(f"find {dataDynamicPath}/dynamic/stellarPopulations -name '{pattern}' -ctime +14 -exec rm {{}} \\;", shell=True)
    # CAMB transfer function files (older than 14 days).
    subprocess.run(f"find {dataDynamicPath}/dynamic/largeScaleStructure -name 'transfer_function_CAMB_*.xml' -ctime +14 -exec rm {{}} \\;", shell=True)

# Simply run the models.
subprocess.run(f"cd ..; scripts/aux/launch.pl testSuite/test-methods.xml {launchOptions}", shell=True)

# Check for failed models.
logFiles = glob.glob("outputs/test-methods/galacticus_*/galacticus.log")
failures = []
for logFile in logFiles:
    result = subprocess.run(f"grep -q -i -e fatal -e aborted {logFile}", shell=True)
    if result.returncode == 0:
        failures.append(logFile)

if failures:
    for failure in failures:
        print(f"FAILED: log from {failure}:")
        with open(failure) as f:
            print(f.read())
else:
    print("SUCCESS!")
