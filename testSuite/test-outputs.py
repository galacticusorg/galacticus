#!/usr/bin/env python3
import subprocess
import sys
import glob
import argparse

# Run a set of Galacticus models to test various output options.
# Andrew Benson (ported to Python)

# Parse command line options.
parser = argparse.ArgumentParser()
parser.add_argument("--launchMethod",  type=str, default="local")
parser.add_argument("--threadMaximum", type=int, default=1)
parser.add_argument("--ompThreads",    type=int, default=4)
parser.add_argument("--instance",      type=str, default=None)
args, _ = parser.parse_known_args()

launchOptions = f"--launchMethod {args.launchMethod} --threadMaximum {args.threadMaximum}"
if args.instance:
    launchOptions += f" --instance {args.instance}"

# Simply run the models.
subprocess.run(f"cd ..; scripts/aux/launch.pl testSuite/test-outputs.xml {launchOptions}", shell=True)

# Check for failed models.
logFiles = glob.glob("outputs/test-outputs/galacticus_*/galacticus.log")
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
