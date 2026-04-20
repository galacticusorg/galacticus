#!/usr/bin/env python3
import subprocess
import sys
import glob
import shutil

# Export trees from Galacticus and check that they are written correctly.
# Andrew Benson (ported to Python)

# Run models.
subprocess.run(
    "cd ..; mkdir -p testSuite/outputs/test-merger-tree-write; "
    "./scripts/aux/launch.py testSuite/parameters/test-merger-tree-write.xml --launchMethod local --threadMaximum 1 --ompThreads 4; "
    "./scripts/aux/launch.py testSuite/parameters/test-merger-tree-write-secondary.xml --launchMethod local --threadMaximum 1 --ompThreads 4",
    shell=True
)

# Check for failed models.
logFiles = glob.glob("outputs/test-merger-tree-write/galacticus_*/galacticus.log")
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
    print("SUCCESS: model run")

# Validate the IRATE-format output.
validator = shutil.which("iratevalidate")
if validator:
    status = subprocess.run("iratevalidate outputs/test-merger-tree-write/exportedTreesIRATE.hdf5", shell=True)
    if status.returncode != 0:
        print("FAILED: IRATE-format file output by Galacticus did not validate")
else:
    print("SKIP: iratevalidate is not installed - validation of IRATE-format file will be skipped")
