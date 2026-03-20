#!/usr/bin/env python3
import subprocess
import sys
import os
import glob

# Run a set of short Galacticus models which test cases that have failed before,
# in order to catch regressions.
# Andrew Benson (ported to Python)

# Indicate that this test can manage its own jobs.
# selfManage: true

outputDirectory = "outputs/regressions"
subprocess.run(f"mkdir -p {outputDirectory}", shell=True)

# Find all regression parameter files and run them.
overallStatus = "SUCCESS"
for filePath in sorted(glob.glob("regressions/**/*.xml", recursive=True) + glob.glob("regressions/**/*.py", recursive=True)):
    status      = None
    logFilePath = None
    if filePath.endswith(".xml"):
        print(f"Running regression: {filePath}")
        logFilePath = f"{outputDirectory}/{os.path.basename(filePath).replace('.xml', '')}.log"
        with open(logFilePath, "w") as logFile:
            status = subprocess.run(
                f"cd ..; ./Galacticus.exe testSuite/{filePath}",
                shell=True, stdout=logFile, stderr=subprocess.STDOUT
            )
    elif filePath.endswith(".py"):
        print(f"Running regression script: {filePath}")
        logFilePath = f"{outputDirectory}/{os.path.basename(filePath).replace('.py', '')}.log"
        with open(logFilePath, "w") as logFile:
            status = subprocess.run(
                f"cd ..; python3 testSuite/{filePath}",
                shell=True, stdout=logFile, stderr=subprocess.STDOUT
            )
    result1 = subprocess.run(f"grep -q -i -e fatal -e aborted {logFilePath}", shell=True)
    result2 = subprocess.run(f"grep -q FAIL {logFilePath}"                  , shell=True)
    if result1.returncode == 0 or result2.returncode == 0 or status.returncode != 0:
        print(f"FAILED: regression '{filePath}'")
        with open(logFilePath) as f:
            print(f.read())
            overallStatus = "FAILED"
    else:
        print(f"SUCCESS: regression '{filePath}'")

print(f"{overallStatus}: regressions")
