import os

# Ensure the outputs directory exists before any log files are opened.
os.makedirs("outputs", exist_ok=True)
#!/usr/bin/env python3
import subprocess
import sys
import os
import glob
import argparse
import time

# Master test orchestrator - runs all test scripts.
# Andrew Benson (ported to Python)

parser = argparse.ArgumentParser()
parser.add_argument("--galacticusPath", type=str, default="..")
args, _ = parser.parse_known_args()

# Change to testSuite directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Find all test scripts.
testScripts = sorted(glob.glob("test-*.py"))
# Exclude self.
testScripts = [s for s in testScripts if s != "test-all.py"]

overallStatus = "SUCCESS"
results = []

for script in testScripts:
    print(f"\n==> Running {script}...")
    startTime = time.time()
    with open(f"outputs/{os.path.basename(script).replace('.py', '.log')}", "w") as logFile:
        status = subprocess.run(
            f"python3 {script}",
            shell=True, stdout=logFile, stderr=subprocess.STDOUT
        )
    elapsed = time.time() - startTime
    # Judge the outcome from the log. A test signals failure by printing a line containing "FAILED", and a whole-test skip
    # (an unmet prerequisite, such as too few MPI processes) by printing a line beginning "SKIPPED". A non-zero exit status is
    # also a failure - it means the script died before it could print its own marker.
    logPath = f"outputs/{os.path.basename(script).replace('.py', '.log')}"
    failed  = subprocess.run(f"grep -q -e FAIL -e FAILED {logPath}"  , shell=True).returncode == 0
    skipped = subprocess.run(f"grep -q '^SKIPPED' {logPath}"         , shell=True).returncode == 0
    if failed or status.returncode != 0:
        testStatus = "FAILED"
        overallStatus = "FAILED"
    elif skipped:
        testStatus = "SKIPPED"
    else:
        testStatus = "PASSED"
    results.append((script, testStatus, elapsed))
    print(f"    {testStatus}: {script} ({elapsed:.1f}s)")

print("\n\n=== Test Summary ===")
for script, testStatus, elapsed in results:
    print(f"  {testStatus:8s}: {script} ({elapsed:.1f}s)")

skippedCount = sum(1 for _, testStatus, _ in results if testStatus == "SKIPPED")
if skippedCount > 0:
    print(f"\n{skippedCount} test(s) were skipped because a prerequisite was not met - they tested nothing.")

print(f"\nOverall: {overallStatus}")
