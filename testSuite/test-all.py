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
    # Check log for FAIL.
    logPath = f"outputs/{os.path.basename(script).replace('.py', '.log')}"
    result  = subprocess.run(f"grep -q -e FAIL -e FAILED {logPath}", shell=True)
    if result.returncode == 0 or status.returncode != 0:
        testStatus = "FAILED"
        overallStatus = "FAILED"
    else:
        testStatus = "PASSED"
    results.append((script, testStatus, elapsed))
    print(f"    {testStatus}: {script} ({elapsed:.1f}s)")

print("\n\n=== Test Summary ===")
for script, testStatus, elapsed in results:
    print(f"  {testStatus:8s}: {script} ({elapsed:.1f}s)")

print(f"\nOverall: {overallStatus}")
