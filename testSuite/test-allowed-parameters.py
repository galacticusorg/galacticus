#!/usr/bin/env python3
import subprocess
import sys
import re

# Run a Galacticus model to test allowed parameter functionality.
# Andrew Benson (ported to Python)

# Run the test model.
subprocess.run("mkdir -p outputs", shell=True)
with open("outputs/test-allowed-parameters.log", "w") as logFile:
    subprocess.run(
        "cd ..; ./Galacticus.exe testSuite/parameters/test-allowed-parameters.xml",
        shell=True, stdout=logFile, stderr=subprocess.STDOUT
    )

# Parse the log file looking for disallowed parameters.
disallowed = []
with open("outputs/test-allowed-parameters.log") as log:
    for line in log:
        m = re.search(r'unrecognized parameter \[([a-zA-Z0-9\-]+) in [a-zA-Z\/]+\]', line)
        if m:
            disallowed.append(m.group(1))

# Extract a sorted list of unique disallowed parameters.
disallowed = sorted(set(disallowed))

# Check against expectations.
status = "success" if ":".join(disallowed) == "scaleCutOff" else "FAILURE"
if status == "FAILURE":
    print("  -> unrecognized parameters: " + ", ".join(disallowed))
print("Test allowed parameters functionality: " + status)
