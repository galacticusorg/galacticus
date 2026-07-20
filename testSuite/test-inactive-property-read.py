#!/usr/bin/env python3
import subprocess
import sys
import os

# Test that the run-time check for reads of inactive property values during active-property differential evolution
# (issue #128) fires when such a read occurs. The model uses a Jacobian-based ODE solver with the disk stellar
# luminosities solved as inactive properties, together with a test-only node operator
# ("testInactivePropertyRead") that deliberately reads those (inactive) luminosities during evaluation of the
# derivatives of the active properties. The check is emitted only in debugging (-DDEBUGGING) builds, so the test
# is skipped when the executable was not built with debugging enabled.
# Andrew Benson

# The guard's error message is compiled into the executable only in debugging builds. Detect whether debugging was
# enabled by searching the executable for that (compiled-in) message.
debugging = (
    subprocess.run(
        "grep -q -a 'value of inactive property' ../Galacticus.exe",
        shell=True
    ).returncode == 0
)

# In a non-debugging build the check is compiled out, so there is nothing to test - skip without running the model.
if not debugging:
    print("SKIPPED: executable was not built with debugging enabled, so the inactive-property-read check is absent")
    sys.exit(0)

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Run the model.
with open("outputs/inactivePropertyRead.log", "w") as logFile:
    status = subprocess.run(
        "cd ..; export OMP_NUM_THREADS=1; ./Galacticus.exe testSuite/parameters/inactivePropertyRead.xml",
        shell=True, stdout=logFile, stderr=subprocess.STDOUT
    )

# Check that the expected error message was given.
messageGiven = (
    subprocess.run(
        "grep -q 'value of inactive property \"luminositiesStellar\"' outputs/inactivePropertyRead.log",
        shell=True
    ).returncode == 0
)

if messageGiven:
    print("SUCCESS: read of inactive property during active-property evolution was detected")
elif status.returncode == 0:
    print("FAILED: read of inactive property during active-property evolution was not detected (model ran to completion)")
else:
    print("FAILED: model aborted but the expected inactive-property-read error message was not given")
