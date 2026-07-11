#!/usr/bin/env python3
import subprocess

# Test parameter file differences tool.
# Andrew Benson (07-September-2024)

# `parametersDiff.py` follows the diff exit convention: 0 = no differences, 1 = differences,
# >= 2 = trouble (e.g. the xdiff download failed). Interpreting any non-zero code as
# "differences found" would let an infrastructure error masquerade as a pass (or a failure)
# in the cases below, so match the expected codes exactly.

print("Respect order; '.4f' format")
status = subprocess.run("../scripts/parameters/parametersDiff.py --canonicalizeValues .4f --respectOrder parameters/parametersDiff1.xml parameters/parametersDiff2.xml",shell=True)
if status.returncode == 1:
    result = "SUCCESS"
elif status.returncode == 0:
    result = "FAILED"
else:
    result = "FAILED (parametersDiff.py errored with exit code "+str(status.returncode)+")"
print(result+": differences expected")

print("Ignore order; '.4f' format")
status = subprocess.run("../scripts/parameters/parametersDiff.py --canonicalizeValues .4f parameters/parametersDiff1.xml parameters/parametersDiff2.xml",shell=True)
if status.returncode == 1:
    result = "SUCCESS"
elif status.returncode == 0:
    result = "FAILED"
else:
    result = "FAILED (parametersDiff.py errored with exit code "+str(status.returncode)+")"
print(result+": differences expected")

print("Ignore order; '.3f' format")
status = subprocess.run("../scripts/parameters/parametersDiff.py --canonicalizeValues .3f parameters/parametersDiff1.xml parameters/parametersDiff2.xml",shell=True)
if status.returncode == 0:
    result = "SUCCESS"
elif status.returncode == 1:
    result = "FAILED"
else:
    result = "FAILED (parametersDiff.py errored with exit code "+str(status.returncode)+")"
print(result+": no differences expected")
