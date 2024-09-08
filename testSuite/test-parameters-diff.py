#!/usr/bin/env python3
import subprocess

# Test parameter file differences tool.
# Andrew Benson (07-September-2024)

print("Respect order; '.4f' format")
status = subprocess.run("../scripts/parameters/parametersDiff.py --canonicalizeValues .4f --respectOrder parameters/parametersDiff1.xml parameters/parametersDiff2.xml",shell=True)
result = "FAILED" if status.returncode == 0 else "SUCCESS"
print(result+": differences expected")

print("Ignore order; '.4f' format")
status = subprocess.run("../scripts/parameters/parametersDiff.py --canonicalizeValues .4f parameters/parametersDiff1.xml parameters/parametersDiff2.xml",shell=True)
result = "FAILED" if status.returncode == 0 else "SUCCESS"
print(result+": differences expected")

print("Ignore order; '.3f' format")
status = subprocess.run("../scripts/parameters/parametersDiff.py --canonicalizeValues .3f parameters/parametersDiff1.xml parameters/parametersDiff2.xml",shell=True)
result = "FAILED" if status.returncode != 0 else "SUCCESS"
print(result+": no differences expected")
