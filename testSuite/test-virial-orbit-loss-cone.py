#!/usr/bin/env python3
"""Regression test for the ``lossCone`` virial orbit class (issue #1321).

The class deadlocked and could not be run at all. The defects it carried were
reachable only once a model actually drove a tabulation, and they failed in
three distinct ways---an unbreakable hang, a segmentation fault, and floating
point exceptions---so this test runs a small model which uses the class and
checks for all three:

* the run is given a wall-clock timeout, so a re-introduced deadlock fails the
  test promptly rather than hanging the CI job until its own (multi-hour)
  timeout;
* the run must exit cleanly, with no error markers in its log;
* the resulting galaxy properties must all be finite, so that a re-introduced
  divide-by-zero which is *not* trapped cannot pass silently as a NaN.

The model is run multi-threaded because the class tabulates inside a nested
OpenMP parallel region, and that nesting is what surfaced both the original
deadlock and a cache-file collision between threads.

Andrew Benson (03-August-2026)
"""

import os
import subprocess
import sys

import h5py
import numpy as np

# The tabulation is coarse (see the parameter file), but it is still by far the dominant cost of this model, which takes
# around 15s locally. Allow ample headroom over that so that a slow runner does not produce a spurious failure, while still
# failing far sooner than the CI job timeout if the class deadlocks.
timeoutSeconds = 900

pathParameters      = os.path.abspath("parameters/test-virial-orbit-loss-cone.xml" )
pathOutputDirectory = os.path.abspath("outputs"                                    )
pathOutputLog       = os.path.abspath("outputs/test-virial-orbit-loss-cone.log"     )
pathOutputModel     = os.path.abspath("outputs/test-virial-orbit-loss-cone.hdf5"    )

os.makedirs(pathOutputDirectory,exist_ok=True)
if os.path.exists(pathOutputModel):
    os.remove(pathOutputModel)

# Run the model, under a timeout so that a deadlock is reported as a failure.
print("Running model...")
environment                    = os.environ.copy()
environment["OMP_NUM_THREADS"] = environment.get("OMP_NUM_THREADS","4")
with open(pathOutputLog,"w") as log:
    try:
        status = subprocess.run(
            f"cd ..; ./Galacticus.exe {pathParameters}",
            stdout  = log,
            stderr  = log,
            shell   = True,
            env     = environment,
            timeout = timeoutSeconds,
        )
    except subprocess.TimeoutExpired:
        print(f"FAILED: model run did not complete within {timeoutSeconds}s - the lossCone class may be deadlocked")
        sys.exit(0)
print("...done ("+str(status)+")")
if status.returncode != 0:
    print("FAILED: model run:")
    subprocess.run(f"cat {pathOutputLog}",shell=True)
    sys.exit(0)

# Check the log for errors. A floating point exception (which is how several of the defects in this class manifested) is
# reported through this path.
print("Checking for errors...")
status = subprocess.run(
    "grep -q -i"
    " -e fatal"
    " -e aborted"
    " -e \"floating point exception\""
    " -e \"Galacticus experienced an error in the GSL library\""
    f" {pathOutputLog}",
    shell=True,
)
if status.returncode == 0:
    print("FAILED: model run (errors):")
    subprocess.run(f"cat {pathOutputLog}",shell=True)
    sys.exit(0)
print("SUCCESS: model run")

# Validate the output. The orbital parameters drawn from this class feed the merging satellite calculations, so a
# non-finite orbit propagates into the galaxy properties. Checking for finiteness catches a re-introduced division by
# zero which happens not to be trapped.
countNodes = 0
countBad   = 0
with h5py.File(pathOutputModel,"r") as model:
    for outputName, output in model["Outputs"].items():
        if "nodeData" not in output:
            continue
        for propertyName, property_ in output["nodeData"].items():
            values = property_[:]
            if not np.issubdtype(values.dtype,np.floating):
                continue
            if propertyName == "nodeIndex":
                continue
            countNodes = max(countNodes,values.size)
            if not np.all(np.isfinite(values)):
                countBad += 1
                print(f"FAILED: non-finite values in {outputName}/nodeData/{propertyName}")

if countNodes == 0:
    print("FAILED: model produced no galaxies - the test is not exercising the lossCone class")
    sys.exit(0)
if countBad > 0:
    sys.exit(0)

print(f"SUCCESS: virial orbit lossCone ({countNodes} galaxies, all properties finite)")
