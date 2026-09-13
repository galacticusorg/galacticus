#!/usr/bin/env python3
"""Regression test for the ``lossCone`` virial orbit class (issues #1321 and #1426).

The class deadlocked and could not be run at all. The defects it carried were
reachable only once a model actually drove a tabulation, and they failed in
three distinct ways---an unbreakable hang, a segmentation fault, and floating
point exceptions---so this test runs small models which use the class and
checks for all three:

* each run is given a wall-clock timeout, so a re-introduced deadlock fails the
  test promptly rather than hanging the CI job until its own (multi-hour)
  timeout;
* each run must exit cleanly, with no error markers in its log;
* the resulting galaxy properties must all be finite, so that a re-introduced
  divide-by-zero which is *not* trapped cannot pass silently as a NaN.

Two models are run. The first (issue #1321) tabulates at the epochs reached by
trees based at z=0, with a velocity grid chosen to sample the degenerate orbits
which must be special-cased. The second (issue #1426) is based at z=10, which
forces the tabulation to be built at an epoch so early that the halo mass
function and branching rate underflow to zero for the most massive hosts on the
(epoch-independent) mass lattice, making the environmental boost factor a 0/0
invalid operation.

Each model is run multi-threaded because the class tabulates inside a nested
OpenMP parallel region, and that nesting is what surfaced both the original
deadlock and a cache-file collision between threads.

Andrew Benson (03-August-2026)
"""

import os
import subprocess
import sys

import h5py
import numpy as np

# The tabulations are coarse (see the parameter files), but they are still by far the dominant cost of these models, which
# take around 15s and 105s respectively locally. Allow ample headroom over that so that a slow runner does not produce a
# spurious failure, while still failing far sooner than the CI job timeout if the class deadlocks.
timeoutSeconds = 900

models = [
    {
        "label"     : "default epoch",
        "parameters": "parameters/test-virial-orbit-loss-cone.xml",
        "model"     : "outputs/test-virial-orbit-loss-cone.hdf5",
        "log"       : "outputs/test-virial-orbit-loss-cone.log",
    },
    {
        "label"     : "early epoch",
        "parameters": "parameters/test-virial-orbit-loss-cone-early-epoch.xml",
        "model"     : "outputs/test-virial-orbit-loss-cone-early-epoch.hdf5",
        "log"       : "outputs/test-virial-orbit-loss-cone-early-epoch.log",
    },
]

pathOutputDirectory = os.path.abspath("outputs")
os.makedirs(pathOutputDirectory,exist_ok=True)


def modelRun(model):
    """Run one model, returning True if it completed cleanly."""
    pathParameters = os.path.abspath(model["parameters"])
    pathOutputLog  = os.path.abspath(model["log"       ])
    print("Running model ("+model["label"]+")...")
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
            print(f"FAILED: model run ({model['label']}) did not complete within {timeoutSeconds}s - the lossCone class may be deadlocked")
            return False
    print("...done ("+str(status)+")")
    if status.returncode != 0:
        print("FAILED: model run ("+model["label"]+"):")
        subprocess.run(f"cat {pathOutputLog}",shell=True)
        return False
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
        print("FAILED: model run ("+model["label"]+") (errors):")
        subprocess.run(f"cat {pathOutputLog}",shell=True)
        return False
    print("SUCCESS: model run ("+model["label"]+")")
    return True


def modelValidate(model):
    """Validate the output of one model, returning True if all properties are finite."""
    # The orbital parameters drawn from this class feed the merging satellite calculations, so a non-finite orbit propagates
    # into the galaxy properties. Checking for finiteness catches a re-introduced division by zero which happens not to be
    # trapped.
    pathOutputModel = os.path.abspath(model["model"])
    countNodes = 0
    countBad   = 0
    with h5py.File(pathOutputModel,"r") as modelFile:
        for outputName, output in modelFile["Outputs"].items():
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
                    print(f"FAILED: non-finite values in {outputName}/nodeData/{propertyName} ({model['label']})")
    if countNodes == 0:
        print("FAILED: model ("+model["label"]+") produced no galaxies - the test is not exercising the lossCone class")
        return False
    if countBad > 0:
        return False
    print(f"SUCCESS: virial orbit lossCone ({model['label']}; {countNodes} galaxies, all properties finite)")
    return True


for model in models:
    pathOutputModel = os.path.abspath(model["model"])
    if os.path.exists(pathOutputModel):
        os.remove(pathOutputModel)
    if modelRun(model):
        modelValidate(model)

sys.exit(0)
