#!/usr/bin/env python3
"""Characterization test for the *legitimate* (bounded) construction recursion
handled by ``recursive="yes"`` (issue #695).

The complementary test ``test-objectBuilder-recursion.py`` pins the #397
*abort* path for unintended, unbounded cycles. This test pins the opposite:
that a construction cycle a class explicitly opts in to (via ``recursive="yes"``)
builds successfully and runs end to end, exercising the machinery that the
issue #695 refactor will replace:

  * the recursive-build short-circuit that returns a shim wired to the object
    already under construction (percolation <-> halo-scale cycle),
  * a deepCopy of the recursive object graph -- forced by running
    ``evolveForests`` multi-threaded, which makes a per-thread deep copy of the
    node-operator tree (which reaches the recursive virial-density-contrast via
    the dark-matter halo scale),
  * a stateStore round-trip of the shim -- enabled by ``stateFileRoot``,
  * serialisation of the shim's descriptor into the output file's parameters.

It asserts the model:

  * runs to completion (return code 0 -- it must NOT hit the #397 abort, which
    would mean the legitimate cycle was mis-detected as unbounded),
  * does NOT die from a segmentation fault (return code 139 / negative signal),
  * does NOT emit the recursion-abort diagnostic, and
  * produces a non-empty output file with a parameters descriptor.

Andrew Benson
"""

import subprocess
import sys
import os

# Ensure output directories exist.
subprocess.run("mkdir -p outputs/recursionSuccessPercolation", shell=True)

parameterFile = "testSuite/parameters/recursionSuccessPercolation.xml"
outputFile    = "outputs/recursionSuccessPercolation.hdf5"

# Remove any stale output so the existence check below is meaningful.
subprocess.run(f"rm -f {outputFile}", shell=True)

# Run the model multi-threaded so that evolveForests makes per-thread deep
# copies of the object graph, exercising the recursive-shim deepCopy fix-up.
process = subprocess.run(
    "export OMP_NUM_THREADS=4; cd ..; "
    f"./Galacticus.exe {parameterFile}",
    shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    universal_newlines=True,
)
output     = process.stdout or ""
returnCode = process.returncode

failed = False

# A segmentation fault (unguarded recursion / dangling shim) reports 139 via the
# shell (128 + SIGSEGV), or a negative signal number if not run through a shell.
if returnCode == 139 or returnCode < 0:
    print(f"FAILED: model crashed (return code {returnCode}) - the recursive "
          f"build/deepCopy path is broken")
    print(output)
    failed = True

# The #397 abort diagnostic must NOT appear: this is a legitimate, bounded
# cycle that the recursive="yes" short-circuit is expected to handle.
if "recursive build of" in output or \
   "composites a member of its own class" in output:
    print("FAILED: the legitimate recursive build was mis-detected as an "
          "unbounded recursion and aborted")
    print(output)
    failed = True

# The model must run to completion.
if returnCode != 0:
    print(f"FAILED: model did not run to completion (return code {returnCode})")
    print(output)
    failed = True

# A successful run must have produced a non-empty output file.
if not failed:
    if not os.path.exists(outputFile):
        print(f"FAILED: output file {outputFile} was not produced")
        failed = True
    elif os.path.getsize(outputFile) == 0:
        print(f"FAILED: output file {outputFile} is empty")
        failed = True

# Confirm the output carries a parameters descriptor (the shim's descriptor was
# serialised without re-entering the cycle).
if not failed:
    try:
        import h5py
        with h5py.File(outputFile, "r") as f:
            if "Parameters" not in f:
                print("FAILED: output file has no 'Parameters' group "
                      "(descriptor serialisation failed)")
                failed = True
    except ImportError:
        print("WARNING: h5py not available; skipping descriptor check")

if failed:
    print("FAILED: object-builder successful-recursion test")
    sys.exit(1)

print("SUCCESS: object-builder successful-recursion test "
      f"(return code {returnCode})")
