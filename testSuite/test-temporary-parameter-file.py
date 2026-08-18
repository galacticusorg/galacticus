#!/usr/bin/env python3
"""Test that the temporary parameter file is removed from the file system (issue #1391).

Every run accumulates its output parameter set in a temporary HDF5 file,
`glcTmpPar.<pid>.<n>`, placed on a memory-backed file system (`/dev/shm`, or
`/tmp` on macOS) and copied into the output file at the end of the run. The file
is removed from the file system as soon as it has been opened - it remains fully
usable through the open file handle, but has no name in the file system, so it
cannot outlive the run even if the run is killed.

When this failed to happen the files accumulated silently, one per run, until the
file system was full - after which *every* subsequent run died at startup while
closing the temporary file, an error which gives no hint that its cause is the
leftovers of earlier runs.

This test runs a model and requires that it leave no temporary parameter file of
its own behind.

Andrew Benson
"""

import glob
import os
import platform
import subprocess
import sys

# The temporary path and file name root must match those used in
# `inputParametersConstructorNode` in `source/utility/input_parameters.F90`.
temporaryPath     = "/tmp" if platform.system() == "Darwin" else "/dev/shm"
temporaryRootName = "glcTmpPar"
temporaryPattern  = os.path.join(temporaryPath,temporaryRootName+".*")

parameterFile = "parameters/quickTest.xml"
changeFile    = "testSuite/parameters/temporaryParameterFile.xml"
outputFile    = "outputs/temporaryParameterFile.hdf5"
logFileName   = "outputs/test-temporary-parameter-file.log"

# Change to the testSuite directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))
subprocess.run("mkdir -p outputs",shell=True)
subprocess.run(f"rm -f {outputFile}",shell=True)

failed = False

# Run the model, recording the temporary files present before and after. The model is run without an intervening shell so that
# the process ID we see is that of the model itself - the temporary file is named for it, which lets us distinguish files left by
# this run from any belonging to other runs sharing the same machine.
filesBefore = set(glob.glob(temporaryPattern))
print("Running model...")
with open(logFileName,"w") as logFile:
    process = subprocess.Popen(
        ["./Galacticus.exe",parameterFile,changeFile],
        cwd="..", stdout=logFile, stderr=subprocess.STDOUT,
    )
    returnCode = process.wait()
print("...done")
filesAfter = set(glob.glob(temporaryPattern))


def processID(fileName):
    """Return the process ID recorded in the name of a temporary parameter file."""
    fields = os.path.basename(fileName).split(".")
    return fields[1] if len(fields) > 1 else None


filesNew    = filesAfter-filesBefore
filesLeaked = [fileName for fileName in filesNew if processID(fileName) == str(process.pid)]
filesOther  = filesNew-set(filesLeaked)

if returnCode != 0:
    print(f"FAILED: model did not run to completion (return code {returnCode}) - see {logFileName}")
    failed = True
elif not os.path.exists(outputFile):
    print(f"FAILED: output file {outputFile} was not produced")
    failed = True

if filesLeaked:
    print(f"FAILED: the model left {len(filesLeaked)} temporary parameter file(s) of its own in {temporaryPath}: {sorted(filesLeaked)}")
    failed = True
    # Do not leave our own mess behind, even on failure.
    for fileName in filesLeaked:
        os.remove(fileName)
else:
    print(f"SUCCESS: the model left no temporary parameter file of its own in {temporaryPath}")

# Files appearing during the run but belonging to some other process are not evidence of a problem here - they are simply other
# models running on the same machine - so they are reported but not judged.
if filesOther:
    print(f"NOTE: {len(filesOther)} temporary parameter file(s) belonging to other processes appeared in {temporaryPath} during this test; these are not the responsibility of this run")

# The output parameter set must survive the removal of the temporary file which carried it.
if not failed:
    import h5py
    with h5py.File(outputFile,"r") as outputFileHDF5:
        if "Parameters" not in outputFileHDF5:
            print("FAILED: the output file contains no `Parameters` group, so the content of the temporary file was lost")
            failed = True
        elif len(outputFileHDF5["Parameters"].attrs) == 0:
            print("FAILED: the `Parameters` group in the output file is empty, so the content of the temporary file was lost")
            failed = True
        else:
            print("SUCCESS: the parameter set accumulated in the temporary file was copied to the output file")

# Failure is signaled by the text above, not by the exit status.
if failed:
    print("FAILED: temporary parameter file handling")
else:
    print("SUCCESS: temporary parameter file handling")
sys.exit(0)
