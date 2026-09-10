#!/usr/bin/env python3
"""Test the recording of parameters which take default values (issue #1353).

Galacticus does not report the use of a default to the terminal. Defaulting is
an extremely common idiom - a small model defaults of order a hundred parameters
and as many classes - so reporting each one swamped the bounded list of warnings
that `Error_Report` replays on a fatal error, crowding genuine warnings out of it
entirely. The information is recorded in the output file instead, where each
defaulted parameter carries a companion attribute named for it plus the
`{defaulted}` suffix.

This test pins that contract from both sides:

  * No default-related message appears on the terminal, even at `warn` verbosity.
  * Both kinds of default are marked in the output file - a parameter whose
    *value* was defaulted, and a `functionClass` parameter whose *class* was not
    specified and so was built from its default.
  * A parameter given explicitly in the parameter file is not marked, so the
    absence of a marker means the value came from the user.
  * Every marker annotates a parameter which is itself present, and the set of
    markers does not depend on the number of threads.

Andrew Benson
"""

import os
import re
import subprocess
import sys

import h5py

# The suffix marking a defaulted parameter in the output file. This must match
# `defaultedMarkerSuffix` in `source/utility/input_parameters.F90`.
defaultedMarkerSuffix = "{defaulted}"

parameterFile = "testSuite/parameters/parametersExtract.xml"
changeFile    = "testSuite/parameters/defaultParameterReporting.xml"
outputFile    = "outputs/defaultParameterReporting.hdf5"
logFileName   = "outputs/test-parameters-default-reporting.log"

subprocess.run("mkdir -p outputs",shell=True)

failed = False


def runModel(threadCount):
    """Run the model on the given number of threads, returning its output."""
    subprocess.run(f"rm -f {outputFile}",shell=True)
    process = subprocess.run(
        f"export OMP_NUM_THREADS={threadCount}; cd ..; "
        f"./Galacticus.exe {parameterFile} {changeFile}",
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True,
    )
    return process.returncode, process.stdout or ""


def readMarkers():
    """Return {group: (values, markers)} for the output file's Parameters group."""
    groups = {}

    def visit(name,node):
        if not isinstance(node,h5py.Group):
            return
        values  = set()
        markers = set()
        for attributeName in node.attrs.keys():
            if attributeName.endswith(defaultedMarkerSuffix):
                markers.add(attributeName[:-len(defaultedMarkerSuffix)])
            else:
                values.add(attributeName)
        groups[name] = (values,markers)

    with h5py.File(outputFile,"r") as outputFileHDF5:
        parameters = outputFileHDF5["Parameters"]
        visit("parameters",parameters)
        parameters.visititems(visit)
    return groups


# Run the model. The thread count matters: the record of defaulted parameters is
# shared rather than copied per thread, so a per-thread copy must not be able to
# change what is recorded.
print("Running model...")
returnCode,output = runModel(4)
with open(logFileName,"w") as logFile:
    logFile.write(output)
print("...done")

if returnCode != 0:
    print(f"FAILED: model did not run to completion (return code {returnCode}) - see {logFileName}")
    failed = True

# No default-related message may be emitted, even at `warn` verbosity.
if not failed:
    stray = re.findall(r"Using default (?:value|class) for parameter '\[[^\]]+\]'",output)
    if stray:
        print(f"FAILED: {len(stray)} default-related messages were emitted to the terminal; these are recorded in the output file and must not be reported as warnings")
        failed = True
    else:
        print("SUCCESS: no default-related messages were emitted to the terminal")

if not failed and not os.path.exists(outputFile):
    print(f"FAILED: output file {outputFile} was not produced")
    failed = True

if not failed:
    groups = readMarkers()

    # Every marker must annotate a parameter which is itself present.
    markerCount = sum(len(markers) for _,markers in groups.values())
    orphaned    = [f"{group}/{marker}"
                   for group,(values,markers) in groups.items()
                   for marker in markers if marker not in values]
    if markerCount == 0:
        print("FAILED: no defaulted parameters were marked in the output file's `Parameters` group")
        failed = True
    elif orphaned:
        print(f"FAILED: {len(orphaned)} markers annotate no parameter of their own. First few: {orphaned[:5]}")
        failed = True
    else:
        print(f"SUCCESS: {markerCount} defaulted parameters are marked, each annotating a parameter that is present")

# A defaulted *class* must be marked. `linearGrowth` is not specified in the
# parameter file, so its class is built from the default - the same default that
# `test-parameters-extract.py` checks appears in an extracted parameter file.
if not failed:
    topLevelValues,topLevelMarkers = groups["parameters"]
    defaultedClass = "linearGrowth"
    if defaultedClass not in topLevelValues:
        print(f"FAILED: expected defaulted class `{defaultedClass}` is absent from the output file")
        failed = True
    elif defaultedClass not in topLevelMarkers:
        print(f"FAILED: class `{defaultedClass}` was built from its default but is not marked as defaulted")
        failed = True
    else:
        print(f"SUCCESS: defaulted class `{defaultedClass}` is marked in the output file")

# A parameter given explicitly in the parameter file must carry no marker.
if not failed:
    explicit = "componentBasic"
    if explicit not in topLevelValues:
        print(f"FAILED: expected parameter `{explicit}` is absent from the output file, so the marking of explicitly-set parameters cannot be checked")
        failed = True
    elif explicit in topLevelMarkers:
        print(f"FAILED: parameter `{explicit}` is set explicitly in the parameter file but is marked as having taken a default value")
        failed = True
    else:
        print(f"SUCCESS: explicitly-set parameter `{explicit}` carries no defaulted marker")

# The set of markers must not depend on the thread count.
if not failed:
    markersMultiThreaded = {(group,marker)
                            for group,(_,markers) in groups.items()
                            for marker in markers}
    returnCodeSingle,_ = runModel(1)
    if returnCodeSingle != 0:
        print(f"FAILED: single-threaded model did not run to completion (return code {returnCodeSingle})")
        failed = True
    else:
        markersSingleThreaded = {(group,marker)
                                 for group,(_,markers) in readMarkers().items()
                                 for marker in markers}
        if markersSingleThreaded == markersMultiThreaded:
            print(f"SUCCESS: the {len(markersMultiThreaded)} markers are identical on 1 and 4 threads")
        else:
            difference = markersSingleThreaded ^ markersMultiThreaded
            print(f"FAILED: the markers differ between 1 and 4 threads ({len(difference)} differences). First few: {sorted(difference)[:5]}")
            failed = True

# Failure is signaled by the text above, not by the exit status.
if failed:
    print("FAILED: default parameter recording")
else:
    print("SUCCESS: default parameter recording")
sys.exit(0)
