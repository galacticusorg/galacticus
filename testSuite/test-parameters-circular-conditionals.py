#!/usr/bin/env python3
"""Check that circular `active` conditionals in a parameter file are reported, rather than looping forever.

A parameter may be made conditional on the value of another through its `active` attribute. Galacticus evaluates these
conditionals in passes over the parameter tree until every parameter's active state is known. With circular conditionals (A
active only if [B] has some value, and B active only if [A] has some value) no state can ever be settled, and the evaluation used
to loop forever (fixed October 2026). It now stops with an error.

This script runs quickTest.xml with:

  * two mutually conditional parameters, and checks that the model stops with the error (within a time limit, as the old
    behavior was to hang);
  * a chain of three conditionals, each depending on the next, listed in the reverse of the order in which they can be
    evaluated, so that several passes are needed - this must not be mistaken for a circular case, and the model must run.

Andrew Benson, Claude
"""
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

# This script runs with the working directory set to testSuite/, while models are run from the repository root.
pathOutputs = os.path.abspath("outputs/parametersCircularConditionals")
os.makedirs(pathOutputs, exist_ok=True)

# Time limit (in seconds) for each model. Evaluating the conditionals takes a negligible time, so a model which does not
# finish within this limit is taken to be looping.
timeLimit = 600

def buildParameters(name, conditionals):
    # Build a parameter file from quickTest.xml with the given conditional parameters added. Each conditional is given as a
    # tuple (name, value, condition), with condition None for an unconditional parameter. These parameters are read by nothing,
    # so are marked to suppress warnings about unrecognized parameters.
    tree = ET.parse("../parameters/quickTest.xml")
    root = tree.getroot()
    for element in root.findall("outputFileName"):
        root.remove(element)
    ET.SubElement(root, "outputFileName", value=os.path.join(pathOutputs, f"{name}.hdf5"))
    for nameParameter, value, condition in conditionals:
        attributes = {"value": value, "ignoreWarnings": "true"}
        if condition is not None:
            attributes["active"] = condition
        ET.SubElement(root, nameParameter, **attributes)
    fileParameters = os.path.join(pathOutputs, f"{name}.xml")
    tree.write(fileParameters)
    return fileParameters

def runModel(name, fileParameters):
    # Run a model, returning its exit status (None if the time limit was reached) and its log.
    try:
        status = subprocess.run(["./Galacticus.exe", fileParameters], cwd="..", capture_output=True, timeout=timeLimit)
        returnCode = status.returncode
        log        = status.stdout.decode()+status.stderr.decode()
    except subprocess.TimeoutExpired as e:
        returnCode = None
        log        = (e.stdout.decode() if e.stdout else "")+(e.stderr.decode() if e.stderr else "")
    fileLog = os.path.join(pathOutputs, f"{name}.log")
    with open(fileLog, "w") as logFile:
        logFile.write(log)
    return returnCode, log, fileLog

# Circular conditionals must be reported.
fileParameters = buildParameters(
    "circular",
    [
        ("conditionalTestA", "1", "[conditionalTestB] == 1"),
        ("conditionalTestB", "1", "[conditionalTestA] == 1")
    ]
)
returnCode, log, fileLog = runModel("circular", fileParameters)
if returnCode is None:
    print(f"FAILED: the model with circular conditionals did not stop within {timeLimit} s (log in {fileLog})")
elif returnCode == 0:
    print(f"FAILED: the model with circular conditionals ran, but should have been refused (log in {fileLog})")
elif "check for circular conditionals" in log:
    print("SUCCESS: circular conditionals are reported")
else:
    print(f"FAILED: the model with circular conditionals did not stop with the expected message (log in {fileLog}):")
    print(log[-1500:])

# A chain of conditionals needing several passes to evaluate must not be mistaken for circular conditionals.
fileParameters = buildParameters(
    "chain",
    [
        ("conditionalTestA", "1", "[conditionalTestB] == 1"),
        ("conditionalTestB", "1", "[conditionalTestC] == 1"),
        ("conditionalTestC", "1", None                     )
    ]
)
returnCode, log, fileLog = runModel("chain", fileParameters)
if returnCode is None:
    print(f"FAILED: the model with a chain of conditionals did not finish within {timeLimit} s (log in {fileLog})")
elif returnCode != 0 or "check for circular conditionals" in log:
    print(f"FAILED: the model with a chain of conditionals did not run (log in {fileLog}):")
    print(log[-1500:])
else:
    print("SUCCESS: a chain of conditionals is evaluated")
sys.exit(0)
