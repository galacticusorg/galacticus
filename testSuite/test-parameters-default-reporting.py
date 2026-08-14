#!/usr/bin/env python3
"""Test the reporting of parameters which take default values (issue #1353).

Two properties are pinned here, corresponding to the two defects that this test
was written alongside the fix for. They partly masked one another, so both must
be checked:

  * Each defaulted parameter is reported *exactly once*, however many threads the
    model runs on. The record of which parameters have already been reported is
    shared, rather than being deep-copied into the per-thread copies of the
    parameter set - were it copied, each thread would re-report parameters that
    another thread had already reported, and the message count would scale with
    the thread count. The model is therefore run multi-threaded here; on a single
    thread this check could not fail.

  * Parameters which share a container are reported *separately*. The record is
    keyed on the full parameter path including the parameter name; were it keyed
    on the containing node's path alone, all but the first parameter defaulted
    within any one container would be silently suppressed.

The marking of defaulted parameters in the output file's `Parameters` group is
also checked: a parameter which took a default value carries a companion
attribute named for it plus the `{defaulted}` suffix, and one given explicitly in
the parameter file carries no such attribute.

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
subprocess.run(f"rm -f {outputFile}",shell=True)

failed = False

# Run the model on several threads. The thread count matters: the duplicate
# reporting this test guards against was produced by the per-thread copies of the
# parameter set, so it is invisible on a single thread.
print("Running model...")
process = subprocess.run(
    "export OMP_NUM_THREADS=4; cd ..; "
    f"./Galacticus.exe {parameterFile} {changeFile}",
    shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    universal_newlines=True,
)
output = process.stdout or ""
with open(logFileName,"w") as logFile:
    logFile.write(output)
print("...done")

if process.returncode != 0:
    print(f"FAILED: model did not run to completion (return code {process.returncode}) - see {logFileName}")
    failed = True

# Collect the reported parameters. The message names the parameter by its full
# path, in brackets.
reported = re.findall(r"Using default value for parameter '\[([^\]]+)\]'",output)

if not failed and len(reported) == 0:
    print("FAILED: no default parameter values were reported - the test model is expected to default many parameters")
    failed = True

# Each parameter must be reported exactly once.
if not failed:
    repeated = sorted({parameter for parameter in reported if reported.count(parameter) > 1})
    if repeated:
        print(f"FAILED: {len(repeated)} parameters were reported more than once; the record of already-reported parameters is not being shared between copies of the parameter set. First few: {repeated[:5]}")
        failed = True
    else:
        print(f"SUCCESS: each of the {len(reported)} defaulted parameters was reported exactly once")

# Parameters sharing a container must be reported separately. Split each reported
# path into its containing path and the parameter name, and require that at least
# one container contributes more than one parameter - which cannot happen if the
# containing path alone is used to identify a parameter.
if not failed:
    containers = {}
    for parameter in reported:
        container,_,name = parameter.rpartition("/")
        containers.setdefault(container,[]).append(name)
    shared = {container: names for container,names in containers.items() if len(names) > 1}
    if shared:
        print(f"SUCCESS: parameters sharing a container are reported separately ({len(shared)} containers report more than one parameter)")
    else:
        print("FAILED: no container reported more than one parameter, which suggests that parameters sharing a container are colliding and all but the first are being suppressed")
        failed = True

# Check the marking of defaulted parameters in the output file.
if not failed:
    if not os.path.exists(outputFile):
        print(f"FAILED: output file {outputFile} was not produced")
        failed = True

if not failed:
    with h5py.File(outputFile,"r") as outputFileHDF5:
        parameters = outputFileHDF5["Parameters"]

        # Every marker must annotate a parameter which is itself present.
        markerCount = 0
        orphaned    = []

        def checkGroup(name,node):
            global markerCount
            if not isinstance(node,h5py.Group):
                return
            for attributeName in node.attrs.keys():
                if attributeName.endswith(defaultedMarkerSuffix):
                    markerCount += 1
                    marked = attributeName[:-len(defaultedMarkerSuffix)]
                    if marked not in node.attrs:
                        orphaned.append(f"{name}/{attributeName}")

        checkGroup("parameters",parameters)
        parameters.visititems(checkGroup)

        if markerCount == 0:
            print("FAILED: no defaulted parameters were marked in the output file's `Parameters` group")
            failed = True
        elif orphaned:
            print(f"FAILED: {len(orphaned)} defaulted markers annotate no parameter of their own. First few: {orphaned[:5]}")
            failed = True
        else:
            print(f"SUCCESS: {markerCount} defaulted parameters are marked in the output file")

        # A parameter given explicitly in the parameter file must carry no marker.
        explicit = "componentBasic"
        if explicit not in parameters.attrs:
            print(f"FAILED: expected parameter `{explicit}` is absent from the output file, so the marking of explicitly-set parameters cannot be checked")
            failed = True
        elif explicit+defaultedMarkerSuffix in parameters.attrs:
            print(f"FAILED: parameter `{explicit}` is set explicitly in the parameter file but is marked as having taken a default value")
            failed = True
        else:
            print(f"SUCCESS: explicitly-set parameter `{explicit}` carries no defaulted marker")

# Failure is signaled by the text above, not by the exit status.
if failed:
    print("FAILED: default parameter reporting")
else:
    print("SUCCESS: default parameter reporting")
sys.exit(0)
