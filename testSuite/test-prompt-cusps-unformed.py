#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np
import lxml.etree as ET

# Test that halos in which no prompt cusp can form are pruned from merger trees.
#
# Prompt cusp formation is, by definition, the formation event of a halo. A halo for which the cusp collapse epoch would lie in
# its own future can therefore not have formed at all: the `darkMatterProfilePromptCusps` node operator assigns no cusp to such
# halos and labels them `promptCuspUnformed`, and a `pruneByFilter` merger tree operator removes them from the tree prior to
# evolution. See https://github.com/galacticusorg/galacticus/issues/1362.
#
# The model used here has a 3 keV thermal warm dark matter particle and a merger tree mass resolution of 5e6 M☉ - conditions
# under which such halos do arise. The model is run twice: once with the pruning operator removed, to establish that the
# condition is actually triggered under these conditions, and once with it in place, to establish that no halo lacking a cusp
# survives into the evolved tree.
# Andrew Benson

# Ensure the output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

parameterFileBase = "parameters/promptCuspsUnformed.xml"


def runModel(label, prune):
    """Run the model with pruning enabled or disabled, and return the cusp amplitudes of all output halos."""
    tree = ET.parse(parameterFileBase)
    root = tree.getroot()
    root.find("outputFileName").set("value", "testSuite/outputs/promptCuspsUnformed_" + label + ".hdf5")
    if not prune:
        # Remove the pruning operator, leaving the remainder of the operator sequence intact.
        operator = root.find(".//mergerTreeOperator[@value='pruneByFilter']")
        if operator is None:
            print("FAILED: could not find the `pruneByFilter` operator in " + parameterFileBase)
            sys.exit(0)
        operator.getparent().remove(operator)
    parameterFile = "outputs/promptCuspsUnformed_" + label + ".xml"
    logFile = "outputs/promptCuspsUnformed_" + label + ".log"
    tree.write(parameterFile)
    status = subprocess.run(
        "cd ..; ./Galacticus.exe testSuite/" + parameterFile + " > testSuite/" + logFile + " 2>&1",
        shell=True
    )
    if status.returncode != 0:
        print("FAILED: model '" + label + "' failed to run - see " + logFile)
        sys.exit(0)
    # Check the model log for failure markers.
    with open(logFile, "r") as log:
        logText = log.read()
    for marker in ("fatal", "aborted", "ODE integration failed", "unrecognized parameter"):
        if marker in logText:
            print("FAILED: model '" + label + "' reported '" + marker + "' - see " + logFile)
            sys.exit(0)
    # Extract, for every output halo, whether it carries the `promptCuspUnformed` label. The label is used in preference to the
    # cusp amplitude because labels persist through node promotion, whereas the cusp properties of a promoted node are overwritten
    # with those of its parent.
    labels = []
    counts = 0
    with h5py.File("outputs/promptCuspsUnformed_" + label + ".hdf5", "r") as model:
        outputs = model["Outputs"]
        for outputName in outputs:
            nodeData = outputs[outputName]["nodeData"]
            if "nodeLabelPromptCuspUnformed" not in nodeData:
                print("FAILED: model '" + label + "' output contains no `nodeLabelPromptCuspUnformed` dataset")
                sys.exit(0)
            labels.append(nodeData["nodeLabelPromptCuspUnformed"][:])
            counts += nodeData["nodeLabelPromptCuspUnformed"].shape[0]
    return np.concatenate(labels), counts


# Run the model with pruning disabled, and with it enabled.
labelsUnpruned, countUnpruned = runModel("unpruned", prune=False)
labelsPruned  , countPruned   = runModel("pruned"  , prune=True )

# With pruning disabled, halos in which no cusp can form remain in the tree and so must appear in the output carrying the label.
# If none appear then the model is not probing the regime in which this occurs, and the test below would be vacuous.
countUnformedUnpruned = int(np.count_nonzero(labelsUnpruned != 0))
if countUnformedUnpruned == 0:
    print("FAILED: no halo lacking a prompt cusp arose with pruning disabled - the model does not probe the relevant regime")
    sys.exit(0)

# With pruning enabled, no such halo may survive into the evolved tree.
countUnformedPruned = int(np.count_nonzero(labelsPruned != 0))
if countUnformedPruned > 0:
    print("FAILED: " + str(countUnformedPruned) + " halos lacking a prompt cusp survived pruning")
    sys.exit(0)

print("SUCCESS: halos in which no prompt cusp can form are pruned from merger trees")
print("  halos labeled as unable to form, pruning disabled: " + str(countUnformedUnpruned))
print("  halos in output, pruning disabled/enabled: " + str(countUnpruned) + "/" + str(countPruned))
