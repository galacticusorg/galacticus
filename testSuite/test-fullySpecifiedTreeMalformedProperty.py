#!/usr/bin/env python3
"""Check that a malformed fully-specified merger tree is reported, not silently ignored.

Two ways of getting such a tree wrong used to leave the run to do nothing while reporting nothing (issue #1512):

  * a rank-1 node component property is read one XML element per value, and writing it as a single element holding several
    values left the run to be stopped by the XML layer with a zero exit status and a two-line message on standard error naming
    neither the file, the node, nor the property;
  * a node which is listed in the file but is not the child, sibling, or satellite of any other node is simply absent from the
    constructed tree, which was left quietly incomplete.

Andrew Benson
"""

import os
import subprocess
import sys

treeFileName     ="parameters/fullySpecifiedTreeMalformedPropertyTree.xml"
parameterFileName="parameters/fullySpecifiedTreeMalformedProperty.xml"

# Ensure output directory exists.
os.makedirs("outputs",exist_ok=True)

def runModel(name,parameterFile):
    """Run the model and return its exit status and log."""
    logFileName="outputs/"+name+".log"
    with open(logFileName,"w") as logFile:
        status = subprocess.run(
            "cd ..; export OMP_NUM_THREADS=1; ./Galacticus.exe testSuite/"+parameterFile,
            shell=True, stdout=logFile, stderr=subprocess.STDOUT
        )
    return status.returncode, open(logFileName).read()

def check(description,returnCode,log,phrases):
    """Require that the run failed, with a message naming what went wrong and where."""
    if returnCode == 0:
        print("FAILED: "+description+" was not detected")
        return
    missing=[phrase for phrase in phrases if phrase not in log]
    if missing:
        print("FAILED: error message for "+description+" omits: "+", ".join(missing))
    else:
        print("SUCCESS: "+description+" was detected and reported")

# A rank-1 property written as a single element holding three values.
returnCode, log = runModel("fullySpecifiedTreeMalformedProperty",parameterFileName)
check(
    "a malformed rank-1 property in a fully-specified merger tree",
    returnCode, log,
    ["'position'","'orbiting'","'satellite'","node 2","too many values found"]
)

# A node which is specified in the file but is linked to by no other node. Build this case from the files above, so that the
# two differ only in the `firstSatellite` element which is removed here, and in the `position` property which is repaired.
tree=open(treeFileName).read()
tree=tree.replace(
    "      <position>+0.310934597000000 +0.000000000000000 +0.000000000000000</position>\n",
    "      <position>+0.310934597000000</position>\n"
    "      <position>+0.000000000000000</position>\n"
    "      <position>+0.000000000000000</position>\n"
)
tree=tree.replace("    <firstSatellite>2</firstSatellite>\n","")
open("outputs/fullySpecifiedTreeUnreachableNodeTree.xml","w").write(tree)
parameters=open(parameterFileName).read().replace(
    "testSuite/"+treeFileName,
    "testSuite/outputs/fullySpecifiedTreeUnreachableNodeTree.xml"
)
open("outputs/fullySpecifiedTreeUnreachableNode.xml","w").write(parameters)
returnCode, log = runModel("fullySpecifiedTreeUnreachableNode","outputs/fullySpecifiedTreeUnreachableNode.xml")
check(
    "a node unreachable from the root of a fully-specified merger tree",
    returnCode, log,
    ["node 2","not reachable"]
)

sys.exit(0)
