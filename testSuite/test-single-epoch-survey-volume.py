#!/usr/bin/env python3
"""Check that an analysis needing several output times reports a single output time clearly.

`Output_Analysis_Output_Weight_Survey_Volume` refuses a single output time unless the analysis allows it, and aborts with a
message that gives the redshift range in which outputs are needed. To find that range it used to ask for the time at the
survey's maximum distance without passing the analysis' mass limit. For a survey with no upper redshift limit (such as the Li &
White (2009) SDSS geometry used by the galaxySizesSDSS analyses) that distance is the particle horizon, which no time reaches,
so the model aborted instead with a GSL interpolation error ("x values must be strictly increasing") from the comoving distance
table (fixed October 2026).

This script runs quickTest.xml, which has a single output at z=0, with a galaxySizesSDSS analysis, and checks that the model
stops with the intended message.

Andrew Benson, Claude
"""
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

# This script runs with the working directory set to testSuite/, while models are run from the repository root.
pathOutputs = os.path.abspath("outputs/singleEpochSurveyVolume")
os.makedirs(pathOutputs, exist_ok=True)

tree = ET.parse("../parameters/quickTest.xml")
root = tree.getroot()
for name in ("outputFileName", "outputTimes", "mergerTreeOutputter", "outputAnalysis"):
    for element in root.findall(name):
        root.remove(element)
ET.SubElement(ET.SubElement(root, "outputTimes", value="list"), "redshifts", value="0.0")
ET.SubElement(root, "mergerTreeOutputter", value="analyzer")
ET.SubElement(ET.SubElement(root, "outputAnalysis", value="galaxySizesSDSS"), "distributionNumber", value="10")
ET.SubElement(root, "outputFileName", value=os.path.join(pathOutputs, "model.hdf5"))
fileParameters = os.path.join(pathOutputs, "singleEpochSurveyVolume.xml")
tree.write(fileParameters)

status = subprocess.run(["./Galacticus.exe", fileParameters], cwd="..", capture_output=True)
log = status.stdout.decode()+status.stderr.decode()
with open(os.path.join(pathOutputs, "model.log"), "w") as fileLog:
    fileLog.write(log)
if status.returncode == 0:
    print("FAILED: the model ran, but a single output time should have been refused")
elif "single epoch output not permitted" in log and "add at least two outputs in this interval" in log:
    print("SUCCESS: a single output time is refused with the expected message")
else:
    print("FAILED: the model did not stop with the expected message (log in testSuite/outputs/singleEpochSurveyVolume/model.log):")
    print(log[-1500:])
sys.exit(0)
