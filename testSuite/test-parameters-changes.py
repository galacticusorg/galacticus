#!/usr/bin/env python3
import sys
import subprocess
import lxml.etree as ET

# Check that changes are applied to a parameter file correctly.
# Andrew Benson (26-February-2025)

# Run the model and check for completion.
print("Running model...")
status = subprocess.run("mkdir -p outputs",shell=True)
log = open("outputs/test-parameters-changes.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/changesBase.xml testSuite/parameters/changes.xml; ./scripts/parameters/parametersExtract.py testSuite/outputs/parametersChanged.hdf5 testSuite/outputs/parametersChanged.xml",stdout=log,stderr=log,shell=True)
log.close()
if status.returncode != 0:
    print("...done ("+str(status)+")")
    print("FAILED: model run:")
    subprocess.run("cat outputs/test-parameters-changes.log",shell=True)
    sys.exit()
else:
    print("...done")
    print("Checking for errors...")
    status = subprocess.run("grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" outputs/test-parameters-changes.log",shell=True)
    if status.returncode == 0:
        print("...done ("+str(status)+")")
        print("FAILED: model run (errors):")
        subprocess.run("cat outputs/test-parameters-changes.log",shell=True)
        sys.exit()
    else:
        print("...done")
        print("SUCCESS: model run")

# Check that all changes were applied.
status          = "SUCCESS"
parameters      = ET.parse("outputs/parametersChanged.xml")
galacticFilter  = parameters.findall("./galacticFilter[@value='haloIsolated']")
indexShift      = parameters.findall("./nodeOperator/nodeOperator[@value='indexShift']")
null            = parameters.findall("./nodeOperator/nodeOperator[@value='null']")
blackHolesWinds = parameters.findall("./nodeOperator/nodeOperator[@value='blackHolesWinds']")
feedbackOld     = parameters.findall("./nodeOperator/nodeOperator[@value='stellarFeedbackSpheroids']/stellarFeedbackOutflows/stellarFeedbackOutflows[@value='powerLaw']")
feedbackNew     = parameters.findall("./nodeOperator/nodeOperator[@value='stellarFeedbackSpheroids']/stellarFeedbackOutflows/stellarFeedbackOutflows[@value='vlctyMxSclng']")
stability       = parameters.findall("./nodeOperator/nodeOperator[@value='barInstability']/galacticDynamicsBarInstability/stabilityThresholdGaseous")
lastHost        = parameters.findall("./nodeOperator/nodeOperator[@value='indexLastHost']")
branchTip       = parameters.findall("./nodeOperator/nodeOperator[@value='indexBranchTip']")
efficiency      = parameters.findall("./starFormationRateSpheroids/starFormationTimescale/efficiency")
timescale       = parameters.findall("./starFormationRateSpheroids/starFormationTimescale/timescaleMinimum")
if len(galacticFilter ) != 1:
    print('`<galacticFilter value="haloIsolated">` is not present')
    status="FAILED"
if len(indexShift     ) != 1:
    print('`<nodeOperator value="indexShift">` is not present')
    status="FAILED"
if len(null           ) != 1:
    print('`<nodeOperator value="null">` is not present')
    status="FAILED"
if len(lastHost       ) != 1:
    print('`<nodeOperator value="indexLastHost">` is not present')
    status="FAILED"
if len(branchTip      ) != 1:
    print('`<nodeOperator value="indexBranchTip">` is not present')
    status="FAILED"
if len(blackHolesWinds) != 0:
    print('`<nodeOperator value="blackHolesWinds">` is present')
    status="FAILED"
if len(feedbackOld    ) != 0:
    print('`<stellarFeedbackOutflows value="powerLaw">` is present')
    status="FAILED"
if len(feedbackNew    ) != 1:
    print('`<stellarFeedbackOutflows value="vlctyMxSclng">` is not present')
    status="FAILED"
if len(timescale      ) != 1:
    print('`<timescaleMinimum>` is not present')
    status="FAILED"
if len(stability      ) != 1:
    print('`<stabilityThresholdGaseous>` is not present')
    status="FAILED"
else:
    if stability[0].get('value') != "0.75":
        print('`<stabilityThresholdGaseous>` has incorrect value')
        status="FAILED"
if len(efficiency     ) != 1:
    print('`<efficiency>` is not present')
    status="FAILED"
else:
    if efficiency[0].get('value') != "0.06":
        print('`<efficiency>` has incorrect value')
        status="FAILED"

# Report status.
print(status+": changes to parameters")
