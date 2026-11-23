#!/usr/bin/env python3
"""Check the model-scoped fact registry (the `modelFact` and `modelRequire` directives).

Both `galacticDiskHeatingImpulsive` and `darkMatterProfileHeatingTidal` require the fact
`satellite:tidalHeatingNormalized:evolved`, which is established by the `satelliteTidalHeating` node operator. The dark
matter heating object is built with the task, so its requirement is checked when the task is constructed. The disk
heating object is owned by the disk component, so it is built (and records its requirement) only when node components
are initialized, after every node operator already exists, and is checked then. It is checked that:

* with the node operator removed, the model stops with an error that gives the reason for the dark matter heating
  requirement and names the node operator that would satisfy it, before any tree is evolved;
* with the node operator removed and the dark matter heating switched to `null`, the model stops with the
  corresponding error for the disk heating requirement, again before any tree is evolved;
* with the node operator present, the requirements are satisfied and the model runs to completion;
* with the node operator present, a `taskMulti` running two `evolveForests` tasks of the model runs to completion. The facts
  are recorded for the whole model, so they must survive the reset of the parameters at the end of the first task, which
  is followed by the second task's own check after its node components are initialized.

Andrew Benson, Claude
"""
import os
import re
import subprocess
import sys

# This script runs with the working directory set to testSuite/, while the model is run from the repository root.
parameterFileBase = "testSuite/parameters/tidalStrippingGalacticComponentsHeated.xml"
operatorElement   = '<nodeOperator value="satelliteTidalHeating"/>'
outputFileBase    = "testSuite/outputs/tidalStrippingGalacticComponentsHeated.hdf5"
heatingPattern    = re.compile(r'<darkMatterProfileHeating value="tidal">.*?</darkMatterProfileHeating>', re.DOTALL)
fact              = "the fact 'satellite:tidalHeatingNormalized:evolved' is established by"

os.makedirs("outputs", exist_ok=True)

with open(os.path.join("..", parameterFileBase)) as fileBase:
    parametersBase = fileBase.read()
if parametersBase.count(operatorElement) != 1 or len(heatingPattern.findall(parametersBase)) != 1:
    print(f"FAILED: expected exactly one {operatorElement} and one tidal darkMatterProfileHeating in {parameterFileBase}")
    sys.exit(0)

failures = []

def checkMissingProvider(label, reason, nullHeating):
    """Run a copy of the model without the satelliteTidalHeating node operator, and check that it stops with the
    generated error, giving `reason`, before evolving any tree."""
    parameterFileTest = f"testSuite/outputs/modelFactsMissingProvider_{label}.xml"
    parameters        = parametersBase.replace(operatorElement, "")
    if nullHeating:
        parameters = heatingPattern.sub('<darkMatterProfileHeating value="null"/>', parameters)
    parameters = re.sub(r'<outputFileName\s+value="[^"]*"\s*/>', f'<outputFileName value="testSuite/outputs/modelFactsMissingProvider_{label}.hdf5"/>', parameters)
    with open(os.path.join("..", parameterFileTest), "w") as fileTest:
        fileTest.write(parameters)
    status   = subprocess.run(f"cd ..; ./Galacticus.exe {parameterFileTest}", shell=True, capture_output=True)
    log      = status.stdout.decode() + status.stderr.decode()
    failed   = len(failures)
    if status.returncode == 0:
        failures.append(f"{label}: the model without the satelliteTidalHeating node operator ran to completion")
    else:
        for expected, description in (
                (reason         , "the reason for the requirement"),
                (fact           , "the required fact"             ),
                (operatorElement, "the XML selecting the provider"),
        ):
            if expected not in log:
                failures.append(f"{label}: the error report does not give {description}")
        if re.search(r"evolving (?:forest|tree)", log, re.IGNORECASE):
            failures.append(f"{label}: the requirement was not checked before trees were evolved")
    if len(failures) > failed:
        print(log[-3000:])

# The dark matter heating requirement is recorded, and checked, when the task is built.
checkMissingProvider("halo", "tidal heating of dark matter halos requires that the tidal heating of satellites is evolved", nullHeating=False)
# The disk heating requirement is recorded only when node components are initialized, and is checked then.
checkMissingProvider("disk", "heating of disks requires that the tidal heating of satellites is evolved"                  , nullHeating=True )

# With the provider, the requirement is satisfied and the model runs.
status = subprocess.run(f"cd ..; ./Galacticus.exe {parameterFileBase}", shell=True, capture_output=True)
if status.returncode != 0 or not os.path.exists(os.path.join("..", outputFileBase)):
    failures.append("the model with the satelliteTidalHeating node operator did not run")
    print((status.stdout.decode() + status.stderr.decode())[-3000:])

# A taskMulti with two evolveForests tasks. The facts must not be cleared by the parameter reset at the end of the first.
parameterFileMulti = "testSuite/outputs/modelFactsTaskMulti.xml"
parameters         = parametersBase.replace('<task value="evolveForests"/>', '<task value="multi"><task value="evolveForests"/><task value="evolveForests"/></task>')
parameters         = re.sub(r'<treesPerDecade\s+value="[^"]*"\s*/>', '<treesPerDecade value="1"/>', parameters)
parameters         = re.sub(r'<outputFileName\s+value="[^"]*"\s*/>', '<outputFileName value="testSuite/outputs/modelFactsTaskMulti.hdf5"/>', parameters)
if parameters.count('<task value="evolveForests"/>') != 2:
    failures.append("could not construct the taskMulti model")
else:
    with open(os.path.join("..", parameterFileMulti), "w") as fileMulti:
        fileMulti.write(parameters)
    status = subprocess.run(f"cd ..; ./Galacticus.exe {parameterFileMulti}", shell=True, capture_output=True)
    log    = status.stdout.decode() + status.stderr.decode()
    if status.returncode != 0 or log.count("Done task: merger tree evolution") != 2:
        failures.append("the taskMulti model with two evolveForests tasks did not run both tasks")
        print(log[-3000:])

if failures:
    for failure in failures:
        print(f"FAILED: {failure}")
else:
    print("SUCCESS: unmet model fact requirements (checked at task build and at component initialization) are reported, with their provider, before evolution; met ones are accepted, including across the tasks of a taskMulti")
sys.exit(0)
