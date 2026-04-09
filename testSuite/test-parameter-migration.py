#!/usr/bin/env python3
import subprocess
import sys
import xml.etree.ElementTree as ET

# Test migration of parameter files.
# Andrew Benson (ported to Python)

subprocess.run("mkdir -p outputs", shell=True)

# Migrate the test parameter file.
status = subprocess.run(
    "cd ..; ./scripts/aux/parametersMigrate.py testSuite/parameters/parameterMigration.xml testSuite/outputs/parameterMigrated.xml --lastModifiedRevision 6eab8997cd73cb0a474228ade542d133890ad138^",
    shell=True
)
if status.returncode == 0:
    print("PASSED: migration of parameter file")
else:
    print("FAILED: migration of parameter file")

# Parse the migrated parameter file.
tree       = ET.parse("outputs/parameterMigrated.xml")
root       = tree.getroot()

# Check expected state.
nodeOperators = root.findall(".//nodeOperator/nodeOperator")
if nodeOperators:
    firstOperator = nodeOperators[0]
    if firstOperator.find("massDestructionAbsolute") is None and "massDestructionAbsolute" not in firstOperator.attrib:
        print("FAILED: missing parameter 'massDestructionAbsolute'")

if root.find(".//spheroidVerySimpleTrackLuminosities") is not None:
    print("FAILED: unremoved parameter 'spheroidVerySimpleTrackLuminosities'")
