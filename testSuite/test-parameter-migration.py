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
    print("SUCCESS: migration of parameter file")
else:
    print("FAILED: migration of parameter file")
    sys.exit(0)

# Parse the migrated parameter file.
tree       = ET.parse("outputs/parameterMigrated.xml")
root       = tree.getroot()

# Check expected state.
failures      = 0
nodeOperators = root.findall(".//nodeOperator/nodeOperator")
if nodeOperators:
    firstOperator = nodeOperators[0]
    if firstOperator.find("massDestructionAbsolute") is None and "massDestructionAbsolute" not in firstOperator.attrib:
        print("FAILED: missing parameter 'massDestructionAbsolute'")
        failures += 1

if root.find(".//spheroidVerySimpleTrackLuminosities") is not None:
    print("FAILED: unremoved parameter 'spheroidVerySimpleTrackLuminosities'")
    failures += 1

# The birth cloud coefficient must have become a column of gas. A coefficient of 0.5 is half the column for which a cloud
# of local interstellar medium metallicity has unit optical depth with the default dust properties, 22.8404 M☉/pc².
birthCloud = root.find(".//dustAttenuation[@value='birthCloud']")
if birthCloud is None:
    print("FAILED: the birthCloud attenuator is missing after migration")
    failures += 1
else:
    if birthCloud.find("coefficient") is not None:
        print("FAILED: unmigrated parameter 'coefficient' remains on birthCloud")
        failures += 1
    column = birthCloud.find("densitySurfaceGas")
    if column is None:
        print("FAILED: missing parameter 'densitySurfaceGas' on birthCloud")
        failures += 1
    elif abs(float(column.get("value")) / (0.5 * 22.8404) - 1.0) > 1.0e-5:
        print(f"FAILED: birthCloud 'densitySurfaceGas' is {column.get('value')}, expected {0.5 * 22.8404:.6e}")
        failures += 1

# The compendium's dust-to-metals ratio must have moved into a dust properties object of its own.
compendium = root.find(".//dustAttenuation[@value='atlasCompendium']")
if compendium is None:
    print("FAILED: the atlasCompendium attenuator is missing after migration")
    failures += 1
else:
    if compendium.find("dustToMetalsRatio") is not None:
        print("FAILED: unmigrated parameter 'dustToMetalsRatio' remains on atlasCompendium")
        failures += 1
    ratio = compendium.find("dustProperties[@value='simple']/dustToMetalsRatio")
    if ratio is None or ratio.get("value") != "0.3":
        print("FAILED: atlasCompendium has no dustProperties object carrying its dust-to-metals ratio")
        failures += 1

if failures == 0:
    print("SUCCESS: migrated parameter file has the expected content")
