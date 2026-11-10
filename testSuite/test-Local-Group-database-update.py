#!/usr/bin/env python3
"""Test the conversions made by the `localGroupDatabase` task, which fills in missing fields of the Local Group database.

This test writes a small database to a temporary data directory, runs the task on it, and checks the fields which the task
adds against values computed here. It covers two cases which the task has previously got wrong: a declination between 0 and -1
degrees, whose sign is held only by the "degrees" attribute of the sexagesimal form as "-0"; and the uncertainty of a
heliocentric distance computed from a distance modulus.

The task overwrites the database in place, so the real database must never be used here. The path to the temporary
directory is passed to Galacticus explicitly, and the test also confirms that the task wrote to it.

Andrew Benson, Claude
"""

import math
import os
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

toleranceDegrees = 1.0e-9
toleranceDistance = 1.0e-9

database = """<?xml version="1.0" encoding="UTF-8"?>
<localGroup>
  <galaxies>
    <galaxy name="The Galaxy">
      <rightAscensionDecimal value="266.0"/>
      <declinationDecimal value="-29.0"/>
      <distanceHeliocentric value="0.008" uncertainty="0.0001"/>
    </galaxy>
    <galaxy name="Andromeda">
      <rightAscensionDecimal value="10.0"/>
      <declinationDecimal value="41.0"/>
      <distanceHeliocentric value="0.8" uncertainty="0.01"/>
    </galaxy>
    <galaxy name="declinationDecimalSmallNegative">
      <rightAscensionDecimal value="10.0"/>
      <declinationDecimal value="-0.64"/>
      <distanceHeliocentric value="0.1"/>
    </galaxy>
    <galaxy name="declinationDecimalSmallPositive">
      <rightAscensionDecimal value="10.0"/>
      <declinationDecimal value="0.64"/>
      <distanceHeliocentric value="0.1"/>
    </galaxy>
    <galaxy name="declinationDecimalNegative">
      <rightAscensionDecimal value="10.0"/>
      <declinationDecimal value="-31.5"/>
      <distanceHeliocentric value="0.1"/>
    </galaxy>
    <galaxy name="declinationSexagesimalSmallNegative">
      <rightAscensionDecimal value="10.0"/>
      <declinationSexagesimal degrees="-0" arcminutes="32" arcseconds="0"/>
      <distanceHeliocentric value="0.1"/>
    </galaxy>
    <galaxy name="declinationSexagesimalNegative">
      <rightAscensionDecimal value="10.0"/>
      <declinationSexagesimal degrees="-31" arcminutes="30" arcseconds="0"/>
      <distanceHeliocentric value="0.1"/>
    </galaxy>
    <galaxy name="distanceFromModulus">
      <rightAscensionDecimal value="10.0"/>
      <declinationDecimal value="10.0"/>
      <distanceModulus value="20.0" uncertaintyLow="0.1" uncertaintyHigh="0.2"/>
    </galaxy>
    <galaxy name="modulusFromDistance">
      <rightAscensionDecimal value="10.0"/>
      <declinationDecimal value="10.0"/>
      <distanceHeliocentric value="0.1" uncertaintyLow="0.002" uncertaintyHigh="0.004"/>
    </galaxy>
  </galaxies>
</localGroup>
"""

directoryData = tempfile.mkdtemp(prefix="test-Local-Group-database-update-")
directoryDatabase = os.path.join(directoryData, "static/observations/localGroup")
os.makedirs(directoryDatabase)
fileDatabase = os.path.join(directoryDatabase, "localGroupSatellites.xml")
with open(fileDatabase, "w") as output:
    output.write(database)

# Run the task on the temporary database.
os.makedirs("outputs", exist_ok=True)
environment = dict(os.environ, GALACTICUS_DATA_PATH=directoryData + "/")
with open("outputs/test-Local-Group-database-update.log", "w") as log:
    status = subprocess.run(
        ["./Galacticus.exe", "parameters/localGroupDatabaseUpdate.xml"], cwd="..", env=environment, stdout=log, stderr=subprocess.STDOUT
    )
if status.returncode != 0:
    print("FAILED: the localGroupDatabase task did not complete (see outputs/test-Local-Group-database-update.log)")
    shutil.rmtree(directoryData)
    sys.exit(0)

with open(fileDatabase) as updated:
    if updated.read() == database:
        print("FAILED: the localGroupDatabase task did not update the temporary database")
        shutil.rmtree(directoryData)
        sys.exit(0)
galaxies = {galaxy.get("name"): galaxy for galaxy in ET.parse(fileDatabase).getroot().iter("galaxy")}
shutil.rmtree(directoryData)

failures = []


def check(label, value, expected, tolerance):
    """Record a failure if a value differs from that expected."""
    if abs(value - expected) > tolerance:
        failures.append(f"{label}: found {value}, expected {expected}")


def checkSexagesimal(name, negative, arcminutes, arcseconds):
    """Check the sexagesimal declination added by the task, including the sign of its degrees."""
    element = galaxies[name].find("declinationSexagesimal")
    degrees = element.get("degrees")
    if degrees.strip().startswith("-") != negative:
        failures.append(f"{name}: degrees of '{degrees}' has the wrong sign")
    check(f"{name} arcminutes", float(element.get("arcminutes")), arcminutes, toleranceDegrees)
    check(f"{name} arcseconds", float(element.get("arcseconds")), arcseconds, 1.0e-6)


# Declination, from decimal to sexagesimal.
checkSexagesimal("declinationDecimalSmallNegative", True, 38.0, 24.0)
checkSexagesimal("declinationDecimalSmallPositive", False, 38.0, 24.0)
checkSexagesimal("declinationDecimalNegative", True, 30.0, 0.0)

# Declination, from sexagesimal to decimal.
check("declinationSexagesimalSmallNegative", float(galaxies["declinationSexagesimalSmallNegative"].find("declinationDecimal").get("value")), -32.0 / 60.0, toleranceDegrees)
check("declinationSexagesimalNegative", float(galaxies["declinationSexagesimalNegative"].find("declinationDecimal").get("value")), -31.5, toleranceDegrees)

# Heliocentric distance, from distance modulus: d = 10^((mu-25)/5) Mpc, and sigma_d = d ln(10) sigma_mu / 5.
element = galaxies["distanceFromModulus"].find("distanceHeliocentric")
check("distanceFromModulus value", float(element.get("value")), 0.1, toleranceDistance)
check("distanceFromModulus uncertaintyLow", float(element.get("uncertaintyLow")), 0.1 * math.log(10.0) / 5.0 * 0.1, toleranceDistance)
check("distanceFromModulus uncertaintyHigh", float(element.get("uncertaintyHigh")), 0.1 * math.log(10.0) / 5.0 * 0.2, toleranceDistance)

# Distance modulus, from heliocentric distance: sigma_mu = 5 sigma_d / (d ln(10)).
element = galaxies["modulusFromDistance"].find("distanceModulus")
check("modulusFromDistance value", float(element.get("value")), 20.0, toleranceDistance)
check("modulusFromDistance uncertaintyLow", float(element.get("uncertaintyLow")), 5.0 / math.log(10.0) * 0.002 / 0.1, toleranceDistance)
check("modulusFromDistance uncertaintyHigh", float(element.get("uncertaintyHigh")), 5.0 / math.log(10.0) * 0.004 / 0.1, toleranceDistance)

for failure in failures:
    print(f"  -> {failure}")
if failures:
    print("FAILED: Local Group database conversions")
else:
    print("SUCCESS: Local Group database conversions")
sys.exit(0)
