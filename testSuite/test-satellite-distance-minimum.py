#!/usr/bin/env python3
import subprocess
import sys
import xml.etree.ElementTree as ET
import h5py
import numpy as np

# Check calculations of minimum satellite distance from host center.
# Andrew Benson (ported to Python)

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Run the model.
with open("outputs/satelliteDistanceMinimum.log", "w") as logFile:
    status = subprocess.run(
        "cd ..; ./Galacticus.exe testSuite/parameters/satelliteDistanceMinimum.xml",
        shell=True, stdout=logFile, stderr=subprocess.STDOUT
    )
if status.returncode != 0:
    print("FAILED:  model run:")
    with open("outputs/satelliteDistanceMinimum.log") as f:
        print(f.read())
else:
    print("SUCCESS: model run")

# Read the model data and check for consistency.
with h5py.File("outputs/satelliteDistanceMinimum.hdf5", "r") as model:
    outputs = model["Outputs"]
    for outputName in outputs.keys():
        output          = outputs[outputName]
        nodeData        = output["nodeData"]
        isIsolated      = nodeData["nodeIsIsolated"][:]
        distanceMinimum = nodeData["satelliteDistanceMinimum"][:]
        positionX       = nodeData["satellitePositionX"][:]
        positionY       = nodeData["satellitePositionY"][:]
        positionZ       = nodeData["satellitePositionZ"][:]
        distance        = np.sqrt(positionX**2 + positionY**2 + positionZ**2)
        satellites      = np.where(isIsolated == 0)[0]
        tolerance       = 1.0e-6
        status_str      = "SUCCESS" if np.all(distanceMinimum[satellites] <= distance[satellites] * (1.0 + tolerance)) else "FAILED"
        offset          = distanceMinimum[satellites] / distance[satellites] - 1.0
        print(f"Maximum fractional offset: {offset.max()}")
        print(f"{status_str}: minimum distance from host center")

# Repeat with distances measured relative to the isolated host, and with pre-infall orbits tracked, so that isolated halos
# off the main branch carry positions relative to their own future hosts. For a satellite whose immediate host is isolated,
# the distance to the isolated host is just its own orbital radius, so its minimum distance can not exceed its current
# orbital radius. (Until October 2026 the isolated host's own position was added to the distance.) Outputs at several
# redshifts are used so that isolated hosts off the main branch are present.
tree = ET.parse("parameters/satelliteDistanceMinimum.xml")
root = tree.getroot()
root.find("outputFileName").set("value", "testSuite/outputs/satelliteDistanceMinimumIsolatedHost.hdf5")
root.find("outputTimes").find("redshifts").set("value", "0.5 1.0 1.5 2.0 3.0")
for operator in root.find("nodeOperator").findall("nodeOperator"):
    if operator.get("value") == "satelliteMinimumDistance":
        ET.SubElement(operator, "relativeTo", value="isolatedHost")
    if operator.get("value") == "satelliteOrbit":
        ET.SubElement(operator, "trackPreInfallOrbit", value="true")
tree.write("outputs/satelliteDistanceMinimumIsolatedHost.xml")
with open("outputs/satelliteDistanceMinimumIsolatedHost.log", "w") as logFile:
    status = subprocess.run(
        "cd ..; ./Galacticus.exe testSuite/outputs/satelliteDistanceMinimumIsolatedHost.xml",
        shell=True, stdout=logFile, stderr=subprocess.STDOUT
    )
if status.returncode != 0:
    print("FAILED: model run (isolated host)")
    sys.exit(0)
print("SUCCESS: model run (isolated host)")
countChecked = 0
countFailed  = 0
with h5py.File("outputs/satelliteDistanceMinimumIsolatedHost.hdf5", "r") as model:
    for outputName, output in model["Outputs"].items():
        nodeData        = output["nodeData"]
        isIsolated      = nodeData["nodeIsIsolated"][:]
        nodeIndex       = nodeData["nodeIndex"     ][:]
        parentIndex     = nodeData["parentIndex"   ][:]
        distanceMinimum = nodeData["satelliteDistanceMinimum"][:]
        distance        = np.sqrt(nodeData["satellitePositionX"][:]**2+nodeData["satellitePositionY"][:]**2+nodeData["satellitePositionZ"][:]**2)
        isolatedIndices = set(nodeIndex[isIsolated == 1])
        firstLevel      = np.array([isIsolated[i] == 0 and parentIndex[i] in isolatedIndices for i in range(nodeIndex.size)], dtype=bool)
        countChecked   += np.count_nonzero(firstLevel)
        countFailed    += np.count_nonzero(distanceMinimum[firstLevel] > distance[firstLevel]*(1.0+1.0e-6))
if countChecked == 0:
    print("FAILED: no satellites of isolated hosts were found")
elif countFailed > 0:
    print(f"FAILED: minimum distance from the isolated host exceeds the current orbital radius for {countFailed} of {countChecked} satellites of isolated hosts")
else:
    print(f"SUCCESS: minimum distance from the isolated host (checked for {countChecked} satellites of isolated hosts)")
