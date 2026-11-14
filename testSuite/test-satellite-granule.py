#!/usr/bin/env python3
import sys
import os
import subprocess
import xml.etree.ElementTree as ET
import h5py
import numpy as np

# Test the satellite granule operator for soliton+NFW, soliton-only, and NFW-only host halos.
# Yu Zhao (10-November-2026)
#
# The model is run three times: with both the core and outer granule contributions, with only the outer contribution
# (alphaCore=0), and with neither (alphaCore=alphaOuter=0). Switching a contribution off by setting its amplitude to zero, rather
# than by removing the operator, keeps the sequence of random numbers drawn unchanged, so that the host halos are identical in all
# three runs. Comparing the satellite velocities between runs then shows which contributions acted on each satellite. The test
# also checks that each host halo is in its intended state, since that state is determined by the dark matter profile, not by the
# merger tree file.
#
# The comparison is made over the interval between the two outputs, by which time every host has reached its final state. A
# satellite is evolved through its first timestep before its host is, so a satellite of the soliton-only host feels the outer
# contribution briefly at the start, while that host is still in its initial soliton+NFW state. Comparing the change in the
# velocity difference over the interval excludes that initial kick.

# Soliton status values (see source/dark_matter_profiles_DMO/soliton_status.F90).
statusSolitonNFW  = 1
statusSolitonOnly = 2
statusNFWOnly     = 3

# Host halos, with their intended states, and the satellites of each one. For each satellite, "core" and "outer" give whether the
# core and outer granule contributions should act on it.
hosts = [
    {"index": 101, "status": statusSolitonNFW , "label": "soliton+NFW" , "satellites": {
        103: {"region": "inside the core radius"            , "core": True , "outer": True },
        104: {"region": "between core and soliton radii"    , "core": False, "outer": True },
        105: {"region": "outside the soliton radius"        , "core": False, "outer": True }}},
    {"index": 203, "status": statusSolitonOnly, "label": "soliton-only", "satellites": {
        204: {"region": "inside the core radius"            , "core": True , "outer": False}}},
    {"index": 301, "status": statusNFWOnly    , "label": "NFW-only"    , "satellites": {
        303: {"region": "within the halo"                   , "core": False, "outer": True }}},
]

# Fractional changes in velocity, accumulated over the interval between the two outputs, below which a contribution is
# considered to be absent, and above which it is considered to be present. Typical changes are 10^-3 or more when a contribution
# is present, and either exactly zero or about 10^-10 when it is absent.
changeAbsent  = 1.0e-8
changePresent = 1.0e-5

os.makedirs("outputs", exist_ok=True)

# Construct parameter files for the three runs.
runs = {
    "full"     : {"alphaCore": None , "alphaOuter": None },
    "outerOnly": {"alphaCore": "0.0", "alphaOuter": None },
    "none"     : {"alphaCore": "0.0", "alphaOuter": "0.0"},
}
for name, alphas in runs.items():
    tree = ET.parse("parameters/satelliteGranule.xml")
    root = tree.getroot()
    for operator in root.find("nodeOperator").findall("nodeOperator"):
        if operator.get("value") == "satelliteGranule":
            for parameter, value in alphas.items():
                if value is not None:
                    operator.find(parameter).set("value", value)
    root.find("outputFileName").set("value", "testSuite/outputs/satelliteGranule_" + name + ".hdf5")
    tree.write("outputs/satelliteGranule_" + name + ".xml")

# Run the models.
for name in runs.keys():
    outputFile = "outputs/satelliteGranule_" + name + ".hdf5"
    if os.path.exists(outputFile):
        os.remove(outputFile)
    status = subprocess.run("cd ..; ./Galacticus.exe testSuite/outputs/satelliteGranule_" + name + ".xml", shell=True)
    if status.returncode != 0:
        print("FAILED: satellite granule model '" + name + "' failed to run")
        sys.exit(0)
    if not os.path.exists(outputFile):
        print("FAILED: satellite granule model '" + name + "' produced no output file")
        sys.exit(0)

def readOutputs(fileName):
    # Read the soliton status and satellite velocity of each node in each output, keyed by output name and then node index.
    data = {}
    with h5py.File(fileName, "r") as model:
        for outputName, output in model["Outputs"].items():
            nodeData = output["nodeData"]
            velocity = np.stack([nodeData["satelliteVelocity" + axis][:] for axis in ("X", "Y", "Z")], axis=1)
            data[outputName] = {
                int(index): {"status": int(status), "velocity": velocity[i]}
                for i, (index, status) in enumerate(zip(nodeData["nodeIndex"][:], nodeData["solitonStatus"][:]))
            }
    return data

results = {name: readOutputs("outputs/satelliteGranule_" + name + ".hdf5") for name in runs.keys()}

outputNames = sorted(results["full"].keys(), key=lambda name: int(name[len("Output"):]))
if len(outputNames) != 2:
    print(f"FAILED: expected two outputs, found {len(outputNames)}")
    sys.exit(0)

def checkContribution(index, label, contribution, runWith, runWithout, expected):
    # Check whether a contribution acted on a satellite during the interval between the two outputs, from the change over that
    # interval in the difference in its velocity between runs with and without the contribution.
    difference = [results[runWith][outputName][index]["velocity"]-results[runWithout][outputName][index]["velocity"] for outputName in outputNames]
    change     = np.linalg.norm(difference[1]-difference[0])/np.linalg.norm(results[runWithout][outputNames[1]][index]["velocity"])
    print(f"satellite {index} ({label}): {contribution} contribution: fractional change in velocity = {change:.3e}")
    if expected and change < changePresent:
        print(f"FAILED: the {contribution} granule contribution did not act on satellite {index} ({label})")
        return False
    if not expected and change > changeAbsent:
        print(f"FAILED: the {contribution} granule contribution acted on satellite {index} ({label}), but should not")
        return False
    return True

success = True
for host in hosts:
    # Check that the host is in the intended state at both outputs.
    for outputName in outputNames:
        if any(host["index"] not in results[name][outputName] for name in runs.keys()):
            print(f"FAILED: {host['label']} host {host['index']} is missing from {outputName}")
            success = False
            continue
        statusHost = results["full"][outputName][host["index"]]["status"]
        if statusHost != host["status"]:
            print(f"FAILED: host {host['index']} has soliton status {statusHost} in {outputName}, but should be {host['label']} ({host['status']})")
            success = False
    for index, satellite in host["satellites"].items():
        label = f"{satellite['region']} of the {host['label']} host"
        if any(index not in results[name][outputName] for name in runs.keys() for outputName in outputNames):
            print(f"FAILED: satellite {index} ({label}) is missing from an output")
            success = False
            continue
        if not all(np.all(np.isfinite(results[name][outputName][index]["velocity"])) for name in runs.keys() for outputName in outputNames):
            print(f"FAILED: satellite {index} ({label}) has a non-finite velocity")
            success = False
            continue
        # The core contribution is the difference between the full and outer-only runs; the outer contribution is the difference
        # between the outer-only run and the run with neither contribution.
        success &= checkContribution(index, label, "core" , "full"     , "outerOnly", satellite["core" ])
        success &= checkContribution(index, label, "outer", "outerOnly", "none"     , satellite["outer"])

if success:
    print("SUCCESS: the core and outer granule contributions acted as expected on satellites of soliton+NFW, soliton-only, and NFW-only hosts")
