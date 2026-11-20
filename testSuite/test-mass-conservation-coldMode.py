#!/usr/bin/env python3
import subprocess
import sys
import glob
import h5py
import numpy as np
import xml.etree.ElementTree as ET

# Run a cold mode model to test mass conservation.
# Andrew Benson (ported to Python); CGM starvation variant added by Claude (19-November-2026)

# Run the model.
subprocess.run(
    "cd ..; ./scripts/aux/launch.py testSuite/parameters/test-mass-conservation-coldMode.xml --launchMethod local --threadMaximum 1 --ompThreads 4",
    shell=True
)

# Check for failed models.
logFiles = glob.glob("outputs/test-mass-conservation-coldMode/galacticus_*/galacticus.log")
failures = []
for logFile in logFiles:
    result = subprocess.run(f"grep -q -i -e fatal -e aborted {logFile}", shell=True)
    if result.returncode == 0:
        failures.append(logFile)

if failures:
    for failure in failures:
        print(f"FAILED: log from {failure}:")
        with open(failure) as f:
            print(f.read())
else:
    print("SUCCESS: model run")

# Extract masses and check they add up to the expected value.
def checkConservation(fileName,label):
  with h5py.File(fileName, "r") as galacticus:
      cosmologyParams = galacticus["Parameters/cosmologyParameters"]
      OmegaBaryon     = cosmologyParams.attrs["OmegaBaryon"]
      OmegaMatter     = cosmologyParams.attrs["OmegaMatter"]
      baryonFraction  = OmegaBaryon / OmegaMatter

      nodeData = galacticus["Outputs/Output1/nodeData"]
      props    = {}
      for key in ("mergerTreeWeight", "blackHoleMass", "diskMassStellar", "diskMassGas",
                  "spheroidMassStellar", "spheroidMassGas",
                  "hotHaloMass", "hotHaloMassCold", "hotHaloOutflowedMass",
                  "hotHaloUnaccretedMass", "nodeIsIsolated", "basicMass", "mergerTreeIndex"):
          props[key] = nodeData[key][:]

  centrals   = np.where(props["nodeIsIsolated"] == 1)[0]
  satellites = np.where(props["nodeIsIsolated"] == 0)[0]

  massSatellites = np.zeros(len(centrals))
  for i, ci in enumerate(centrals):
      treeIdx = props["mergerTreeIndex"][ci]
      inTree  = np.where(props["mergerTreeIndex"][satellites] == treeIdx)[0]
      satMass = (
          props["diskMassStellar"][satellites][inTree].sum()
          + props["diskMassGas"][satellites][inTree].sum()
          + props["spheroidMassStellar"][satellites][inTree].sum()
          + props["spheroidMassGas"][satellites][inTree].sum()
          + props["hotHaloMass"][satellites][inTree].sum()
          + props["hotHaloMassCold"][satellites][inTree].sum()
          + props["hotHaloOutflowedMass"][satellites][inTree].sum()
          + props["hotHaloUnaccretedMass"][satellites][inTree].sum()
          + props["blackHoleMass"][satellites][inTree].sum()
      )
      massSatellites[i] = satMass / baryonFraction / props["basicMass"][ci]

  normFactor = baryonFraction * props["basicMass"]
  massTotal = (
      props["diskMassStellar"][centrals]         / normFactor[centrals]
      + props["diskMassGas"][centrals]           / normFactor[centrals]
      + props["spheroidMassStellar"][centrals]   / normFactor[centrals]
      + props["spheroidMassGas"][centrals]       / normFactor[centrals]
      + props["hotHaloMass"][centrals]           / normFactor[centrals]
      + props["hotHaloMassCold"][centrals]       / normFactor[centrals]
      + props["hotHaloOutflowedMass"][centrals]  / normFactor[centrals]
      + props["hotHaloUnaccretedMass"][centrals] / normFactor[centrals]
      + props["blackHoleMass"][centrals]         / normFactor[centrals]
      + massSatellites
  )

  if np.any(np.abs(massTotal - 1.0) > 1.0e-3):
      print(f"FAILED: mass conservation [{label}] failure -> {massTotal - 1.0}")
  else:
      print(f"SUCCESS: mass conservation [{label}]")

checkConservation("outputs/test-mass-conservation-coldMode/galacticus_0:1/galacticus.hdf5","coldMode")

# Repeat with the CGMStarvation node operator, set to trim the CGM to the universal baryon fraction at node mergers. This removes
# hot and cold mode CGM in proportion, moving it to the unaccreted reservoir, so mass must still be conserved. Starvation is
# limited to outflowed gas, which is then moved to the host's CGM (rather than being removed from the model, as is all CGM of
# satellites otherwise).
tree       = ET.parse("parameters/test-mass-conservation-coldMode.xml")
grid       = tree.getroot()
grid.find("modelRootDirectory").text = "testSuite/outputs/test-mass-conservation-coldMode-starvation"
parameters = grid.find("parameters")
operators  = [operator for operator in parameters.findall("nodeOperator") if operator.get("value") == "multi"][0]
starvation = ET.SubElement(operators ,"nodeOperator"                   ,{"value": "CGMStarvation", "iterable": "no"})
ET.SubElement(             starvation,"fractionBaryonLimitInNodeMerger",{"value": "true"         })
ET.SubElement(             starvation,"starveOutflowsOnly"             ,{"value": "true"         })
tree.write("outputs/test-mass-conservation-coldMode-starvation.xml")
subprocess.run(
    "cd ..; ./scripts/aux/launch.py testSuite/outputs/test-mass-conservation-coldMode-starvation.xml --launchMethod local --threadMaximum 1 --ompThreads 4",
    shell=True
)
logFiles = glob.glob("outputs/test-mass-conservation-coldMode-starvation/galacticus_*/galacticus.log")
failures = [logFile for logFile in logFiles if subprocess.run(f"grep -q -i -e fatal -e aborted {logFile}",shell=True).returncode == 0]
if not logFiles or failures:
    print("FAILED: model with CGM starvation failed")
    for failure in failures:
        with open(failure) as f:
            print(f.read())
else:
    print("SUCCESS: model with CGM starvation run")
    checkConservation("outputs/test-mass-conservation-coldMode-starvation/galacticus_0:1/galacticus.hdf5","coldMode with CGM starvation")
