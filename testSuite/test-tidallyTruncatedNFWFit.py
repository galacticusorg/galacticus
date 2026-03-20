#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Test tidally truncated NFW profile fits.
# Andrew Benson (ported to Python)

PI = np.pi

# Run the model.
subprocess.run("mkdir -p outputs", shell=True)
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/tidallyTruncatedNFWFit.xml", shell=True)
if status.returncode != 0:
    print("FAILED: model failed to run")
    sys.exit(0)

# Read model data.
with h5py.File("outputs/tidallyTruncatedNFWFit.hdf5", "r") as model:
    nodes = model["Outputs/Output1/nodeData"]
    data  = {}
    for key in ("basicMass", "darkMatterProfileScale", "darkMatterOnlyRadiusVirial",
                "radiusTidalTruncationNFW", "densityProfile", "densityProfileRadius"):
        data[key] = nodes[key][:]

data["concentration"]        = data["darkMatterOnlyRadiusVirial"] / data["darkMatterProfileScale"]
data["densityNormalization"] = (
    data["basicMass"]
    / data["darkMatterProfileScale"]**3
    / 4.0 / PI
    / (np.log(1.0 + data["concentration"]) - data["concentration"] / (1.0 + data["concentration"]))
)
nHalos = len(data["basicMass"])
data["metricUntruncated"] = np.zeros(nHalos)
data["metricTruncated"  ] = np.zeros(nHalos)

# Iterate over all halos/subhalos.
for i in range(nHalos):
    radii              = data["densityProfileRadius"][i, :]
    densityTarget      = data["densityProfile"      ][i, :]
    xs                 = radii / data["darkMatterProfileScale"  ][i]
    xt                 = radii / data["radiusTidalTruncationNFW"][i]
    densityUntruncated = data["densityNormalization"][i] / xs / (1.0 + xs)**2
    densityTruncated   = data["densityNormalization"][i] / xs / (1.0 + xs)**2 / (1.0 + xt**2)
    data["metricUntruncated"][i] = np.sum(np.log10(densityTarget / densityUntruncated)**2) / len(radii)
    data["metricTruncated"][i]   = np.sum(np.log10(densityTarget / densityTruncated)**2)   / len(radii)

status_str = "succeeded" if np.median(data["metricTruncated"]) < 0.0125 else "FAILED"
print(f"{status_str}: tidally truncated NFW fit")
