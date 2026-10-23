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

# The threshold was raised from 0.0125 to 0.015 when the tidal heating specific energy gradient was corrected to include
# the radial derivative of its second-order coefficient (which depends on the density logarithmic slope). That correction
# makes the heated density consistent with the Jacobian of the profile's own initial-to-final radius mapping - previously
# the class' density did not equal (1/4 pi r^2) dM/dr - and the analytic truncated-NFW form fits the corrected profile
# slightly less well: the median metric moved from 0.011909 to 0.013424, measured by zeroing only that term. In RMS terms
# the residual grew from 0.109 to 0.116 dex. The threshold measures how well the fitting form describes the profile, not
# whether the profile is right, so it is the threshold which moves.
#
# The threshold was then raised from 0.015 to 0.0165 because the median depends on the tree realization: over eight seeds it
# averages 0.0126 (range 0.0111-0.0139) with the Cole2000 builder as it was, and 0.0130 (range 0.0114-0.0150) after the builder
# stopped removing sub-resolution accretion twice from the secondary progenitor - no significant shift, but this test's seed
# gave 0.01504, the highest of all sixteen runs. 0.0165 is about three times the seed-to-seed scatter above the mean.
status_str = "SUCCESS" if np.median(data["metricTruncated"]) < 0.0165 else "FAILED"
print(f"{status_str}: tidally truncated NFW fit")
