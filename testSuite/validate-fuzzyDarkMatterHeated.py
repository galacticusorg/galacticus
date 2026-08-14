#!/usr/bin/env python3
import sys
import os
import subprocess
import numpy as np
import h5py

# Run a model to validate the solitonNFWHeated fuzzy dark matter model, and check the structural state of the resulting halos.
# Yu Zhao (10-August-2026)

# Create output path.
try:
    os.mkdir("outputs")
except FileExistsError:
    pass

fileName = "outputs/validate_darkMatterOnlySubhalos_Symphony_MilkyWay_resolutionX1_CDM.hdf5"

# Run the validation model. The final parameter file additionally writes per-node data, so that the structural state of the fuzzy
# dark matter halos can be checked below - the output analyses alone do not expose it.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/validate_darkMatterOnlySubhalos_Symphony_resolutionX1_CDM.xml parameters/reference/powerSpectraSuppressed.xml parameters/reference/fuzzyDarkMatterHeated.xml testSuite/parameters/resolutionM1e9.xml testSuite/parameters/fuzzyDarkMatterHeatedValidate.xml",shell=True)
if status.returncode != 0:
    print("FAILED: Fuzzy Dark Matter validation model using solitonNFWHeated failed to run")
    sys.exit(0)
if not os.path.exists(fileName):
    print("FAILED: Fuzzy Dark Matter validation model using solitonNFWHeated produced no output file")
    sys.exit(0)

# Members of the `solitonStatus` enumeration (see source/dark_matter_profiles_DMO/soliton_status.F90). `uninitialized` is
# necessarily zero, as meta-properties read as zero before they are first set.
statusUninitialized = 0
statusSolitonNFW    = 1
statusSolitonOnly   = 2
statusNFWOnly       = 3

failures = []
model    = h5py.File(fileName,"r")

# Check that the output analyses were evaluated.
for analysisName in ( "subhaloMassFunction", "subhaloRadialDistribution", "subhaloVelocityMaximumMean" ):
    if analysisName not in model['analyses']:
        failures.append("output analysis '"+analysisName+"' is missing")

nodeData      = model['Outputs/Output1/nodeData']
solitonStatus = np.array(nodeData['solitonStatus'   ]).astype(int)
radiusSoliton = np.array(nodeData['solitonRadiusSoliton'])
massCore      = np.array(nodeData['solitonMassCore' ])
countNodes    = len(solitonStatus)

if countNodes == 0:
    failures.append("model contains no nodes")
else:
    # Every halo must be in a state which is a member of the enumeration.
    statusesValid = ( statusUninitialized, statusSolitonNFW, statusSolitonOnly, statusNFWOnly )
    statusesFound = set(solitonStatus.tolist())
    statusesOther = statusesFound-set(statusesValid)
    if statusesOther:
        failures.append("halos found in unrecognized soliton states: "+str(sorted(statusesOther)))

    # Solitons must form. If none do, the core-halo mass relation or the solution for the soliton transition radius is broken.
    countSolitonNFW = np.count_nonzero(solitonStatus == statusSolitonNFW)
    if countSolitonNFW == 0:
        failures.append("no halo has a solitonic core")

    # Tidal stripping must be able to remove the NFW envelope from a satellite and leave it soliton-only. If none do, the tidal
    # radius is not being limited to the soliton radius, or the transition to the soliton-only state is never triggered.
    countSolitonOnly = np.count_nonzero(solitonStatus == statusSolitonOnly)
    if countSolitonOnly == 0:
        failures.append("no halo reached the soliton-only state")

    # A halo with a solitonic core embedded in an NFW envelope must have a transition radius between the two.
    selection = solitonStatus == statusSolitonNFW
    if np.count_nonzero(selection) > 0 and np.any(radiusSoliton[selection] <= 0.0):
        failures.append("soliton+NFW halos found with a non-positive soliton transition radius")

    # A halo stripped down to its core has no NFW envelope, and so no transition radius.
    selection = solitonStatus == statusSolitonOnly
    if np.count_nonzero(selection) > 0 and np.any(radiusSoliton[selection] > 0.0):
        failures.append("soliton-only halos found with a positive soliton transition radius")

    # Any halo with a solitonic core must have a positive core mass.
    selection = (solitonStatus == statusSolitonNFW) | (solitonStatus == statusSolitonOnly)
    if np.count_nonzero(selection) > 0 and np.any(massCore[selection] <= 0.0):
        failures.append("halos with a solitonic core found with a non-positive core mass")

    print("Halo states: %d nodes; %d soliton+NFW, %d soliton-only, %d NFW-only, %d uninitialized" %
          (
              countNodes                                                    ,
              countSolitonNFW                                               ,
              countSolitonOnly                                              ,
              np.count_nonzero(solitonStatus == statusNFWOnly      )         ,
              np.count_nonzero(solitonStatus == statusUninitialized)
          ))

model.close()

if failures:
    for failure in failures:
        print("FAILED: Fuzzy Dark Matter solitonNFWHeated validation: "+failure)
    sys.exit(0)

print("SUCCESS: Fuzzy Dark Matter validation model using solitonNFWHeated")
