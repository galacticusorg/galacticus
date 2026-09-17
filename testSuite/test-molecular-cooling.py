#!/usr/bin/env python3
import subprocess
import sys
import h5py
import re
import os
import numpy as np
import warnings

# Run a simple molecular cooling model and verify results against prior known values.
# Andrew Benson (04-April-2025)

# Run the model and check for completion.
print("   Running model...")
status = subprocess.run("mkdir -p outputs",shell=True)
log = open("outputs/test-molecular-cooling.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/molecularCooling.xml",stdout=log,stderr=log,shell=True)
log.close()
if status.returncode != 0:
    print("   ...done ("+str(status)+")")
    print("   FAILED: model run:")
    subprocess.run("cat outputs/test-molecular-cooling.log",shell=True)
    sys.exit(0)
else:
    print("   ...done")
    print("   Checking for errors...")
    status = subprocess.run("grep -q -i -e fatal -e aborted -e \"task failed\" -e \"Galacticus experienced an error in the GSL library\" outputs/test-molecular-cooling.log",shell=True)
    if status.returncode == 0:
        print("   ...done ("+str(status)+")")
        print("   FAILED: model run (errors):")
        subprocess.run("cat outputs/test-molecular-cooling.log",shell=True)
        sys.exit(0)
    else:
        print("   ...done")
        print("   SUCCESS: model run")

# Open the model and extract the molecular masses and cooling functions.
model                  = h5py.File('outputs/molecularCooling.hdf5','r')
outputs                = model['Outputs'                                                 ]
# Iterate over outputs.
massesHydrogenMolecular   = np.zeros(4)
coolingFunctionsMolecular = np.zeros(4)
for (outputName, output) in outputs.items():
    match = re.match(r'^Output(\d+)',outputName) 
    if not match:
        continue
    outputIndex                            = int(match.group(1))-1
    nodes                                  = output                  ['nodeData'                         ]
    massHydrogenMolecular                  = nodes                   ['hotHaloChemicalsMolecularHydrogen'][:]
    coolingFunctionMolecular               = nodes                   ['cgmMolecularCoolingFunction'      ][:]
    massesHydrogenMolecular  [outputIndex] = massHydrogenMolecular                                        [0]
    coolingFunctionsMolecular[outputIndex] = coolingFunctionMolecular                                     [0][0]

# Target values. Outputs are ordered by increasing cosmic time, so these correspond to z = 12, 10, 8, and 6 respectively.
#
# The z=8 values depend on how much molecular hydrogen is destroyed between z=10 and z=8, and so on how the photodissociation
# rate coefficient rises as the ultraviolet background switches on. That switch-on is a step: the tabulated background used here
# extends only to z=10 (a cosmic time of 0.47226 Gyr), below which it is extrapolated to zero. Any tabulation of the rate
# coefficient smears that step over one of its intervals, and formation and destruction of molecular hydrogen are finely
# balanced in the interval which follows it - the molecular hydrogen mass rises from z=10 to z=8 here, where a slightly earlier
# switch-on would have it fall. These two values are therefore sensitive to where the tabulation's points fall relative to the
# step, and were updated when those points were placed on an absolute lattice. The remaining outputs are insensitive to it.
#
# The values were updated again when the Blitz & Rosolowsky (2006) star formation rate surface density law, which this model
# uses, was corrected: its characteristic pressure had been a factor of some 7,700 too small, and its molecular fraction was
# min(R_mol,1) rather than R_mol/(1+R_mol). That model's disk is consequently far less molecular, its star formation rate lower,
# and its feedback weaker, which feeds through to the hot halo. The masses fell by 2.0%, 6.4% and 8.4% at z = 12, 10 and 8, and
# the cooling functions by 0.3%, 3.0% and 7.0%. Restoring the two superseded expressions - with the tabulation of the disk
# integral still removed - reproduces the previous values exactly, which is what identifies the star formation law rather than
# the removal of that tabulation as the cause.
#
# Repeated builds of the same source give these to about one part in 10^5, so the 1% tolerance below is not tight.
massesHydrogenMolecularTarget   = np.array([94177.45499300, 221804.24175394, 336649.62031359, 0.            ])
coolingFunctionsMolecularTarget = np.array([2.75853507e-24, 5.34553994e-24 , 2.32584862e-24 , 0.00000000e+00])

# Report on status.
massesAgree  = np.isclose(massesHydrogenMolecular  ,massesHydrogenMolecularTarget  ,rtol=1.0e-2,atol=1.0e3 )
coolingAgree = np.isclose(coolingFunctionsMolecular,coolingFunctionsMolecularTarget,rtol=1.0e-2,atol=1.0e-26)
status       = np.all(massesAgree) and np.all(coolingAgree)
if status:
    print("   SUCCESS: H₂ mass and cooling function")
else:
    print("   FAILED: H₂ mass and cooling function")
    # Report which outputs disagree, and by how much, so that a failure says whether the model has moved a little or a lot.
    redshifts = ( 12, 10, 8, 6 )
    for i in range(4):
        if not massesAgree[i]:
            print("     z = {0:2d}: H₂ mass          = {1:18.8f}, expected {2:18.8f} ({3:+7.2f}%)".format(redshifts[i],massesHydrogenMolecular[i],massesHydrogenMolecularTarget[i],100.0*(massesHydrogenMolecular[i]/massesHydrogenMolecularTarget[i]-1.0) if massesHydrogenMolecularTarget[i] != 0.0 else np.nan))
        if not coolingAgree[i]:
            print("     z = {0:2d}: cooling function = {1:18.8e}, expected {2:18.8e} ({3:+7.2f}%)".format(redshifts[i],coolingFunctionsMolecular[i],coolingFunctionsMolecularTarget[i],100.0*(coolingFunctionsMolecular[i]/coolingFunctionsMolecularTarget[i]-1.0) if coolingFunctionsMolecularTarget[i] != 0.0 else np.nan))
