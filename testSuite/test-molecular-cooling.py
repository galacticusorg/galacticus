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
    sys.exit()
else:
    print("   ...done")
    print("   Checking for errors...")
    status = subprocess.run("grep -q -i -e fatal -e aborted -e \"task failed\" -e \"Galacticus experienced an error in the GSL library\" outputs/test-molecular-cooling.log",shell=True)
    if status.returncode == 0:
        print("   ...done ("+str(status)+")")
        print("   FAILED: model run (errors):")
        subprocess.run("cat outputs/test-molecular-cooling.log",shell=True)
        sys.exit()
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

# Target values.
massesHydrogenMolecularTarget   = np.array([96105.91566116, 237085.66451119, 119155.89135269, 0.            ])
coolingFunctionsMolecularTarget = np.array([2.76678102e-24, 5.50976858e-24 , 7.04508046e-25 , 0.00000000e+00])

# Report on status.
status = np.allclose(massesHydrogenMolecular,massesHydrogenMolecularTarget,rtol=1.0e-2,atol=1.0e3) and np.allclose(coolingFunctionsMolecular,coolingFunctionsMolecularTarget,rtol=1.0e-2,atol=1.0e-26)
if status:
    print("   SUCCESS: H₂ mass and cooling function")
else:
    print("   FAILED: H₂ mass and cooling function")
