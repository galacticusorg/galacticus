#!/usr/bin/env python3
import subprocess
import os
import sys
import h5py
import numpy as np
import xml.etree.ElementTree as ET

# Test spheroidRadiusPowerLaw nodeOperator.
# Charles Gannon (31-August-2024)
path_in_param = os.path.abspath(f"parameters/test-spheroid-power-law.xml"                       )
path_out_dir  = os.path.abspath(f"outputs/test-spheroid-power-law"                              )
path_out_log  = os.path.abspath(f"outputs/test-spheroid-power-law/test-spheroid-power-law.log"  )
path_out_hdf  = os.path.abspath(f"outputs/test-spheroid-power-law/test-spheroid-power-law.hdf5" )

status = subprocess.run(f"mkdir -p {path_out_dir}",shell=True)

# Run the model and check for completion 
print("Running model...")
 
log = open(path_out_log,"w")
status = subprocess.run(f"cd ..; ./Galacticus.exe {path_in_param}",stdout=log,stderr=log,shell=True)
log.close()

print("...done ("+str(status)+")")

if status.returncode != 0:
    print("FAILED: model run:")
    subprocess.run(f"cat {path_out_log}",shell=True)
    sys.exit()
print("Checking for errors...")
status = subprocess.run(f"grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" {path_out_log}",shell=True)
print("...done ("+str(status)+")")
if status.returncode == 0:
    print("FAILED: model run (errors):")
    subprocess.run(f"cat {path_out_log}",shell=True)
    sys.exit()
print("SUCCESS: model run")
 
# Open the model and extract the recycled fraction.
model         = h5py.File(path_out_hdf,'r')
massStellar   = model["Outputs/Output1/nodeData/spheroidMassStellar" ][:]
radiusStellar = model["Outputs/Output1/nodeData/spheroidRadius"      ] 
isIsolated    = model["Outputs/Output1/nodeData/nodeIsIsolated"      ][:].astype(bool)
 
massStellarHost, radiusStellarHost =  massStellar[isIsolated], radiusStellar[isIsolated]
 
alpha = 0.56
beta  = 1.19E-9

radiusStellarHostPython = (
                            +beta
                            *massStellarHost**alpha
                          )  
if np.allclose(np.log10(radiusStellarHost), np.log10(radiusStellarHostPython), rtol=1e-6):
    print(f"SUCCESS: results do agree"   )
else:
    print(f"FAILED: results do not agree")
