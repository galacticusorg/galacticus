#!/usr/bin/env python3
import subprocess
import sys
import h5py
import math
import numpy as np

# Test profiles of decaying dark matter models.
# Andrew Benson (28-October-2024)

# Set target densities for these models. Each model is a single built tree, so these values depend on its realization: they
# were updated when the Cole2000 builder stopped removing sub-resolution accretion twice from the secondary progenitor, which
# changes the tree and so the halo's concentration (by 1.7-4.9%).
densitiesTarget = {
    "CDM":           3.81736478e+16,
    "massLossFalse": 2.74626651e+15,
    "gamma0.0":      2.58908136e+15,
    "gamma0.5":      1.67163002e+15
    }

# Run the model and check for completion
print("Running models...")
status = subprocess.run("mkdir -p outputs",shell=True)
densities = {}
failed = False
for modelName in ( "CDM", "massLossFalse", "gamma0.0", "gamma0.5" ):
    log = open("outputs/test-decayingDM-profile-"+modelName+".log","w")
    status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/decayingDM-profile-"+modelName+".xml",stdout=log,stderr=log,shell=True)
    log.close()
    print("...done ("+str(status)+")")
    if status.returncode != 0:
        print("FAILED: model run ("+modelName+"):")
        subprocess.run("cat outputs/test-decayingDM-profile-"+modelName+".log",shell=True)
        sys.exit(0)
    print("Checking for errors...")
    status = subprocess.run("grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" outputs/test-decayingDM-profile-"+modelName+".log",shell=True)
    print("...done ("+str(status)+")")
    if status.returncode == 0:
        print("FAILED: model run (errors):")
        subprocess.run("cat outputs/test-decayingDM-profile-"+modelName+".log",shell=True)
        sys.exit(0)
    print("SUCCESS: model run ("+modelName+")")
    # Open the models and extract the mass profile and satellite bound mass.
    model                = h5py.File('outputs/decayingDM-profile-'+modelName+'.hdf5','r')
    nodes                = model['Outputs/Output1/nodeData']
    densities[modelName] = nodes['densityProfile'][:][0][0]
    if not math.isclose(densities[modelName],densitiesTarget[modelName],rel_tol=1e-2):
        print("FAILED: density mismatch for model '"+modelName+"': "+str(densities[modelName])+" vs. "+str(densitiesTarget[modelName]))
        failed = True
        
# No failures occurred, therefore, success.
if not failed:
    print("SUCCESS: density profiles as expected" )
    
