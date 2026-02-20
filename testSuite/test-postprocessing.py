#!/usr/bin/env python3
import subprocess
import sys
import h5py
import re
import os
import numpy as np

# Run a simple test of postprocssing.
# Andrew Benson (19-November-2025)

# Run the original model and check for completion.
print("   Running original model...")
status = subprocess.run("mkdir -p outputs",shell=True)
log = open("outputs/test-postprocessing-original.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/test-postprocessing-original.xml",stdout=log,stderr=log,shell=True)
log.close()
if status.returncode != 0:
    print("   ...done ("+str(status)+")")
    print("   FAILED: model run:")
    subprocess.run("cat outputs/test-postprocessing-original.log",shell=True)
    sys.exit()
else:
    print("   ...done")
    print("   Checking for errors...")
    status = subprocess.run("grep -q -i -e fatal -e aborted -e \"task failed\" -e \"Galacticus experienced an error in the GSL library\" outputs/test-postprocessing-original.log",shell=True)
    if status.returncode == 0:
        print("   ...done ("+str(status)+")")
        print("   FAILED: model run (errors):")
        subprocess.run("cat outputs/test-postprocessing-original.log",shell=True)
        sys.exit()
    else:
        print("   ...done")
        print("   SUCCESS: model run")

# Run the postprocessing model and check for completion.
print("   Running postprocessed model...")
status = subprocess.run("mkdir -p outputs",shell=True)
log = open("outputs/test-postprocessing-postprocessed.log","w")
status = subprocess.run("cd ..; export OMP_NUM_THREADS=1; ./Galacticus.exe testSuite/parameters/test-postprocessing-original.xml testSuite/parameters/test-postprocessing-perform.xml",stdout=log,stderr=log,shell=True)
log.close()
if status.returncode != 0:
    print("   ...done ("+str(status)+")")
    print("   FAILED: model run:")
    subprocess.run("cat outputs/test-postprocessing-postprocessed.log",shell=True)
    sys.exit()
else:
    print("   ...done")
    print("   Checking for errors...")
    status = subprocess.run("grep -q -i -e fatal -e aborted -e \"task failed\" -e \"Galacticus experienced an error in the GSL library\" outputs/test-postprocessing-postprocessed.log",shell=True)
    if status.returncode == 0:
        print("   ...done ("+str(status)+")")
        print("   FAILED: model run (errors):")
        subprocess.run("cat outputs/test-postprocessing-postprocessed.log",shell=True)
        sys.exit()
    else:
        print("   ...done")
        print("   SUCCESS: model run")

# Open the model and compare.
modelOriginal        = h5py.File('outputs/postprocessing-original.hdf5'     ,'r')
modelPostprocessed   = h5py.File('outputs/postprocessing-postprocessed.hdf5','r')
outputsOriginal      = modelOriginal     ['Outputs']
outputsPostprocessed = modelPostprocessed['Outputs']
# Iterate over outputs.
for (outputName, outputOriginal) in outputsOriginal.items():
    match = re.match(r'^Output(\d+)',outputName) 
    if not match:
        continue
    outputPostprocessed = outputsPostprocessed[outputName]
    nodesOriginal       = outputOriginal      ['nodeData']
    nodesPostprocessed  = outputPostprocessed ['nodeData']
    # Iterate over datasets.
    for (datasetName, datasetOriginal) in nodesOriginal.items():
        if not datasetName in nodesPostprocessed:
            print(f"FAILED: dataset '{datasetName}' does not exist in postprocessed file output '{outputName}'")
            sys.exit()
        datasetPostprocessed = nodesPostprocessed[datasetName]
        if len(datasetPostprocessed[:]) != len(datasetOriginal[:]):
            print(f"FAILED: dataset '{datasetName}' differs in length in postprocessed file output '{outputName}'")
            sys.exit()
        if not np.allclose(datasetOriginal[:],datasetPostprocessed[:],rtol=1.0e-3,atol=1.0e-30):
            print(f"FAILED: dataset '{datasetName}' differs in postprocessed file output '{outputName}'")
            sys.exit()
    # Check that the new dataset appears in the postprocessed file.
    if not 'luminosityEmissionLineDisk:oxygenII3727' in nodesPostprocessed:
        print(f"FAILED: new dataset 'luminosityEmissionLineDisk:oxygenII3727' not present in postprocessed file output '{outputName}'")
        sys.exit()
  
# Report success.
print("SUCCESS: postprocessing")
