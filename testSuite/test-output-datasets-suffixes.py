#!/usr/bin/env python3
import subprocess
import h5py
import re

# Check that suffixes are correctly added to dataset names.
# Andrew Benson (06-September-2024)

# Run the model and check for completion.
print("Running model...")
status = subprocess.run("mkdir -p outputs",shell=True)
log = open("outputs/test-output-datasets-suffixes.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/outputDatasetSuffixes.xml",stdout=log,stderr=log,shell=True)
log.close()
if status.returncode != 0:
    print("   ...done ("+str(status)+")")
    print("   FAILED: model run:")
    subprocess.run("cat outputs/test-output-datasets-suffixes.log",shell=True)
    sys.exit()
else:
    print("   ...done")
    print("   Checking for errors...")
    status = subprocess.run("grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" outputs/test-output-datasets-suffixes.log",shell=True)
    if status.returncode == 0:
        print("   ...done ("+str(status)+")")
        print("   FAILED: model run (errors):")
        subprocess.run("cat outputs/test-output-datasets-suffixes.log",shell=True)
        sys.exit()
    else:
        print("   ...done")
        print("   SUCCESS: model run")

# Open the model and look for duplicated datasets.
model  = h5py.File('outputs/outputDatasetSuffixes.hdf5','r')
nodes  = model['Outputs/Output1/nodeData']
status = "SUCCESS"
for name in nodes:
    matched = re.match(r'(.*)_duplicate',name)
    if matched:
        nameOriginal = matched.group(1)
        if not nameOriginal in nodes:
            status = "FAILED"

# Report status.
print(status+": duplicated dataset names found")
