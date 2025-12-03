#!/usr/bin/env python3
import subprocess
import re

# Check that strict handling of parameter files works as intended.
# Andrew Benson (03-December-2025)

# Iterate over models.
for model in ( "strictOutdated", "unstrictOutdated", "strictUnrecognized", "unstrictUnrecognized" ):
    # Run the model and check for completion.
    print(f"Running model '{model}'...")
    status = subprocess.run("mkdir -p outputs",shell=True)
    log = open(f"outputs/test-{model}.log","w")
    status = subprocess.run(f"cd ..; ./Galacticus.exe testSuite/parameters/{model}.xml",stdout=log,stderr=log,shell=True)
    log.close()
    print("...done ("+str(status)+")")
    if status.returncode == 0:
        if re.match(r'strict',model):
            # Strict model, succeeded - not expected.
            print(f"FAIL: model {model} succeeded")
        else:
            # Unstrict model, succeeded - expected.
            print(f"SUCCESS: model {model} succeeded")
    else:
        if re.match(r'strict',model):
            # Strict model, failed - expected.
            print(f"SUCCESS: model {model} failed")
        else:
            # Unstrict model, failed - not expected.
            print(f"FAIL: model {model} failed")

