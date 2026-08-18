#!/usr/bin/env python3
import subprocess
import re
import os

# Check that strict handling of parameter files works as intended.
# Andrew Benson (03-December-2025)

# Specify models.
models = [ "strictUnrecognized", "unstrictUnrecognized" ]
## Determine if we have the Git2 library and so can test outdated parameter files.
log = open(f"outputs/test-strict-git2.log","w")
status = subprocess.run(f"cd ..; ./Galacticus.exe parameters/report.xml | grep -q GIT2AVAIL",stdout=log,stderr=log,shell=True)
log.close()
haveGit2 = status.returncode == 0
if haveGit2:
    models.append("strictOutdated"  )
    models.append("unstrictOutdated")

# Report the model log on an unexpected outcome.  The `testSuite/outputs`
# directory is not retained by the pull-request checks, so a model that fails
# in CI leaves nothing to diagnose it with unless its log is echoed here.
def reportLog(fileName):
    print(f"vvvvv begin {fileName} vvvvv")
    try:
        lines = open(fileName,"r",errors="replace").readlines()
    except OSError as error:
        lines = [f"(unable to read log: {error})\n"]
    # Keep both ends: an unexpected *success* is explained by what the model
    # did not complain about at startup, an unexpected *failure* by how it
    # ended.
    if len(lines) > 300:
        lines = lines[:50]+[f"... {len(lines)-300} lines omitted ...\n"]+lines[-250:]
    for line in lines:
        print(line.rstrip("\n"))
    print(f"^^^^^ end {fileName} ^^^^^")

# Iterate over models.
for model in models:
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
            reportLog(f"outputs/test-{model}.log")
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
            reportLog(f"outputs/test-{model}.log")

