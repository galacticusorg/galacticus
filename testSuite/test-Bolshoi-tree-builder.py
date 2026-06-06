#!/usr/bin/env python3
import subprocess
import sys
import shutil

# Take trees from the Bolshoi simulation and process into both Galacticus and IRATE formats.
# Andrew Benson (ported to Python)

# Create output folder.
subprocess.run("cd ..; mkdir -p testSuite/outputs", shell=True)

# Create a file in Galacticus format.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/mergerTreeFileBuildBolshoi.xml", shell=True)
if status.returncode != 0:
    print("FAILED: failed to make Galacticus-format merger tree file from Bolshoi merger tree")
    sys.exit(0)

# Create a file in IRATE format.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/mergerTreeFileBuildBolshoiIRATE.xml", shell=True)
if status.returncode != 0:
    print("FAILED: failed to make IRATE-format merger tree file from Bolshoi merger tree")
    sys.exit(0)

# Validate the IRATE format file.
validator = shutil.which("iratevalidate")
if validator:
    status = subprocess.run("cd ..; iratevalidate testSuite/outputs/bolshoiTestTreesIRATE_in.hdf5", shell=True)
    if status.returncode != 0:
        print("FAILED: IRATE-format file generated from Bolshoi merger tree did not validate")
        sys.exit(0)
else:
    print("SKIP: iratevalidate is not installed - validation of IRATE-format file will be skipped")

# Run Galacticus on file in Galacticus format.
subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/bolshoiTestTreesGLC.xml", shell=True)
