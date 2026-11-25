#!/usr/bin/env python3
import subprocess
import sys
import shutil
import h5py

# Take trees from the Bolshoi simulation and process into both Galacticus and IRATE formats.
# Andrew Benson (ported to Python)

# Create output folder.
subprocess.run("cd ..; mkdir -p testSuite/outputs", shell=True)

# Create a file in Galacticus format.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/mergerTreeFileBuildBolshoi.xml", shell=True)
if status.returncode != 0:
    print("FAILED: failed to make Galacticus-format merger tree file from Bolshoi merger tree")
    sys.exit(0)

# Check that the value of sigma_8 recorded in the file matches that requested. (Formerly, the cosmological model was given inside
# the task, and so was not seen by objects built by default - such as the transferred power spectrum from which sigma_8 is
# computed - so a different value was recorded.)
with h5py.File("../testSuite/outputs/bolshoiTestTreesGLC_in.hdf5", "r") as treeFile:
    sigma8 = treeFile["cosmology"].attrs["sigma_8"]
if abs(sigma8 / 0.82 - 1.0) > 1.0e-3:
    print(f"FAILED: sigma_8 recorded in the Bolshoi merger tree file is {sigma8}, but 0.82 was requested")
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
    print("iratevalidate is not installed - validation of the IRATE-format file will be skipped")

# Run Galacticus on file in Galacticus format.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/bolshoiTestTreesGLC.xml", shell=True)
if status.returncode != 0:
    print("FAILED: failed to run Galacticus on the Galacticus-format Bolshoi merger tree file")
    sys.exit(0)

print("SUCCESS: Bolshoi merger tree file build")
