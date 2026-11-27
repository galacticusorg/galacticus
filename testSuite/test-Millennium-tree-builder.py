#!/usr/bin/env python3
import subprocess
import sys
import shutil
import h5py

# Take trees from the Millennium Database and process into both Galacticus and IRATE formats.
# Andrew Benson (ported to Python)

# Create a file in Galacticus format.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/mergerTreeFileBuildMillennium.xml", shell=True)
if status.returncode != 0:
    print("FAILED: failed to make Galacticus-format merger tree file from Millennium database output")
    sys.exit(0)

# Check that the value of sigma_8 recorded in the file matches that requested. (Formerly, the cosmological model was given inside
# the task, and so was not seen by objects built by default - such as the transferred power spectrum from which sigma_8 is
# computed - so a different value was recorded.)
with h5py.File("../testSuite/outputs/millenniumTestTreesGLC.hdf5", "r") as treeFile:
    sigma8 = treeFile["cosmology"].attrs["sigma_8"]
    simulation = dict(treeFile["simulation"].attrs)
if abs(sigma8 / 0.9 - 1.0) > 1.0e-3:
    print(f"FAILED: sigma_8 recorded in the Millennium merger tree file is {sigma8}, but 0.9 was requested")
    sys.exit(0)
# Check that each simulation metadata entry was recorded. (Formerly, a malformed entry in the parameter file caused the starting
# redshift to be silently dropped.)
if simulation.get("startRedshift") != 127.0:
    print(f"FAILED: startRedshift recorded in the Millennium merger tree file is {simulation.get('startRedshift')}, but 127 was given")
    sys.exit(0)

# Extract a single tree from the file we just created.
status = subprocess.run("cd ..; ./scripts/aux/extractSingleTree.py testSuite/outputs/millenniumTestTreesGLC.hdf5 testSuite/outputs/millenniumTestTreesGLCSingle.hdf5 79000000", shell=True)
if status.returncode != 0:
    print("FAILED: failed to extract a single merger tree")
    sys.exit(0)

# Run the single tree to verify that the file is valid.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/extractSingleTreeRun.xml", shell=True)
if status.returncode != 0:
    print("FAILED: failed to run single merger tree")
    sys.exit(0)

# Create a file in IRATE format.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/mergerTreeFileBuildMillenniumIRATE.xml", shell=True)
if status.returncode != 0:
    print("FAILED: failed to make IRATE-format merger tree file from Millennium database output")
    sys.exit(0)

# Validate the IRATE format file.
validator = shutil.which("iratevalidate")
if validator:
    status = subprocess.run("cd ..; iratevalidate testSuite/outputs/millenniumTestTreesIRATE.hdf5", shell=True)
    if status.returncode != 0:
        print("FAILED: IRATE-format file generated from Millennium database output did not validate")
        sys.exit(0)
else:
    print("iratevalidate is not installed - validation of the IRATE-format file will be skipped")

print("SUCCESS: Millennium merger tree file build")
