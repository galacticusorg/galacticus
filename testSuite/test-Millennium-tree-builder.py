#!/usr/bin/env python3
import subprocess
import sys
import shutil

# Take trees from the Millennium Database and process into both Galacticus and IRATE formats.
# Andrew Benson (ported to Python)

# Create a file in Galacticus format.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/mergerTreeFileBuildMillennium.xml", shell=True)
if status.returncode != 0:
    print("FAILED: failed to make Galacticus-format merger tree file from Millennium database output")
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
    print("SKIPPED: iratevalidate is not installed - validation of IRATE-format file will be skipped")
