#!/usr/bin/env python3
"""Check the table of projected half-light radii of inclined exponential disks used by the galaxySizesSDSS analyses.

The `outputAnalysisDistributionOperatorDiskSizeInclntn` class tabulates, as a function of cos(i), the half-light radius of an
exponential disk seen at inclination i relative to that seen face-on, and stores the table in the dynamic data cache. The
tabulation runs under OpenMP, and every thread used to share the one GSL root solver, so with an empty cache the table came out
corrupted, and was then stored. A corrupted table is not monotonic, and so could not be reversed: the model aborted with
"reversed table would not be monotonic", in this run and every later one (fixed October 2026).

This script runs a small model containing a galaxySizesSDSS analysis, several times, each with a private cache that lacks the
table (all other cached files are taken from the real cache), and compares the table computed by each run with an independent
calculation. It then corrupts the stored table and checks that a further run detects this and recomputes it.

Andrew Benson, Claude
"""
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

import h5py
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

# Number of runs with an empty cache, and the number of OpenMP threads used in each.
COUNT_RUNS   = 8
COUNT_THREAD = 16
# Tolerance on log10 of the radius ratio. The table is computed with a relative tolerance of 1e-6.
TOLERANCE    = 1.0e-5
# Grid of cos(i) used by the class.
COS_INCLINATION = np.linspace(1.0, 1.0e-3, 100)

# Independent calculation of the table. For an exponential disk of unit scale length seen at inclination i, the light inside the
# isophote of semi-major axis X is, as a fraction of the total,
#   F(X) = [2π cos(i)]⁻¹ ∫_0^{2π} dφ [1-exp(-aX)(1+aX)]/a², where a(φ) = √(sin²φ/cos²i + cos²φ),
# after integrating over the radial coordinate analytically.
def fractionLight(X, c):
    def integrand(phi):
        aSquared = (np.sin(phi)/c)**2+np.cos(phi)**2
        a        = np.sqrt(aSquared)
        return (-np.expm1(-a*X)-a*X*np.exp(-a*X))/aSquared
    return quad(integrand, 0.0, 2.0*np.pi, epsabs=0.0, epsrel=1.0e-12, limit=500)[0]/(2.0*np.pi*c)

radiusHalf = np.array([brentq(lambda X: fractionLight(X, c)-0.5, 1.0e-3, 1.0e3, xtol=1.0e-14, rtol=1.0e-14) for c in COS_INCLINATION])
tableExpected = np.log10(radiusHalf/radiusHalf[0])

# Locate the real cache. This script runs with the working directory set to testSuite/, while models are run from the
# repository root.
if "GALACTICUS_DYNAMIC_DATA_PATH" in os.environ:
    pathCache = os.environ["GALACTICUS_DYNAMIC_DATA_PATH"]
elif "GALACTICUS_DATA_PATH" in os.environ:
    pathCache = os.path.join(os.environ["GALACTICUS_DATA_PATH"], "dynamic")
else:
    print("SKIPPED: neither GALACTICUS_DYNAMIC_DATA_PATH nor GALACTICUS_DATA_PATH is set")
    sys.exit(0)
pathOutputs = os.path.abspath("outputs/diskSizeInclinationTable")
os.makedirs(pathOutputs, exist_ok=True)
nameTable = "diskExponentialInclinedHalfMassRadii.hdf5"

# Make a private cache linking to every entry of the real cache except the table (and its lock file).
def makeCache(path):
    subprocess.run(["rm", "-rf", path], check=True)
    os.makedirs(os.path.join(path, "galacticStructure"))
    if os.path.isdir(pathCache):
        for name in os.listdir(pathCache):
            if name != "galacticStructure":
                os.symlink(os.path.join(pathCache, name), os.path.join(path, name))
        pathStructure = os.path.join(pathCache, "galacticStructure")
        if os.path.isdir(pathStructure):
            for name in os.listdir(pathStructure):
                if not name.startswith(nameTable):
                    os.symlink(os.path.join(pathStructure, name), os.path.join(path, "galacticStructure", name))

# Build the parameter file: quickTest.xml, with outputs spanning the redshift range needed by the analysis.
tree = ET.parse("../parameters/quickTest.xml")
root = tree.getroot()
for name in ("outputFileName", "outputTimes", "mergerTreeOutputter", "outputAnalysis"):
    for element in root.findall(name):
        root.remove(element)
ET.SubElement(ET.SubElement(root, "outputTimes", value="list"), "redshifts", value="0.00 0.02 0.04 0.06 0.08 0.10")
ET.SubElement(root, "mergerTreeOutputter", value="analyzer")
ET.SubElement(ET.SubElement(root, "outputAnalysis", value="galaxySizesSDSS"), "distributionNumber", value="10")
ET.SubElement(root, "outputFileName", value=os.path.join(pathOutputs, "model.hdf5"))
fileParameters = os.path.join(pathOutputs, "diskSizeInclinationTable.xml")
tree.write(fileParameters)

def runModel(label, pathCachePrivate):
    environment = dict(os.environ, GALACTICUS_DYNAMIC_DATA_PATH=pathCachePrivate, OMP_NUM_THREADS=str(COUNT_THREAD))
    status = subprocess.run(["./Galacticus.exe", fileParameters], cwd="..", env=environment, capture_output=True)
    log = status.stdout.decode()+status.stderr.decode()
    with open(os.path.join(pathOutputs, f"{label}.log"), "w") as fileLog:
        fileLog.write(log)
    if status.returncode != 0:
        print(f"FAILED: model run '{label}' failed (log in testSuite/outputs/diskSizeInclinationTable/{label}.log):")
        print(log[-1500:])
        return None
    fileTable = os.path.join(pathCachePrivate, "galacticStructure", nameTable)
    if not os.path.exists(fileTable):
        print(f"FAILED: model run '{label}' did not store the table")
        return None
    with h5py.File(fileTable, "r") as table:
        return table["halfMassRadii"][:]

def check(label, table):
    if table is None:
        return
    deviation = np.abs(table-tableExpected)
    if deviation.max() < TOLERANCE:
        print(f"SUCCESS: table from run '{label}' agrees with the independent calculation (maximum deviation {deviation.max():.1e})")
    else:
        print(f"FAILED: table from run '{label}' deviates from the independent calculation by up to {deviation.max():.3g} (at index {deviation.argmax()+1})")

# Runs with an empty cache.
for i in range(COUNT_RUNS):
    label = f"empty{i+1}"
    pathCachePrivate = os.path.join(pathOutputs, f"cache{i+1}")
    makeCache(pathCachePrivate)
    check(label, runModel(label, pathCachePrivate))

# Corrupt the table stored by the last run (making it non-monotonic), and check that it is recomputed.
fileTable = os.path.join(pathCachePrivate, "galacticStructure", nameTable)
if os.path.exists(fileTable):
    with h5py.File(fileTable, "r+") as table:
        values = table["halfMassRadii"][:]
        values[[40, 41]] = values[[41, 40]]
        table["halfMassRadii"][...] = values
    check("corrupted", runModel("corrupted", pathCachePrivate))
sys.exit(0)
