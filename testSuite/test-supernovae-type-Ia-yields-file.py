#!/usr/bin/env python3
import os
import re
import subprocess
import sys
import h5py
import numpy as np

# Check that the Type Ia supernova yield table can be selected at run time via the `fileName` parameter.
# Andrew Benson (09-August-2026); generated with assistance from Claude.

# The `supernovaeTypeIa` classes which inherit from `supernovaeTypeIaFixedYield` read their yields from a file
# named by the `fileName` parameter. This test verifies that the parameter is actually honored -- that the named
# file is read and affects the metals produced -- rather than being silently ignored.
#
# Two otherwise identical models are run: one using the default yield table, and one using a copy of that table
# in which every isotope yield has been doubled. Doubling the Type Ia yields must increase the metals produced
# while leaving the stellar mass untouched, since the yields do not feed back on star formation in these models.
# That second condition is the important one: it distinguishes "the yield table was read" from "the two models
# differ for some unrelated reason".
#
# Note that this needs a model which actually evaluates the yield integral. Parameter files which set
# `recycledFraction` and `metalYield` explicitly on `stellarPopulation` -- `parameters/quickTest.xml` among them
# -- never read the yield table at all, and so cannot exercise this path.

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Build the doubled-yield table from the default table.
dataPath = os.environ.get("GALACTICUS_DATA_PATH")
if dataPath is None:
    print("FAILED: the GALACTICUS_DATA_PATH environment variable is not set")
    sys.exit(0)
yieldsDefault = os.path.join(dataPath, "static", "stellarAstrophysics", "Supernovae_Type_Ia_Yields.xml")
if not os.path.isfile(yieldsDefault):
    print(f"FAILED: default Type Ia yields file '{yieldsDefault}' not found")
    sys.exit(0)
with open(yieldsDefault) as yieldsFile:
    yieldsText = yieldsFile.read()
yieldsDoubled, countYields = re.subn(
    r"<yield>([^<]+)</yield>",
    lambda match: f"<yield>{2.0*float(match.group(1)):.6e}</yield>",
    yieldsText,
)
if countYields == 0:
    print(f"FAILED: no yields found in '{yieldsDefault}'")
    sys.exit(0)
with open("outputs/supernovaeTypeIaYieldsDoubled.xml", "w") as doubledFile:
    doubledFile.write(yieldsDoubled)
print(f"SUCCESS: built doubled yield table ({countYields} isotopes)")

# Run both models.
models = {}
for name in "Default", "Doubled":
    logFileName = f"outputs/supernovaeTypeIaYields{name}.log"
    with open(logFileName, "w") as logFile:
        status = subprocess.run(
            f"cd ..; ./Galacticus.exe testSuite/parameters/supernovaeTypeIaYields{name}.xml",
            shell=True, stdout=logFile, stderr=subprocess.STDOUT
        )
    if status.returncode != 0:
        print(f"FAILED: {name.lower()} yields model run:")
        with open(logFileName) as logFile:
            print(logFile.read())
        sys.exit(0)
    print(f"SUCCESS: {name.lower()} yields model run")
    with h5py.File(f"outputs/supernovaeTypeIaYields{name}.hdf5", "r") as model:
        outputs   = model["Outputs"]
        nodeData  = outputs[sorted(outputs.keys())[-1]]["nodeData"]
        models[name] = {
            "metals"      : np.array(nodeData["diskAbundancesGasMetals"    ][:]).sum()
                           +np.array(nodeData["diskAbundancesStellarMetals"][:]).sum(),
            "massStellar" : np.array(nodeData["diskMassStellar"            ][:]).sum(),
        }

# The default model must actually produce metals -- otherwise the comparison below is vacuous.
if models["Default"]["metals"] > 0.0:
    print("SUCCESS: default model produces metals")
else:
    print(f"FAILED: default model produces no metals ({models['Default']['metals']})")
    sys.exit(0)

# Doubling the Type Ia yields must increase the metals produced. The increase is the Type Ia share of all metal
# production, so it is well below a factor of two; the bounds here are loose enough to tolerate changes in the
# default initial mass function or stellar yields, but far above any numerical noise (the two models are
# otherwise identical, and deterministic).
increase = (models["Doubled"]["metals"]-models["Default"]["metals"])/models["Default"]["metals"]
if 0.005 < increase < 0.5:
    print(f"SUCCESS: doubling Type Ia yields increases metals (by {100.0*increase:.3f}%)")
else:
    print(f"FAILED: doubling Type Ia yields changed metals by {100.0*increase:.3f}%, "
          "which is outside the expected range of 0.5% to 50%")

# Stellar mass must be unaffected: the yields do not feed back on star formation here, so any change would mean
# the two models differ in something other than the yield table.
massDifference = abs(models["Doubled"]["massStellar"]-models["Default"]["massStellar"]) \
                 /models["Default"]["massStellar"]
if massDifference < 1.0e-6:
    print("SUCCESS: stellar mass is unchanged by the yield table")
else:
    print(f"FAILED: stellar mass changed by {100.0*massDifference:.3e}% between the two models, so they differ "
          "in more than the Type Ia yield table")

sys.exit(0)
