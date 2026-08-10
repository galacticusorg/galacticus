#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check that Type Ia supernova yields can be made to depend on the progenitor metallicity.
# Andrew Benson (09-August-2026); generated with assistance from Claude.

# A Type Ia yields file may wrap each set of isotopes in a `yieldsMetallicity` element carrying the metallicity
# at which that set applies, in which case Galacticus interpolates between them, holding the yield constant
# beyond the tabulated range. Seitenzahl et al. (2013) post-processed their N100 explosion model at four
# progenitor metallicities, so that sequence isolates the effect: only the metallicity differs between the four
# sets, not the explosion model.
#
# Testing this needs care, because the absolute signal is small. A Chandrasekhar-mass white dwarf burns
# essentially completely, so the total Type Ia metal yield varies by only about 0.2% across the sequence, and
# what really changes is the composition. Meanwhile Type Ia supernovae supply only a small share of the metals
# in the model, the rest coming from the stellar yield tables, which are identical between these runs.
#
# The test therefore runs the metallicity-dependent file alongside two single-metallicity files taken from the
# ends of the same sequence, and forms
#
#     f = (dependent - lowest) / (Solar - lowest)
#
# The stellar contribution is common to all three runs and cancels in both differences, so f measures only where
# the metallicity-dependent run sits between the two tabulated extremes. It would be 0 if the reader always took
# the lowest tabulated set, and 1 if it always took the highest; a value strictly between the two is the
# signature of the yields actually being interpolated as the modelled galaxy enriches.

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

models = {}
for name in "MetallicityDependent", "SolarOnly", "LowMetallicity":
    logFileName = f"outputs/supernovaeTypeIaYields{name}.log"
    with open(logFileName, "w") as logFile:
        status = subprocess.run(
            f"cd ..; ./Galacticus.exe testSuite/parameters/supernovaeTypeIaYields{name}.xml",
            shell=True, stdout=logFile, stderr=subprocess.STDOUT
        )
    if status.returncode != 0:
        print(f"FAILED: {name} yields model run:")
        with open(logFileName) as logFile:
            print(logFile.read())
        sys.exit(0)
    print(f"SUCCESS: {name} yields model run")
    with h5py.File(f"outputs/supernovaeTypeIaYields{name}.hdf5", "r") as model:
        outputs  = model["Outputs"]
        nodeData = outputs[sorted(outputs.keys())[-1]]["nodeData"]
        models[name] = np.array(nodeData["diskAbundancesStellarMetals"][:]).sum()

# All three must produce metals. A file whose `yieldsMetallicity` elements failed to parse would supply none.
for name in models:
    if models[name] > 0.0:
        print(f"SUCCESS: {name} yields produce metals")
    else:
        print(f"FAILED: {name} yields produce no metals ({models[name]})")
        sys.exit(0)

# The two single-metallicity files must give different results, or f below is undefined.
separation = models["SolarOnly"]-models["LowMetallicity"]
if abs(separation)/models["SolarOnly"] > 1.0e-9:
    print(f"SUCCESS: the two tabulated extremes give distinguishable results "
          f"(fractional separation {abs(separation)/models['SolarOnly']:.3e})")
else:
    print(f"FAILED: the lowest-metallicity and Solar files give indistinguishable results "
          f"(fractional separation {abs(separation)/models['SolarOnly']:.3e}), so there is nothing to "
          "interpolate between")
    sys.exit(0)

# Where does the metallicity-dependent run sit between them?
fraction = (models["MetallicityDependent"]-models["LowMetallicity"])/separation
if 0.05 < fraction < 0.95:
    print(f"SUCCESS: metallicity-dependent yields interpolate between the tabulated sets (f = {fraction:.4f})")
elif fraction <= 0.05:
    print(f"FAILED: metallicity-dependent yields match the lowest tabulated set (f = {fraction:.4f}), so the "
          "reader appears to ignore metallicity or always take the first set")
else:
    print(f"FAILED: metallicity-dependent yields match the highest tabulated set (f = {fraction:.4f}), so the "
          "reader appears to ignore metallicity or always take the last set")

sys.exit(0)
