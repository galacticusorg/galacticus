#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check that the Limongi & Chieffi (2018) stellar properties compilation can be used in place of the standard one.
# Andrew Benson (09-August-2026); generated with assistance from Claude.

# The compilation combines the Limongi & Chieffi (2018) massive star models (13-120 Msun, extending down to
# Z = 3.2e-5) with the Marigo (2001) asymptotic giant branch models and the Portinari, Chiosi & Bressan (1998)
# lifetimes. It exercises several things which the standard compilation does not: `XInclude` of a newly split
# file, yields tabulated for all 53 elements from hydrogen to bismuth, and a metallicity range extending well
# below that of the standard compilation.
#
# The test runs the same isolated star-forming disk model used by `test-supernovae-type-Ia-yields-file.py`, once
# with each compilation, and checks that both produce metals and that they differ. It deliberately does not
# assert a particular size of difference: the two compilations use different stellar models, so the yields are
# expected to differ, but by an amount which is a scientific result rather than something to pin down here.

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

models = {}
for name, parameterFile in (
        ("standard"      , "supernovaeTypeIaYieldsDefault"       ),
        ("LimongiChieffi", "stellarPropertiesLimongiChieffi2018"  ),
):
    logFileName = f"outputs/stellarProperties{name}.log"
    with open(logFileName, "w") as logFile:
        status = subprocess.run(
            f"cd ..; ./Galacticus.exe testSuite/parameters/{parameterFile}.xml",
            shell=True, stdout=logFile, stderr=subprocess.STDOUT
        )
    if status.returncode != 0:
        print(f"FAILED: {name} compilation model run:")
        with open(logFileName) as logFile:
            print(logFile.read())
        sys.exit(0)
    print(f"SUCCESS: {name} compilation model run")

outputFiles = {
    "standard"      : "outputs/supernovaeTypeIaYieldsDefault.hdf5",
    "LimongiChieffi": "outputs/stellarPropertiesLimongiChieffi2018.hdf5",
}
for name, fileName in outputFiles.items():
    with h5py.File(fileName, "r") as model:
        outputs  = model["Outputs"]
        nodeData = outputs[sorted(outputs.keys())[-1]]["nodeData"]
        models[name] = {
            "metals"      : np.array(nodeData["diskAbundancesStellarMetals"][:]).sum(),
            "massStellar" : np.array(nodeData["diskMassStellar"            ][:]).sum(),
        }

# Both compilations must produce metals. A compilation which silently failed to supply yields -- because an
# `XInclude` did not resolve, say -- would return zero here rather than raising an error.
for name in models:
    if models[name]["metals"] > 0.0:
        print(f"SUCCESS: {name} compilation produces metals")
    else:
        print(f"FAILED: {name} compilation produces no metals ({models[name]['metals']})")

# Both must form a comparable mass of stars: the stellar mass is set by the star formation law and the recycled
# fraction, and while the two compilations give slightly different recycling they should not differ greatly.
massRatio = models["LimongiChieffi"]["massStellar"]/models["standard"]["massStellar"]
if 0.9 < massRatio < 1.1:
    print(f"SUCCESS: stellar masses are comparable (ratio {massRatio:.4f})")
else:
    print(f"FAILED: stellar masses differ by more than 10% (ratio {massRatio:.4f}), which suggests the "
          "compilation is not supplying sensible ejected masses")

# The two compilations use different stellar models, so their yields must differ. Were they identical, the
# compilation would not actually be selecting the Limongi & Chieffi data.
metalRatio = models["LimongiChieffi"]["metals"]/models["standard"]["metals"]
if abs(metalRatio-1.0) > 1.0e-3:
    print(f"SUCCESS: the two compilations give different metal yields (ratio {metalRatio:.4f})")
else:
    print(f"FAILED: the two compilations give the same metal yield (ratio {metalRatio:.4f}), so the Limongi & "
          "Chieffi compilation is not being used")

sys.exit(0)
