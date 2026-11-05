#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check the compilation of stellar properties combining the Karakas (2010) asymptotic giant branch, Doherty et al.
# (2014) super-AGB and Limongi & Chieffi (2018) massive star models.
# Andrew Benson (22-October-2026); generated with assistance from Claude.
#
# The compilation differs from the Limongi & Chieffi (2018) compilation only in its intermediate mass stars: both
# take their massive star yields from the same Limongi & Chieffi models, while this one replaces the Marigo (2001)
# AGB models with those of Karakas (2010) and Doherty et al. (2014). The test runs the isolated star-forming disk
# model used by `test-stellar-properties-LimongiChieffi2018.py` with each of the two compilations, using
# non-instantaneous recycling (so that the element yields of the stellar models are used) and tracking O and Fe.
#
# Since the two compilations share their massive stars and Type Ia supernovae, which are the sources of iron, they
# must produce the same mass of iron: AGB and super-AGB stars neither make nor destroy significant amounts of it
# (they eject iron only at their birth composition, which is not part of the net yield). The two agree to 0.4% in
# practice. A disagreement would indicate that the intermediate mass star data were injecting iron - the failure
# seen when yields were interpolated across gaps in the tabulated stars, which gave 2-8 Msun stars at low
# metallicity up to 0.15 Msun of iron each. The tolerance of 2% is five times the difference expected.
#
# The total metal yields must differ, since the AGB models do, but the size of that difference is a scientific
# result rather than something to pin down here.

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Construct the parameter file for the Limongi & Chieffi (2018) compilation from that for the new compilation, so
# that the two runs differ only in the stellar data used.
nameNew = "stellarPropertiesKarakas2010Doherty2014LimongiChieffi2018"
nameOld = "stellarPropertiesLimongiChieffi2018NonInstantaneous"
with open(f"parameters/{nameNew}.xml") as file:
    parameters = file.read()
parameters = parameters.replace("stellarPropertiesCompilationKarakas2010Doherty2014LimongiChieffi2018.xml",
                                "stellarPropertiesCompilationLimongiChieffi2018.xml"                     )
parameters = parameters.replace(f"testSuite/outputs/{nameNew}.hdf5",f"testSuite/outputs/{nameOld}.hdf5")
with open(f"outputs/{nameOld}.xml","w") as file:
    file.write(parameters)

models = {}
for name, parameterFile in (
        ("KarakasDoherty", f"testSuite/parameters/{nameNew}.xml"),
        ("LimongiChieffi", f"testSuite/outputs/{nameOld}.xml"   ),
):
    logFileName = f"outputs/stellarProperties{name}NonInstantaneous.log"
    with open(logFileName, "w") as logFile:
        status = subprocess.run(
            f"cd ..; ./Galacticus.exe {parameterFile}",
            shell=True, stdout=logFile, stderr=subprocess.STDOUT
        )
    if status.returncode != 0:
        print(f"FAILED: {name} compilation model run:")
        with open(logFileName) as logFile:
            print(logFile.read())
        sys.exit(0)
    print(f"SUCCESS: {name} compilation model run")

outputFiles = {
    "KarakasDoherty": f"outputs/{nameNew}.hdf5",
    "LimongiChieffi": f"outputs/{nameOld}.hdf5",
}
for name, fileName in outputFiles.items():
    with h5py.File(fileName, "r") as model:
        outputs  = model["Outputs"]
        nodeData = outputs[sorted(outputs.keys())[-1]]["nodeData"]
        total    = lambda quantity: np.array(nodeData[f"diskAbundancesStellar{quantity}"][:]).sum() \
                                   +np.array(nodeData[f"diskAbundancesGas{quantity}"    ][:]).sum()
        models[name] = {
            "metals"      : total("Metals"),
            "iron"        : total("Fe"    ),
            "massStellar" : np.array(nodeData["diskMassStellar"][:]).sum(),
        }

# Both compilations must produce metals and iron. A compilation which silently failed to supply yields -- because
# an `XInclude` did not resolve, say -- would return zero here rather than raising an error.
for name in models:
    for quantity in ("metals", "iron"):
        if models[name][quantity] > 0.0:
            print(f"SUCCESS: {name} compilation produces {quantity}")
        else:
            print(f"FAILED: {name} compilation produces no {quantity} ({models[name][quantity]})")

# Both must form a comparable mass of stars.
massRatio = models["KarakasDoherty"]["massStellar"]/models["LimongiChieffi"]["massStellar"]
if 0.9 < massRatio < 1.1:
    print(f"SUCCESS: stellar masses are comparable (ratio {massRatio:.4f})")
else:
    print(f"FAILED: stellar masses differ by more than 10% (ratio {massRatio:.4f}), which suggests the "
          "compilation is not supplying sensible ejected masses")

# The AGB models differ, so the metal yields must differ. Were they identical, the compilation would not actually
# be selecting the Karakas (2010) and Doherty et al. (2014) data.
metalRatio = models["KarakasDoherty"]["metals"]/models["LimongiChieffi"]["metals"]
if abs(metalRatio-1.0) > 1.0e-3:
    print(f"SUCCESS: the two compilations give different metal yields (ratio {metalRatio:.4f})")
else:
    print(f"FAILED: the two compilations give the same metal yield (ratio {metalRatio:.4f}), so the Karakas & "
          "Doherty compilation is not being used")

# The two compilations share their sources of iron, so must produce the same mass of it.
ironRatio = models["KarakasDoherty"]["iron"]/models["LimongiChieffi"]["iron"]
if abs(ironRatio-1.0) < 0.02:
    print(f"SUCCESS: the two compilations produce the same iron (ratio {ironRatio:.4f})")
else:
    print(f"FAILED: the two compilations produce different masses of iron (ratio {ironRatio:.4f}), although they "
          "share the massive star and Type Ia supernova models which make it; the intermediate mass star data "
          "may be injecting iron")

sys.exit(0)
