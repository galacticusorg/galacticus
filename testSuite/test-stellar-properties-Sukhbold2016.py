#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check the Sukhbold et al. (2016) stellar properties compilation.
# Andrew Benson (10-August-2026); generated with assistance from Claude.

# Sukhbold et al. (2016) computed their models at a single metallicity, and the yields are written at two
# bracketing metallicities with identical values to say so explicitly. That encoding is not merely cosmetic: the
# irregular two dimensional interpolation used by `stellarAstrophysicsFile` extrapolates badly when the
# requested metallicity lies outside the range spanned by the tabulated points at the relevant masses, and for a
# table at a single metallicity that is almost always the case. Written at one metallicity these models return
# eleven times the expected metal yield, or zero, or hang the initial mass function integration, depending on
# what they are combined with.
#
# Roughly half the Sukhbold models collapse entirely to a black hole and eject only their wind, so this
# compilation must return substantially fewer metals than the standard one: their initial mass function weighted
# metal yield is 0.0156 against 0.0290 for Portinari, Chiosi & Bressan (1998). Diluted by the asymptotic giant
# branch and Type Ia contributions, which are common to both compilations, that implies a ratio of roughly 0.6
# in the metals a model produces. The test requires the ratio to sit in a band around that, which is what
# distinguishes a working interpolation from the failure modes above -- all of which land far outside it.

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

models = {}
for name, parameterFile in (
        ("standard", "supernovaeTypeIaYieldsDefault"      ),
        ("Sukhbold", "stellarPropertiesSukhbold2016"      ),
):
    logFileName = f"outputs/stellarPropertiesSukhbold{name}.log"
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

outputs = {
    "standard": "outputs/supernovaeTypeIaYieldsDefault.hdf5",
    "Sukhbold": "outputs/stellarPropertiesSukhbold2016.hdf5",
}
for name, fileName in outputs.items():
    with h5py.File(fileName, "r") as model:
        group    = model["Outputs"]
        nodeData = group[sorted(group.keys())[-1]]["nodeData"]
        models[name] = {
            "metals"      : np.array(nodeData["diskAbundancesStellarMetals"][:]).sum(),
            "massStellar" : np.array(nodeData["diskMassStellar"            ][:]).sum(),
        }

for name in models:
    if models[name]["metals"] > 0.0:
        print(f"SUCCESS: {name} compilation produces metals")
    else:
        print(f"FAILED: {name} compilation produces no metals ({models[name]['metals']}); with a "
              "single-metallicity table this is what an interpolation returns when it has nothing to "
              "interpolate between")
        sys.exit(0)

# The explodability prescription must reduce the metals produced, by roughly the factor implied by the initial
# mass function weighted yields. The band is wide enough to tolerate changes in the initial mass function or the
# other components of the compilation, but excludes every failure mode seen with an unbracketed table.
ratio = models["Sukhbold"]["metals"]/models["standard"]["metals"]
if 0.4 < ratio < 0.85:
    print(f"SUCCESS: the explodability prescription reduces the metals produced (ratio {ratio:.4f})")
else:
    print(f"FAILED: the Sukhbold compilation gives {ratio:.4f} times the metals of the standard one, outside "
          "the expected range of 0.4 to 0.85; the metallicity interpolation may be extrapolating")

# The stellar mass should be comparable: less recycling shifts it a little, but not greatly.
massRatio = models["Sukhbold"]["massStellar"]/models["standard"]["massStellar"]
if 0.9 < massRatio < 1.1:
    print(f"SUCCESS: stellar masses are comparable (ratio {massRatio:.4f})")
else:
    print(f"FAILED: stellar masses differ by more than 10% (ratio {massRatio:.4f}), which suggests the "
          "compilation is not supplying sensible ejected masses")

sys.exit(0)
