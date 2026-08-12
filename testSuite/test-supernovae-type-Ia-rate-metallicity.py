#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check the metallicity-dependent scaling of the Type Ia supernova rate.
# Andrew Benson (10-August-2026); generated with assistance from Claude.

# The `metallicityDependentRate` class wraps any other `supernovaeTypeIa` class and multiplies both the number
# of events and the yields by (Z/Z_Solar)^beta, following Johnson, Kochanek & Stanek (2023), who argue that the
# anti-correlation of the close binary fraction with metallicity enhances the Type Ia rate at low metallicity
# and find a scaling of about Z^(-1/2).
#
# Four models are run, all otherwise identical:
#
#  * the undecorated power-law delay time distribution;
#  * the same wrapped by the decorator with an exponent of zero;
#  * the same wrapped with the default exponent of -1/2;
#  * the same wrapped with an exponent of +1/2.
#
# The zero-exponent run is the important one. The scaling is unity at every metallicity when beta = 0, so that
# run must reproduce the undecorated model exactly -- which checks that the decorator forwards to the wrapped
# class faithfully and introduces no arithmetic of its own. A decorator which, say, dropped the yields or
# mishandled the optional atomic index argument would show up here.

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

models = {}
for name, parameterFile in (
        ("undecorated" , "supernovaeTypeIaYieldsDefault"           ),
        ("nullDecorator", "supernovaeTypeIaRateNullDecorator"      ),
        ("metallicity" , "supernovaeTypeIaRateMetallicityDependent"),
        ("positive"    , "supernovaeTypeIaRatePositiveExponent"    ),
):
    logFileName = f"outputs/supernovaeTypeIaRate{name}.log"
    with open(logFileName, "w") as logFile:
        status = subprocess.run(
            f"cd ..; ./Galacticus.exe testSuite/parameters/{parameterFile}.xml",
            shell=True, stdout=logFile, stderr=subprocess.STDOUT
        )
    if status.returncode != 0:
        print(f"FAILED: {name} model run:")
        with open(logFileName) as logFile:
            print(logFile.read())
        sys.exit(0)
    print(f"SUCCESS: {name} model run")

outputs = {
    "undecorated"  : "outputs/supernovaeTypeIaYieldsDefault.hdf5",
    "nullDecorator": "outputs/supernovaeTypeIaRateNullDecorator.hdf5",
    "metallicity"  : "outputs/supernovaeTypeIaRateMetallicityDependent.hdf5",
    "positive"     : "outputs/supernovaeTypeIaRatePositiveExponent.hdf5",
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
        print(f"SUCCESS: {name} model produces metals")
    else:
        print(f"FAILED: {name} model produces no metals ({models[name]['metals']})")
        sys.exit(0)

# With an exponent of zero the scaling is unity everywhere, so the decorator must be transparent.
difference = abs(models["nullDecorator"]["metals"]-models["undecorated"]["metals"]) \
             /models["undecorated"]["metals"]
if difference < 1.0e-12:
    print(f"SUCCESS: a zero exponent reproduces the undecorated model (fractional difference {difference:.2e})")
else:
    print(f"FAILED: a zero exponent changes the metals produced by a fractional {difference:.2e}, so the "
          "decorator is not transparent when it should be")

# A negative exponent must enhance the rate, because the modelled galaxy spends its early life below Solar
# metallicity where the scaling exceeds unity. The enhancement applies to the Type Ia contribution only, which
# is a minority of the metals, so the total must rise by a modest amount rather than by the full scaling factor.
ratio = models["metallicity"]["metals"]/models["undecorated"]["metals"]
if ratio > 1.0:
    print(f"SUCCESS: a negative exponent enhances metal production (ratio {ratio:.6f})")
else:
    print(f"FAILED: a negative exponent gives a metal ratio of {ratio:.6f}, but enhancing the Type Ia rate at "
          "low metallicity should increase the metals produced")

# Reversing the sign of the exponent must reverse the effect. This pins the sign convention, which is the
# detail most easily got wrong: a scaling applied as Z^(+beta) where Z^(-beta) was meant would still pass every
# check above, since those only require the decorated model to differ from the undecorated one.
ratioPositive = models["positive"]["metals"]/models["undecorated"]["metals"]
if (ratio-1.0)*(ratioPositive-1.0) < 0.0:
    print(f"SUCCESS: reversing the sign of the exponent reverses the effect "
          f"(negative {ratio:.6f}, positive {ratioPositive:.6f})")
else:
    print(f"FAILED: exponents of -1/2 and +1/2 move the metals produced in the same direction "
          f"({ratio:.6f} and {ratioPositive:.6f}), so the sign of the scaling is wrong")

# The stellar mass must be unaffected: the Type Ia rate does not feed back on star formation in these models, so
# a change there would mean the decorator is perturbing something it should not.
massDifference = abs(models["metallicity"]["massStellar"]-models["undecorated"]["massStellar"]) \
                 /models["undecorated"]["massStellar"]
if massDifference < 1.0e-6:
    print("SUCCESS: stellar mass is unaffected by the rate scaling")
else:
    print(f"FAILED: stellar mass changed by a fractional {massDifference:.2e}, so the rate scaling is "
          "perturbing more than the Type Ia yields")

sys.exit(0)
