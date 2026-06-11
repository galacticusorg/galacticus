#!/usr/bin/env python3
"""
End-to-end integration test of the parametric self-interacting dark matter (SIDM) model of Yang et al. (2024). A small set of
merger trees is built and evolved, applying the reference parametric SIDM changes on top of the reference cold dark matter
model. As we have no external reference trajectory to compare against, the test verifies that the model runs and that the
evolved gravothermal state is internally consistent, plus the limiting behavior that the model reduces to CDM when the cross
section is taken to zero.
"""
import subprocess
import sys
import h5py
import numpy as np

# SIDM parametric meta-property datasets written by the SIDMParametric property extractor.
tauName  =  "darkMatterProfileSIDMParametricTau"
sidmNames = [
    "darkMatterProfileSIDMParametricVelocityMaximum",
    "darkMatterProfileSIDMParametricRadiusMaximum"  ,
    "darkMatterProfileSIDMParametricDensityScale"   ,
    "darkMatterProfileSIDMParametricRadiusScale"    ,
    "darkMatterProfileSIDMParametricRadiusCore"     ,
]

failures = 0

def runModel(extraChanges=""):
    """Run the parametric SIDM model (base CDM model + SIDM changes [+ extra changes])."""
    command = (
        "cd ..; ./Galacticus.exe"
        " testSuite/parameters/SIDM_parametric.xml"
        " parameters/reference/darkMatterParametricSIDM.xml"
        + (" "+extraChanges if extraChanges else "")
    )
    return subprocess.run(command, shell=True).returncode

def gather(fileName, names):
    """Concatenate the named node datasets across all output times."""
    collected = {name: [] for name in names}
    with h5py.File(fileName, "r") as model:
        for outputName in model["Outputs"]:
            nodeData = model["Outputs"][outputName].get("nodeData")
            if nodeData is None:
                continue
            for name in names:
                if name in nodeData:
                    collected[name].append(nodeData[name][:])
    return {name: (np.concatenate(values) if values else np.array([])) for name, values in collected.items()}

# --- Fiducial run -----------------------------------------------------------------------------------------------------------
if runModel() != 0:
    print("FAIL: parametric SIDM model failed to run")
    sys.exit(1)
data = gather("outputs/SIDM_parametric.hdf5", [tauName]+sidmNames)
tau  = data[tauName]
if tau.size == 0:
    print("FAIL: parametric SIDM model produced no SIDM node data")
    sys.exit(1)
print("SUCCESS: parametric SIDM model ran and produced output")

# tau must lie in [0,1] (a small tolerance allows for numerical round-off).
if np.all((tau >= -1.0e-6) & (tau <= 1.0+1.0e-6)):
    print("SUCCESS: gravothermal time tau lies within [0,1]")
else:
    print(f"FAIL: gravothermal time tau outside [0,1] (min={tau.min():.3e}, max={tau.max():.3e})")
    failures += 1

# Every SIDM quantity must be finite (no NaNs or infinities).
if all(np.all(np.isfinite(data[name])) for name in [tauName]+sidmNames):
    print("SUCCESS: all SIDM quantities are finite")
else:
    nonFinite = [name for name in [tauName]+sidmNames if not np.all(np.isfinite(data[name]))]
    print(f"FAIL: non-finite values in {', '.join(nonFinite)}")
    failures += 1

# The profile-defining quantities (scale density and scale radius, which the SIDM profile is actually built from) must be
# strictly positive, and the core radius must be non-negative (it is zero before a core forms).
densityScale = data["darkMatterProfileSIDMParametricDensityScale"]
radiusScale  = data["darkMatterProfileSIDMParametricRadiusScale" ]
radiusCore   = data["darkMatterProfileSIDMParametricRadiusCore"  ]
if np.all(densityScale > 0.0) and np.all(radiusScale > 0.0) and np.all(radiusCore >= 0.0):
    print("SUCCESS: SIDM scale density and scale radius are positive and the core radius is non-negative")
else:
    print(f"FAIL: unphysical SIDM profile parameters (min density scale={densityScale.min():.3e}, "
          f"min scale radius={radiusScale.min():.3e}, min core radius={radiusCore.min():.3e})")
    failures += 1

# The test should actually exercise the gravothermal evolution (otherwise it is vacuous). The velocityMaximumSIDM dataset is
# the SIDM-CDM difference, so a non-zero shift confirms the maximum-velocity evolution is exercised, not just tau.
if tau.max() > 0.0:
    print(f"SUCCESS: gravothermal evolution occurred (maximum tau = {tau.max():.3g})")
else:
    print("FAIL: no gravothermal evolution occurred (maximum tau = 0); the test is vacuous")
    failures += 1
velocityMaximumShift = data["darkMatterProfileSIDMParametricVelocityMaximum"]
if velocityMaximumShift.max() > 0.0:
    print(f"SUCCESS: the SIDM evolution shifts the maximum circular velocity (max SIDM-CDM shift = {velocityMaximumShift.max():.3g} km/s)")
else:
    print("FAIL: the SIDM-CDM maximum-velocity shift is zero everywhere; the velocity evolution is not exercised")
    failures += 1

# --- No-interaction run (sigma0 -> 0): the model must reduce to CDM ----------------------------------------------------------
if runModel("testSuite/parameters/SIDM_parametric_noInteraction.xml") != 0:
    print("FAIL: no-interaction parametric SIDM model failed to run")
    failures += 1
else:
    noInteraction        = gather("outputs/SIDM_parametric_noInteraction.hdf5",
                                  [tauName,
                                   "darkMatterProfileSIDMParametricVelocityMaximum",
                                   "darkMatterProfileSIDMParametricRadiusMaximum"  ])
    tauZero              = noInteraction[tauName]
    velocityMaximumShift = noInteraction["darkMatterProfileSIDMParametricVelocityMaximum"]
    radiusMaximumShift   = noInteraction["darkMatterProfileSIDMParametricRadiusMaximum"  ]
    # The gravothermal evolution should be suppressed.
    if tauZero.size > 0 and tauZero.max() < 1.0e-2:
        print(f"SUCCESS: with sigma0 -> 0 the gravothermal evolution is suppressed (maximum tau = {tauZero.max():.3g})")
    else:
        maximum = f"{tauZero.max():.3e}" if tauZero.size > 0 else "n/a"
        print(f"FAIL: with sigma0 -> 0 the gravothermal evolution was not suppressed (maximum tau = {maximum})")
        failures += 1
    # The SIDM-CDM Vmax/Rmax differences should vanish, i.e. the SIDM maximum velocity and its radius equal the CDM values.
    if velocityMaximumShift.size > 0 and np.all(np.abs(velocityMaximumShift) < 1.0e-3) and np.all(np.abs(radiusMaximumShift) < 1.0e-8):
        print("SUCCESS: with sigma0 -> 0 the SIDM-CDM Vmax/Rmax differences vanish; the model reduces to CDM")
    else:
        print(f"FAIL: with sigma0 -> 0 the SIDM-CDM differences did not vanish (max|Vmax shift|="
              f"{np.abs(velocityMaximumShift).max():.3e} km/s, max|Rmax shift|={np.abs(radiusMaximumShift).max():.3e} Mpc)")
        failures += 1

sys.exit(1 if failures > 0 else 0)
