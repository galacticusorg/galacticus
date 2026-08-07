#!/usr/bin/env python3
"""Test radius specifiers expressed as fractions of fuzzy dark matter soliton radii.

Radii given as `solitonRadiusCore` or `solitonRadiusSoliton` differ from every
other radius type in that they can be *undefined*: the model may contain no
soliton-forming dark matter profile at all (so nothing creates the
meta-properties in which the radii are stored), or a halo may simply have grown
no soliton. In either case the radius must be reported as the `radiusUndefined`
sentinel rather than as zero - evaluating a profile which diverges at the center
at a radius of zero raises a floating point exception.

Two models are run:

  * one which grows solitons, checking that the resolved radii are the requested
    multiples of the stored radii, and that the sentinel appears exactly where
    the stored radius is non-positive;
  * one with an NFW profile and no soliton anywhere, checking that the run
    completes (rather than aborting because the meta-properties have no creator,
    or trapping on the divergent central density) and that every soliton radius
    is the sentinel.

Andrew Benson (31-July-2026)
"""

import os
import subprocess
import sys

import h5py
import numpy as np

# The `radiusUndefined` sentinel from `Galactic_Structure_Radii_Definitions`.
radiusUndefined = -np.finfo('f8').max

pathOutputDirectory = os.path.abspath("outputs/radiusSpecifiersSoliton")
subprocess.run(f"mkdir -p {pathOutputDirectory}", shell=True)

failures = []


def runModel(name, pathParameters):
    """Run a model, failing the test if it does not complete cleanly."""
    pathLog = os.path.join(pathOutputDirectory, f"{name}.log")
    print(f"Running model '{name}'...")
    with open(pathLog, "w") as log:
        status = subprocess.run(
            f"cd ..; ./Galacticus.exe {pathParameters}",
            stdout=log, stderr=log, shell=True,
        )
    print(f"...done ({status.returncode})")
    if status.returncode != 0:
        failures.append(f"model '{name}' failed to run")
        subprocess.run(f"tail -n 40 {pathLog}", shell=True)
        return False
    grep = subprocess.run(
        "grep -q -i -e fatal -e aborted -e 'floating point'"
        f" -e 'Galacticus experienced an error in the GSL library' {pathLog}",
        shell=True,
    )
    if grep.returncode == 0:
        failures.append(f"model '{name}' reported errors")
        subprocess.run(f"tail -n 40 {pathLog}", shell=True)
        return False
    return True


def check(condition, description):
    print(("SUCCESS: " if condition else "FAILED: ") + description)
    if not condition:
        failures.append(description)


def gather(fileName, countRadii):
    """Concatenate the per-output node data of interest across all outputs."""
    radii, masses, radiiCore, radiiSoliton = [], [], [], []
    with h5py.File(fileName, "r") as model:
        for output in sorted(model["Outputs"].keys()):
            nodes = model["Outputs"][output]["nodeData"]
            radii .append(nodes["massProfileRadius"][:].reshape(-1, countRadii))
            masses.append(nodes["massProfile"      ][:].reshape(-1, countRadii))
            if "solitonRadiusCore" in nodes:
                radiiCore   .append(np.atleast_1d(nodes["solitonRadiusCore"   ][:]))
                radiiSoliton.append(np.atleast_1d(nodes["solitonRadiusSoliton"][:]))
    return (
        np.vstack(radii),
        np.vstack(masses),
        np.concatenate(radiiCore   ) if radiiCore    else None,
        np.concatenate(radiiSoliton) if radiiSoliton else None,
    )


# --- A model which grows solitons -------------------------------------------
if runModel("soliton", "testSuite/parameters/radiusSpecifiersSoliton.xml"):
    radii, masses, radiiCore, radiiSoliton = gather(
        os.path.join(pathOutputDirectory, "soliton.hdf5"), 5)
    # Columns, in the order requested by the parameter file.
    core1, core2, soliton1, soliton2, virial = 0, 1, 2, 3, 4
    defined = radiiCore > 0.0
    check(bool(np.all(np.isclose(radii[defined, core1], radiiCore[defined], rtol=1.0e-14))),
          "a radius of 1.0 core radii equals the stored core radius")
    check(bool(np.all(np.isclose(radii[defined, core2], 2.0 * radiiCore[defined], rtol=1.0e-14))),
          "a radius of 2.0 core radii equals twice the stored core radius")
    definedSoliton = radiiSoliton > 0.0
    check(bool(np.all(np.isclose(radii[definedSoliton, soliton1],
                                 radiiSoliton[definedSoliton], rtol=1.0e-14))),
          "a radius of 1.0 soliton radii equals the stored soliton radius")
    check(bool(np.all(np.isclose(radii[definedSoliton, soliton2],
                                 0.5 * radiiSoliton[definedSoliton], rtol=1.0e-14))),
          "a radius of 0.5 soliton radii equals half the stored soliton radius")
    check(bool(np.array_equal(radii[:, core1   ] == radiusUndefined, radiiCore    <= 0.0)),
          "the core radius is undefined exactly where no soliton core exists")
    check(bool(np.array_equal(radii[:, soliton1] == radiusUndefined, radiiSoliton <= 0.0)),
          "the soliton radius is undefined exactly where no soliton exists")
    check(bool(np.all(masses[radii == radiusUndefined] == radiusUndefined)),
          "the enclosed mass is undefined wherever the radius is undefined")
    check(bool(np.all(radii[:, virial] != radiusUndefined)),
          "virial radii in the same extractor remain defined")

# --- A model with no soliton-forming profile at all --------------------------
if runModel("noSoliton", "testSuite/parameters/radiusSpecifiersNoSoliton.xml"):
    radii, masses, _, _ = gather(
        os.path.join(pathOutputDirectory, "noSoliton.hdf5"), 5)
    check(bool(np.all(radii[:, 0:4] == radiusUndefined)),
          "all soliton radii are undefined when no class creates them")
    check(bool(np.all(masses[:, 0:4] == radiusUndefined)),
          "all masses at soliton radii are undefined when no class creates them")
    check(bool(np.all(radii[:, 4] != radiusUndefined)),
          "virial radii remain defined when soliton radii are not")

if failures:
    print(f"\nFAILED: {len(failures)} check(s) failed")
    for failure in failures:
        print(f"   {failure}")
    # Always exit with status 0 - failure is signaled by "FAILED" in the output above.
    sys.exit(0)
print("\nSUCCESS: all checks passed")
sys.exit(0)
