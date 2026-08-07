#!/usr/bin/env python3
"""Test radius specifiers which evaluate to zero.

A radius specifier can resolve to zero - `diskRadius:all:all:1.0` in a node whose
disk is absent or empty, for example. The radius is perfectly well defined, but
not every property exists there: the density of an NFW halo does not, while the
mass enclosed by a sphere or cylinder of zero radius is simply zero.

Two models are run:

  * one which asks for the density at such a radius, checking that it aborts with
    an error naming the offending specifier rather than raising a floating point
    exception;
  * one which sets `zeroRadiusIsFatal` to false, checking that it completes, that
    the sentinel appears in the density and projected density exactly where the
    radius is zero (while the radius itself is still reported as zero), that the
    enclosed and projected masses are defined and agree there, and that the
    density and projected density of the beta-profile hot halo - which has a
    finite central density, and so also a finite projected density along a line
    of sight through its center - are still evaluated at exactly zero radius. The
    central projected density is checked against the projected mass within a
    small, non-zero radius, which must approach the product of that projected
    density and the area of the cylinder.

Andrew Benson (03-August-2026)
"""

import os
import subprocess
import sys

import h5py
import numpy as np

# The `radiusUndefined` sentinel from `Galactic_Structure_Radii_Definitions`.
radiusUndefined = -np.finfo('f8').max

pathOutputDirectory = os.path.abspath("outputs/radiusSpecifiersZero")
subprocess.run(f"mkdir -p {pathOutputDirectory}", shell=True)

failures = []


def runModel(name, pathParameters):
    """Run a model, returning its exit status and the text of its log."""
    pathLog = os.path.join(pathOutputDirectory, f"{name}.log")
    print(f"Running model '{name}'...")
    with open(pathLog, "w") as log:
        status = subprocess.run(
            f"cd ..; ./Galacticus.exe {pathParameters}",
            stdout=log, stderr=log, shell=True,
        )
    print(f"...done ({status.returncode})")
    with open(pathLog) as log:
        return status.returncode, log.read()


def check(condition, description):
    print(("SUCCESS: " if condition else "FAILED: ") + description)
    if not condition:
        failures.append(description)


def gather(fileName, name, countRadii):
    """Concatenate the per-output node data of the given property across all outputs."""
    values, radii = [], []
    with h5py.File(fileName, "r") as model:
        for output in sorted(model["Outputs"].keys()):
            nodes = model["Outputs"][output]["nodeData"]
            values.append(nodes[name          ][:].reshape(-1, countRadii))
            radii .append(nodes[name+"Radius" ][:].reshape(-1, countRadii))
    return np.concatenate(values), np.concatenate(radii)


# Model 1: the error must be reported, and must name the specifier.
status, log = runModel("fatal", "testSuite/parameters/radiusSpecifiersZeroFatal.xml")
check(status != 0, "model 'fatal' exits with non-zero status")
check(
    "floating point exception" not in log.lower(),
    "model 'fatal' does not raise a floating point exception",
)
check(
    "Radius specifier evaluates to a radius of zero" in log,
    "model 'fatal' reports a zero radius error",
)
check("diskRadius" in log, "model 'fatal' names the offending radius specifier")
check("densityProfile" in log, "model 'fatal' names the property extractor")
if failures:
    # Show the tail of the log to aid diagnosis of any of the above failures.
    subprocess.run(f"tail -n 30 {pathOutputDirectory}/fatal.log", shell=True)

# Model 2: with the error suppressed the model must run to completion.
status, log = runModel("tolerated", "testSuite/parameters/radiusSpecifiersZeroTolerated.xml")
check(status == 0, "model 'tolerated' runs to completion")
check(
    "floating point exception" not in log.lower(),
    "model 'tolerated' does not raise a floating point exception",
)
if status != 0:
    subprocess.run(f"tail -n 30 {pathOutputDirectory}/tolerated.log", shell=True)
    print("FAILED: model 'tolerated' did not complete - skipping output checks")
    # Always exit with status 0 - failure is signaled by "FAILED" in the output above.
    sys.exit(0)

fileName = os.path.join(pathOutputDirectory, "tolerated.hdf5")
density        , radiiDensity         = gather(fileName, "densityProfile" , 3)
densityProjected, radiiDensityProjected = gather(fileName, "projectedDensity", 4)
mass           , radiiMass            = gather(fileName, "massProfile"   , 2)
massProjected  , radiiMassProjected   = gather(fileName, "projectedMass" , 3)

# The first two radii of the `densityProfile` extractor are the disk radius and the virial radius;
# the third is the fixed, zero radius in the hot halo.
isZero = radiiDensity[:, 0] == 0.0
check(np.any(isZero), "some node has a zero disk radius")
check(
    np.all(density[isZero, 0] == radiusUndefined),
    "density is the sentinel where the radius is zero",
)
check(
    np.all(density[~isZero, 0] > 0.0),
    "density is evaluated where the disk radius is non-zero",
)
check(
    np.all(density[:, 1] > 0.0),
    "density is evaluated at the virial radius in every node",
)
check(
    np.all(densityProjected[radiiDensityProjected[:, 0] == 0.0, 0] == radiusUndefined),
    "projected density is the sentinel where the radius is zero in a cusped profile",
)

# Enclosed and projected masses exist at zero radius: both reduce to the mass of any central point
# mass - here the black hole seed - and so must agree with each other, and must never be the
# sentinel.
massZero          = mass         [radiiMass         [:, 0] == 0.0, 0]
massProjectedZero = massProjected[radiiMassProjected[:, 0] == 0.0, 0]
check(massZero.size > 0, "some node has a zero disk radius in the mass profile")
check(
    np.all(massZero >= 0.0),
    "enclosed mass is defined where the radius is zero",
)
check(
    np.all(massProjectedZero == massZero),
    "projected mass at zero radius equals the enclosed mass of the central point mass",
)

# The beta-profile hot halo has a finite central density, so evaluating it at exactly zero radius
# is well defined and must not be rejected - even though this extractor leaves `zeroRadiusIsFatal`
# at its default of true.
check(
    np.all(radiiDensity[:, 2] == 0.0),
    "the fixed hot halo radius is zero",
)
check(
    np.all(density[:, 2] >= 0.0) and np.any(density[:, 2] > 0.0),
    "the central density of the hot halo is evaluated at zero radius",
)

# The line of sight integral through the center of the hot halo converges for the same reason - the
# central logarithmic slope of a beta-profile is zero, which exceeds the -1 that the integral requires
# - and so must also be evaluated, rather than rejected as it is for the NFW halo above.
densityProjectedCentral = densityProjected[:, 2]
check(
    np.all(radiiDensityProjected[:, 2] == 0.0),
    "the fixed hot halo radius is zero in the projected density",
)
check(
    np.all(densityProjectedCentral != radiusUndefined),
    "the central projected density of the hot halo is evaluated at zero radius",
)
check(
    np.all(densityProjectedCentral >= 0.0) and np.any(densityProjectedCentral > 0.0),
    "the central projected density of the hot halo is positive",
)

# The projected mass within a cylinder of radius R approaches pi R^2 times the central projected
# density as R goes to zero. Evaluated at a radius far inside the core radius of the beta-profile
# (0.3 of the virial radius), this is an independent check on the value of the central projected
# density - the two are computed by different extractors, with only the density profile in common.
# Nodes in which the hot halo is empty, or in which it has been stripped back to inside this radius
# (so that the projected density has already fallen away there), are excluded.
inCore = (densityProjectedCentral > 0.0) & (densityProjected[:, 3] > 0.9*densityProjectedCentral)
check(np.any(inCore), "some node has a hot halo resolved well inside its core radius")
densityProjectedImplied = (
    massProjected[inCore, 2]/np.pi/radiiMassProjected[inCore, 2]**2
)
check(
    np.all(np.abs(densityProjectedImplied/densityProjectedCentral[inCore]-1.0) < 5.0e-2),
    "the central projected density of the hot halo matches the projected mass within a small radius",
)

# The central projected density of the beta-profile hot halo is also known analytically. With
# beta = 2/3 (the default), rho(r) = rho_0/[1+(r/r_c)^2], truncated at the outer radius r_o, and
# normalized such that the mass within that radius is the hot halo mass, M:
#
#   Sigma(0) = 2 rho_0 r_c atan(x) = M atan(x)/[2 pi r_c^2 (x - atan x)],   x = r_o/r_c,
#
# where the core radius, r_c, is 0.3 of the virial radius here (see `hotHaloMassDistributionCoreRadius`
# above), and the virial radius is reported as the second radius of the `densityProfile` extractor.
with h5py.File(fileName, "r") as model:
    massHotHalo, radiusOuterHotHalo = [], []
    for output in sorted(model["Outputs"].keys()):
        nodes = model["Outputs"][output]["nodeData"]
        massHotHalo       .append(nodes["hotHaloMass"       ][:])
        radiusOuterHotHalo.append(nodes["hotHaloOuterRadius"][:])
    massHotHalo        = np.concatenate(massHotHalo       )
    radiusOuterHotHalo = np.concatenate(radiusOuterHotHalo)
radiusCore = 0.3*radiiDensity[:, 1]
x          = radiusOuterHotHalo/radiusCore
densityProjectedAnalytic = (
    massHotHalo*np.arctan(x)/(2.0*np.pi*radiusCore**2*(x-np.arctan(x)))
)
check(
    np.all(np.abs(densityProjectedAnalytic[inCore]/densityProjectedCentral[inCore]-1.0) < 1.0e-2),
    "the central projected density of the hot halo matches the analytic beta-profile result",
)

print()
if failures:
    print("FAILED: " + str(len(failures)) + " check(s) failed:")
    for failure in failures:
        print("   " + failure)
    # Always exit with status 0 - failure is signaled by "FAILED" in the output above.
    sys.exit(0)
print("SUCCESS: all checks passed")
