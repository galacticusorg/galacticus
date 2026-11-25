#!/usr/bin/env python3
"""Test the `polarAngles` option of the `densityProfile` node property extractor.

With `polarAngles` the extractor evaluates the density at each radius at each of a list of
polar angles, measured from the symmetry axis of the component, writing one column per
(radius, angle) pair with the radius varying fastest. Two models are run, identical except
that one lists polar angles and the other does not. The checks are:

  * the columns are laid out as documented, with the polar angle of each as its column value;
  * the 90 degree columns are identical to the output of the model without `polarAngles`,
    which evaluates the density in the equatorial plane only;
  * the density of the exponential disk matches the analytic form,
    rho(R,z) = M/(4 pi r_s^2 z_s) exp(-R/r_s) sech^2(z/z_s), with R = r sin(theta),
    z = r cos(theta), and z_s = 0.137 r_s (the default scale height of a dimensionless disk);
  * the density of the (spherical) Hernquist spheroid does not depend on the polar angle;
  * angles symmetric about the equatorial plane (80 and 100 degrees) give the same density.

Andrew Benson (21-November-2026)
"""

import os
import subprocess
import sys

import h5py
import numpy as np

pathOutputDirectory = os.path.abspath("outputs/densityProfilePolarAngles")
subprocess.run(f"mkdir -p {pathOutputDirectory}", shell=True)

failures = []


def check(condition, description):
    print(("SUCCESS: " if condition else "FAILED: ") + description)
    if not condition:
        failures.append(description)


def runModel(name, pathParameters):
    """Run a model, returning whether it completed."""
    pathLog = os.path.join(pathOutputDirectory, f"{name}.log")
    print(f"Running model '{name}'...")
    with open(pathLog, "w") as log:
        status = subprocess.run(f"cd ..; ./Galacticus.exe {pathParameters}", stdout=log, stderr=log, shell=True)
    print(f"...done ({status.returncode})")
    if status.returncode != 0:
        subprocess.run(f"tail -n 30 {pathLog}", shell=True)
    return status.returncode == 0


def read(fileName, names):
    """Read the given node datasets from the (single) output, sorted by tree and node index.

    Trees are evolved in parallel, so their order in the output differs from one run to the next. Node indices are
    unique only within a tree, so nodes are identified by the pair (tree index, node index).
    """
    with h5py.File(fileName, "r") as model:
        output = model["Outputs"][list(model["Outputs"].keys())[0]]
        nodes = output["nodeData"]
        data = {name: nodes[name][:] for name in names}
        data["treeIndex"] = np.repeat(output["mergerTreeIndex"][:], output["mergerTreeCount"][:])
    order = np.lexsort((data["nodeIndex"], data["treeIndex"]))
    return {name: value[order] for name, value in data.items()}


completedAngles = runModel("polarAngles", "testSuite/parameters/densityProfilePolarAngles.xml")
completedEquatorial = runModel("equatorial", "testSuite/parameters/densityProfilePolarAnglesEquatorial.xml")
check(completedAngles, "model with polar angles runs to completion")
check(completedEquatorial, "model without polar angles runs to completion")
if not (completedAngles and completedEquatorial):
    print("FAILED: a model did not complete - skipping output checks")
    # Always exit with status 0 - failure is signaled by "FAILED" in the output above.
    sys.exit(0)

# The radius specifiers and polar angles of the models.
countRadii = 4
angles = np.array([90.0, 80.0, 45.0, 0.0, 100.0])
countAngles = angles.size

properties = ["nodeIndex", "densityProfile", "densityProfileRadius", "diskMassStellar", "diskRadius"]
model = read(os.path.join(pathOutputDirectory, "polarAngles.hdf5"), properties)
equatorial = read(os.path.join(pathOutputDirectory, "equatorial.hdf5"), properties)
with h5py.File(os.path.join(pathOutputDirectory, "polarAngles.hdf5"), "r") as file:
    nodes = file["Outputs"][list(file["Outputs"].keys())[0]]["nodeData"]
    columns = [c.decode() for c in nodes["densityProfileColumns"][:]]
    values = nodes["densityProfileColumnValues"][:]

# Layout: one column per (radius, angle), radius varying fastest, labeled with both.
check(model["densityProfile"].shape[1] == countRadii*countAngles, "one column per radius and polar angle")
check(np.array_equal(values, np.repeat(angles, countRadii)), "column values are the polar angles, radius varying fastest")
check(
    all(columns[i+j*countRadii].endswith(f"θ={angles[j]:.4f}°") for i in range(countRadii) for j in range(countAngles)),
    "column descriptions name the polar angle",
)
check(
    all(columns[i] .split()[0] == columns[i+j*countRadii].split()[0] for i in range(countRadii) for j in range(countAngles)),
    "column descriptions name the same radius specifier at each angle",
)

density = model["densityProfile"].reshape(-1, countAngles, countRadii)
radius = model["densityProfileRadius"].reshape(-1, countAngles, countRadii)

# The equatorial columns must be identical to those of the model without polar angles.
check(
    np.array_equal(model["treeIndex"], equatorial["treeIndex"]) and np.array_equal(model["nodeIndex"], equatorial["nodeIndex"]),
    "both models contain the same nodes",
)
check(
    np.array_equal(density[:, 0, :], equatorial["densityProfile"]),
    "90 degree columns equal the output without polarAngles",
)
check(
    np.array_equal(radius[:, 0, :], equatorial["densityProfileRadius"]),
    "90 degree radii equal the output without polarAngles",
)
# The radius does not depend on the polar angle.
check(np.all(radius == radius[:, :1, :]), "radii are the same at every polar angle")

# The exponential disk: compare with the analytic density.
hasDisk = model["diskMassStellar"] > 0.0
check(np.count_nonzero(hasDisk) >= 5, "several nodes have a stellar disk")
radiusScale = model["diskRadius"][hasDisk, None, None]
heightScale = 0.137*radiusScale
theta = np.deg2rad(angles)[None, :, None]
r = radius[hasDisk, :, :3]
R, z = r*np.sin(theta), r*np.cos(theta)
densityAnalytic = (
    model["diskMassStellar"][hasDisk, None, None]/(4.0*np.pi*radiusScale**2*heightScale)
    * np.exp(-R/radiusScale)/np.cosh(z/heightScale)**2
)
errorDisk = np.max(np.abs(density[hasDisk, :, :3]/densityAnalytic-1.0))
print(f"   maximum fractional error in the disk density: {errorDisk:.3e}")
check(errorDisk < 1.0e-6, "disk density matches the analytic exponential disk at all polar angles")
check(
    np.all(density[hasDisk, 3, :3] < density[hasDisk, 0, :3]),
    "disk density on the axis is below that in the equatorial plane",
)

# Symmetry about the equatorial plane: 80 and 100 degrees.
check(np.allclose(density[:, 1, :], density[:, 4, :], rtol=1.0e-12, atol=0.0), "density is symmetric about the equatorial plane")

# The spheroid is spherical: its density is the same at every polar angle.
densitySpheroid = density[:, :, 3]
hasSpheroid = densitySpheroid[:, 0] > 0.0
check(np.any(hasSpheroid), "some node has a stellar spheroid")
check(
    np.allclose(densitySpheroid[hasSpheroid], densitySpheroid[hasSpheroid, :1], rtol=1.0e-12, atol=0.0),
    "spheroid density does not depend on polar angle",
)

print()
if failures:
    print("FAILED: " + str(len(failures)) + " check(s) failed:")
    for failure in failures:
        print("   " + failure)
    # Always exit with status 0 - failure is signaled by "FAILED" in the output above.
    sys.exit(0)
print("SUCCESS: all checks passed")
