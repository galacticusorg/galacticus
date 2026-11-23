#!/usr/bin/env python3
"""Test that galaxy merger trees record the correct star formation rate for each galaxy at a merger.

When two galaxies merge, the `galaxyMergerTree` node operator records the properties of
the satellite and then of the central, one immediately after the other. The `blitz2006`
star formation rate surface density class memoizes factors computed from the disk of the
last galaxy for which it was evaluated, and its `unchanged` and `intervals` methods
reused those factors without checking that they belonged to the galaxy being evaluated.
The central was then given the star formation rate of the satellite - typically far
smaller, so the central appeared to stop forming stars at the instant of each minor merger.

A model using `blitz2006` is run, and for every merger:

  * the star formation rate recorded for the central must not equal that recorded for the
    satellite (they were identical, to the last bit, whenever the bug struck);
  * in minor mergers (satellite baryonic mass below 1% of the central's), which barely
    change the central, the central's recorded star formation rate must not collapse
    relative to its neighboring samples.

Andrew Benson (21-November-2026)
"""

import os
import subprocess
import sys

import h5py
import numpy as np

pathOutputDirectory = os.path.abspath("outputs/galaxyMergerTreeStarFormation")
subprocess.run(f"mkdir -p {pathOutputDirectory}", shell=True)

failures = []


def check(condition, description):
    print(("SUCCESS: " if condition else "FAILED: ") + description)
    if not condition:
        failures.append(description)


print("Running model...")
with open(os.path.join(pathOutputDirectory, "galaxyMergerTree.log"), "w") as log:
    status = subprocess.run(
        "cd ..; ./Galacticus.exe testSuite/parameters/galaxyMergerTreeStarFormation.xml",
        stdout=log, stderr=log, shell=True,
    )
print(f"...done ({status.returncode})")
if status.returncode != 0:
    subprocess.run(f"tail -n 30 {pathOutputDirectory}/galaxyMergerTree.log", shell=True)
    print("FAILED: model did not complete")
    # Always exit with status 0 - failure is signaled by "FAILED" in the output above.
    sys.exit(0)

# Each node carries the records of every galaxy that merged into it: per branch (one galaxy until it merges), a branch
# index and a count of samples; samples (times and properties) concatenated over branches; and, per merger, the branches of
# the central and satellite and the time.
mergers = []
with h5py.File(os.path.join(pathOutputDirectory, "galaxyMergerTree.hdf5"), "r") as model:
    nodes = model["Outputs/Output1/nodeData"]
    for row in range(nodes["galaxyMergerTreeCount"].shape[0]):
        counts = nodes["galaxyMergerTreeCount"][row]
        if len(counts) < 2:
            continue
        indices = nodes["galaxyMergerTreeNodeIndex"][row]
        time = nodes["galaxyMergerTreeTime"][row]
        rate = nodes["galaxyMergerTreeTotalStarFormationRate"][row]
        mass = nodes["galaxyMergerTreeMassStellarTotal"][row] + nodes["galaxyMergerTreeMassISM"][row]
        end = np.cumsum(counts)
        start = end - counts
        branch = {int(b): slice(s, e) for b, s, e in zip(indices, start, end)}
        for central, satellite, timeMerge in zip(
            nodes["galaxyMergerTreeMergeCentralIndex"][row],
            nodes["galaxyMergerTreeMergeSatelliteIndex"][row],
            nodes["galaxyMergerTreeMergeTime"][row],
        ):
            c, s = branch[int(central)], branch[int(satellite)]
            # The merger time is the central's time, at which it is recorded. (The satellite is recorded at its own time, which
            # can be slightly later.)
            k = np.nonzero(time[c] == timeMerge)[0]
            if k.size != 1:
                continue
            k = int(k[0])
            neighbors = rate[c][[k - 1, k + 1]] if 0 < k < len(time[c]) - 1 else None
            mergers.append(
                dict(
                    rateCentral=rate[c][k],
                    rateSatellite=rate[s][-1],
                    ratioMass=mass[s][-1] / mass[c][k - 1] if k > 0 and mass[c][k - 1] > 0.0 else np.inf,
                    neighbors=neighbors,
                )
            )

print(f"   {len(mergers)} mergers found")
check(len(mergers) >= 50, "the model contains many mergers")

# The central's rate must be its own, not the satellite's.
same = [m for m in mergers if m["rateSatellite"] > 0.0 and m["rateCentral"] == m["rateSatellite"]]
print(f"   {len(same)} mergers in which the central was recorded with the satellite's star formation rate")
check(len(same) == 0, "the central's star formation rate at a merger is never that of the satellite")

# Minor mergers: the central's rate must not collapse.
minor = [m for m in mergers if m["ratioMass"] < 1.0e-2 and m["neighbors"] is not None and np.all(m["neighbors"] > 0.0)]
collapsed = [m for m in minor if m["rateCentral"] < 1.0e-2 * np.min(m["neighbors"])]
print(f"   {len(collapsed)} of {len(minor)} minor mergers in which the central's star formation rate collapses")
check(len(minor) >= 20, "the model contains many minor mergers")
check(len(collapsed) == 0, "the central's star formation rate does not collapse at minor mergers")

print()
if failures:
    print("FAILED: " + str(len(failures)) + " check(s) failed:")
    for failure in failures:
        print("   " + failure)
    # Always exit with status 0 - failure is signaled by "FAILED" in the output above.
    sys.exit(0)
print("SUCCESS: all checks passed")
