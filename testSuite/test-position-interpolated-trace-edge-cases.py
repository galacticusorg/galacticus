#!/usr/bin/env python3
"""Test the positionInterpolated node operator on two kinds of branch common in SUBFIND merger trees.

A model with a lightcone is run on two forests from the milli-Millennium simulation, each of which once made the model fail:

  * forest 6000000 contains a branch that begins as a subhalo with no progenitors and no position history (node 6000387).
    It is followed in its host, and the host's parent is a clone of the host at the same time, so the trace of its position
    once held two steps at the same time, and the model stopped with "duplicated time in position trace";
  * forest 1569000000 is a single node. Its trace has a single step, and so no interval to interpolate across: its position
    came out as NaN (read from beyond the end of an empty array of interpolation coefficients), and the lightcone crossing
    calculation then raised a floating point exception.

The model must now run to completion; the traces of both nodes must be reported (so that the cases were exercised), with
no repeated time in the trace of the first, and the single node's position must be left as given rather than interpolated;
and all positions in the lightcone must be finite.

Andrew Benson (21-November-2026)
"""

import re
import subprocess
import sys

import h5py
import numpy as np

subprocess.run("mkdir -p outputs", shell=True)
failures = []


def check(condition, description):
    print(("SUCCESS: " if condition else "FAILED: ") + description)
    if not condition:
        failures.append(description)


print("Running model...")
with open("outputs/positionTraceEdgeCases.log", "w") as log:
    status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/positionTraceEdgeCases.xml", stdout=log,
                            stderr=log, shell=True)
print(f"...done ({status.returncode})")
with open("outputs/positionTraceEdgeCases.log", errors="replace") as log:
    text = log.read()
check(status.returncode == 0, "model runs to completion")
check("duplicated time" not in text, "no duplicated time in any position trace")
check("floating point exception" not in text.lower(), "no floating point exception")
if status.returncode != 0:
    subprocess.run("tail -n 30 outputs/positionTraceEdgeCases.log", shell=True)
    print("FAILED: model did not complete - skipping output checks")
    # Always exit with status 0 - failure is signaled by "FAILED" in the output above.
    sys.exit(0)

# Both edge cases were traced.
check("Tracing node: 6000387" in text, "the branch beginning as a subhalo was traced")
check("Tracing node: 1569000000" in text, "the single-node tree was traced")
check("No interpolation intervals - position unchanged" in text, "the single node's position is not interpolated")
# The trace report of the first lists each step (node, host, time, satellite status, history index): no time may repeat.
lines = text.splitlines()
start = [i for i, line in enumerate(lines) if "Trace report for node: 6000387" in line]
check(len(start) > 0, "the trace of the branch beginning as a subhalo was reported")
if start:
    step = re.compile(r"^\s*(?:\d+:)?\s+(\d+)\s+(\d+)\s+([0-9.]+)\s+[TF]\s+-?\d+\s*$")
    times = []
    # Steps are interleaved with details of each interpolation; the report ends where the next one (or the next tree) begins.
    for line in lines[start[0] + 1:]:
        if any(marker in line for marker in ("Trace report for node", "Tracing node", "Evolving tree", "Finished tree")):
            break
        match = step.match(line)
        if match:
            times.append(float(match.group(3)))
    check(len(times) >= 2 and len(times) == len(set(times)), f"its trace has {len(times)} steps, at distinct times")

with h5py.File("outputs/positionTraceEdgeCases.hdf5", "r") as model:
    nodes = model["Lightcone/Output1/nodeData"]
    position = np.stack([nodes[f"lightconePosition{axis}"][:] for axis in "XYZ"], axis=1)
check(position.shape[0] > 0, f"the lightcone contains galaxies ({position.shape[0]})")
check(bool(np.all(np.isfinite(position))), "all lightcone positions are finite")

print()
if failures:
    print("FAILED: " + str(len(failures)) + " check(s) failed:")
    for failure in failures:
        print("   " + failure)
    # Always exit with status 0 - failure is signaled by "FAILED" in the output above.
    sys.exit(0)
print("SUCCESS: all checks passed")
