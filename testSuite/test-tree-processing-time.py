#!/usr/bin/env python3
import subprocess
import sys
import os
import h5py
import numpy as np

# Test run-time estimation calibration: the `mergerTreeOperatorTreeProcessingTimer` operator should fit and write a cost model
# (coefficients, residual scatter, and range, versus both tree mass and node count) into the output file, and the `file`
# `metaTreeProcessingTime` class should be able to read that cost model directly from the HDF5 file in a subsequent run.
# Andrew Benson (10-July-2026)

# Ensure the output directory exists.
os.makedirs("outputs", exist_ok=True)

# Run the producer model (tree processing timer operator active).
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/treeProcessingTime.xml", shell=True)
if status.returncode != 0:
    print("FAILED: producer model run")
    sys.exit(0)
print("SUCCESS: producer model run")

# Verify that the fit results were written, with sane shapes and values.
success = True
with h5py.File("outputs/treeProcessingTime.hdf5", "r") as model:
    if "metaData/treeTiming" not in model:
        print("FAILED: no metaData/treeTiming group written")
        sys.exit(0)
    timing = model["metaData/treeTiming"]
    # Check that all expected fit datasets are present with the correct sizes.
    expected = {
        "fitCoefficientMass"      : 3,
        "fitResidualMass"         : 1,
        "fitRangeMass"            : 2,
        "fitCoefficientCountNodes": 3,
        "fitResidualCountNodes"   : 1,
        "fitRangeCountNodes"      : 2,
    }
    for name, size in expected.items():
        if name not in timing:
            print(f"FAILED: fit dataset '{name}' missing")
            success = False
        elif timing[name].shape[0] != size:
            print(f"FAILED: fit dataset '{name}' has size {timing[name].shape[0]}, expected {size}")
            success = False
    if success:
        coefficients = timing["fitCoefficientMass"][:]
        residual     = timing["fitResidualMass"   ][0]
        rangeMass    = timing["fitRangeMass"       ][:]
        # Coefficients and residual must be finite, the residual non-negative, and the mass range ordered and positive.
        if not np.all(np.isfinite(coefficients)):
            print("FAILED: mass fit coefficients are not finite")
            success = False
        if not (np.isfinite(residual) and residual >= 0.0):
            print(f"FAILED: mass fit residual is invalid: {residual}")
            success = False
        if not (0.0 < rangeMass[0] <= rangeMass[1]):
            print(f"FAILED: mass fit range is invalid: {rangeMass}")
            success = False
        # Independently reproduce the fit from the raw timing data and confirm the stored coefficients reproduce it.
        mass        = timing["treeMass"     ][:].astype("float64")
        timeProcess = timing["timeConstruct"][:].astype("float64") + timing["timeEvolve"][:].astype("float64")
        select      = (mass > 0.0) & (timeProcess > 0.0)
        if np.count_nonzero(select) >= 3:
            x        = np.log10(mass[select])
            y        = np.log10(timeProcess[select])
            model_y  = coefficients[0] + coefficients[1] * x + coefficients[2] * x**2
            residualCheck = np.sqrt(np.sum((y - model_y)**2) / (x.size - 3))
            if abs(residualCheck - residual) > 1.0e-3 * max(residualCheck, 1.0e-3):
                print(f"FAILED: stored residual {residual} does not match recomputed residual {residualCheck}")
                success = False
if success:
    print("SUCCESS: tree processing time fit written and validated")

# Run the consumer model, which reads the cost model directly from the producer's HDF5 output file. Success here validates the
# HDF5 read path of the `file` metaTreeProcessingTime class.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/treeProcessingTimeConsume.xml", shell=True)
if status.returncode != 0:
    print("FAILED: consumer model run (reading cost model from HDF5)")
    sys.exit(0)
print("SUCCESS: consumer model run (cost model read from HDF5)")

# Run a model with task-level progress reporting enabled and a task-scope cost model, and confirm that the start-of-run estimate,
# the throttled progress reports, and the end-of-run summary are all emitted.
status = subprocess.run(
    "cd ..; ./Galacticus.exe testSuite/parameters/treeProcessingTimeProgress.xml > testSuite/outputs/treeProcessingTimeProgress.log 2>&1",
    shell=True
)
if status.returncode != 0:
    print("FAILED: progress-reporting model run")
    sys.exit(0)
with open("outputs/treeProcessingTimeProgress.log") as logFile:
    progressLog = logFile.read()
checks = {
    "start-of-run run-time estimate": "Run-time estimate:",
    "throttled progress report"     : "Progress:",
    "end-of-run run-time summary"   : "Run-time summary:",
}
progressSuccess = True
for label, marker in checks.items():
    if marker in progressLog:
        print(f"SUCCESS: {label} emitted")
    else:
        print(f"FAILED: {label} not emitted (marker '{marker}' absent)")
        progressSuccess = False
# The estimate should reference a positive tree count, and the final progress report should reach 100%.
if "Run-time estimate:" in progressLog and " trees," in progressLog and "100.0% complete" not in progressLog:
    print("FAILED: progress never reached 100% complete")
    progressSuccess = False
if progressSuccess:
    print("SUCCESS: task-level progress and run-time estimation reporting")
