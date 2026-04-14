#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Test axion transfer functions.
# Andrew Benson (ported to Python)

models = [
    {"label": "CDM",               "testPeaks": "no"},
    {"label": "AxionCAMB",         "testPeaks": "no"},
    {"label": "AxionMurgia2017",   "testPeaks": "no",  "toleranceCutOff": 0.10},
    {"label": "AxionHu2000",       "testPeaks": "yes", "toleranceCutOff": 0.160, "tolerancePeakWavenumber": 0.090, "tolerancePeakAmplitude": 0.802},
    {"label": "AxionPassaglia2022","testPeaks": "yes", "toleranceCutOff": 0.025, "tolerancePeakWavenumber": 0.017, "tolerancePeakAmplitude": 0.165},
]

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

transferFunctions = {}

for model in models:
    # Run models.
    status = subprocess.run(
        f"cd ..; ./Galacticus.exe testSuite/parameters/powerSpectrum{model['label']}.xml",
        shell=True
    )
    if status.returncode != 0:
        print(f"FAIL: failed to run model 'powerSpectrum{model['label']}.xml'")
        sys.exit(0)

    # Read data.
    with h5py.File(f"outputs/powerSpectrum{model['label']}.hdf5", "r") as galacticus:
        output = galacticus["Outputs/Output1"]
        wavenumber       = output["wavenumber"][:]
        transferFunction = output["transferFunction"][:]

    # Normalize the transfer function.
    transferFunction /= transferFunction[0]
    transferFunction  = np.abs(transferFunction)

    transferFunctions[model["label"]] = {
        "wavenumber":       wavenumber,
        "transferFunction": transferFunction,
    }

# Divide by CDM transfer function for non-CDM models.
for model in models:
    if model["label"] != "CDM":
        transferFunctions[model["label"]]["transferFunction"] /= transferFunctions["CDM"]["transferFunction"]

# Perform tests.
for model in models:
    if model["label"] == "CDM":
        continue

    label  = model["label"]
    tf     = transferFunctions[label]["transferFunction"]
    wn     = transferFunctions[label]["wavenumber"]

    # Find locations of peaks.
    peaks = []
    for i in range(1, len(wn) - 1):
        if tf[i] > tf[i-1] and tf[i] > tf[i+1]:
            peaks.append({"wavenumber": wn[i], "transferFunction": tf[i]})
    transferFunctions[label]["peaks"] = peaks

    # Test amplitudes at ~4 Mpc^-1 against AxionCAMB.
    if "toleranceCutOff" in model:
        wavenumberTarget = 4.0
        indexTarget      = np.argmin(np.abs(wn - wavenumberTarget))
        error            = abs(tf[indexTarget] - transferFunctions["AxionCAMB"]["transferFunction"][indexTarget]) / transferFunctions["AxionCAMB"]["transferFunction"][indexTarget]
        errorPercent     = f"{100.0 * error:4.1f}%"
        status_str       = "FAIL" if error > model["toleranceCutOff"] else "SUCCESS"
        print(f"{label} T(k=4 Mpc\u207b\u00b9): {errorPercent} ({status_str})")

    # Test locations and amplitudes of peaks.
    if "tolerancePeakWavenumber" in model:
        subscripts = ["\u2080", "\u2081", "\u2082", "\u2083", "\u2084", "\u2085", "\u2086", "\u2087", "\u2088", "\u2089"]
        for peak in range(2):
            axionCAMBPeaks = transferFunctions["AxionCAMB"]["peaks"]
            modelPeaks     = transferFunctions[label]["peaks"]
            errorWavenumber = (
                abs(modelPeaks[peak]["wavenumber"] - axionCAMBPeaks[peak]["wavenumber"])
                / axionCAMBPeaks[peak]["wavenumber"]
            )
            errorAmplitude = (
                abs(modelPeaks[peak]["transferFunction"] - axionCAMBPeaks[peak]["transferFunction"])
                / axionCAMBPeaks[peak]["transferFunction"]
            )
            errorWavenumberPercent = f"{100.0 * errorWavenumber:4.1f}%"
            errorAmplitudePercent  = f"{100.0 * errorAmplitude:4.1f}%"
            statusWavenumber = "FAIL" if errorWavenumber > model["tolerancePeakWavenumber"] else "SUCCESS"
            statusAmplitude  = "FAIL" if errorAmplitude  > model["tolerancePeakAmplitude"]  else "SUCCESS"
            print(f"{label} k{subscripts[peak+1]}: {errorWavenumberPercent} ({statusWavenumber})")
            print(f"{label} T{subscripts[peak+1]}: {errorAmplitudePercent} ({statusAmplitude})")
