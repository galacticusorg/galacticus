#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Run models that test that subhalo tidal track evolution by validating against the fitting function of Errani & Navarro (2021).
# Andrew Benson (ported to Python)

# Make output directory.
subprocess.run("mkdir -p outputs/", shell=True)

testCases = [
    {"label": "nonMonotonic", "gamma": 1.00, "fitMetric": 0.0301},
    {"label": "monotonic",    "gamma": 1.00, "fitMetric": 0.0062},
    {"label": "monotonic",    "gamma": 0.50, "fitMetric": 0.0056},
    {"label": "monotonic",    "gamma": 0.00, "fitMetric": 0.0280},
]

for testCase in testCases:
    label    = testCase["label"][0].upper() + testCase["label"][1:]
    gammaStr = f"{testCase['gamma']:3.1f}"
    xmlFile  = f"testSuite/parameters/tidalTracks{label}_gamma{gammaStr}.xml"
    hdf5File = f"outputs/tidalTracks{label}_gamma{gammaStr}.hdf5"

    # Run the tidal tracks model.
    status = subprocess.run(
        f"export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe {xmlFile}",
        shell=True
    )
    if status.returncode != 0:
        print(f"FAIL: tidal track model '{testCase['label']} gamma={gammaStr}' failed to run")
        sys.exit(0)

    # Collect model track data.
    tModel = []
    rModel = []
    vModel = []
    with h5py.File(hdf5File, "r") as model:
        outputs = model["Outputs"]
        for outputName in outputs.keys():
            output   = outputs[outputName]
            nodeData = output["nodeData"]
            timeArr      = nodeData["time"][:]
            isIsolated   = nodeData["nodeIsIsolated"][:]
            vmax         = nodeData["darkMatterProfileDMOVelocityMaximum"][:]
            rmax         = nodeData["darkMatterProfileDMORadiusVelocityMaximum"][:]
            selection    = np.where(isIsolated == 0)[0]
            tModel.extend(timeArr[selection])
            rModel.extend(rmax[selection])
            vModel.extend(vmax[selection])

    tModel = np.array(tModel)
    rModel = np.array(rModel)
    vModel = np.array(vModel)
    order  = np.argsort(tModel)
    r0     = rModel[order][0]
    v0     = vModel[order][0]
    rModel = rModel / r0
    vModel = vModel / v0

    # Construct N-body tidal track.
    la = np.linspace(-6.0, 0.0, 1000)
    a  = 10.0**la

    if testCase["gamma"] == 1.0:
        # Tidal track from Errani & Navarro (2021).
        alpha          = 0.40
        beta           = 0.65
        rNBody         = a
        vNBody         = 2.0**alpha * a**beta / (1.0 + a**2)**alpha
    else:
        # Tidal tracks from Penarrubia et al. (2010).
        if testCase["gamma"] == 0.0:
            mur, etar, muv, etav = -1.30, +0.05, +0.40, +0.37
        elif testCase["gamma"] == 0.5:
            mur, etar, muv, etav = -0.40, +0.27, +0.40, +0.35
        elif testCase["gamma"] == 1.5:
            mur, etar, muv, etav = +0.00, +0.48, +0.40, +0.24
        else:
            raise ValueError(f"unknown gamma: {testCase['gamma']}")
        r      = 2.0**mur * a**etar / (1.0 + a)**mur
        v      = 2.0**muv * a**etav / (1.0 + a)**muv
        sortOrd= np.argsort(r)
        rNBody = r[sortOrd]
        vNBody = v[sortOrd]

    # Compute fit metric.
    selection = np.where((vModel > 0.1) & (rModel > 0.1))[0]
    fitMetric = 0.0
    logRNBody = np.log10(rNBody)
    logVNBody = np.log10(vNBody)
    for i in selection:
        distSq    = (logRNBody - np.log10(rModel[i]))**2 + (logVNBody - np.log10(vModel[i]))**2
        fitMetric += distSq.min()

    status_str = "SUCCESS" if fitMetric < testCase["fitMetric"] else "FAILED"
    print(f"{status_str}: subhalo tidal tracks '{testCase['label']} gamma={gammaStr}'")
