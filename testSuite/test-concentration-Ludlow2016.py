#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Run models that check that the results of Benson, Ludlow, & Cole (2019) can be reproduced.
# Andrew Benson (ported to Python)

np.random.seed(42)

# Make output directory.
subprocess.run("mkdir -p outputs/", shell=True)

# Specify the tests to run. Mean and tolerance targets are taken from Table 2 of Benson, Ludlow, & Cole (2019).
tests = [
    {"suffix": "",              "mean": 1.104, "meanTolerance": 0.040, "scatter": 0.103, "scatterTolerance": 0.008},
    {"suffix": "environmental", "mean": 1.077, "meanTolerance": 0.040, "scatter": 0.114, "scatterTolerance": 0.016},
]

for test in tests:
    suffix        = test["suffix"]
    suffixCapital = suffix[0].upper() + suffix[1:] if suffix else ""

    # Run the model.
    status = subprocess.run(
        f"cd ..; ./Galacticus.exe testSuite/parameters/concentrationDistributionLudlow2016{suffixCapital}.xml",
        shell=True
    )
    if status.returncode != 0:
        print(f"FAILED: Ludlow2016 {suffix} concentration model failed to run")
        sys.exit(0)

    # Select halos in the mass range used for the results in Table 2 of Benson, Ludlow, & Cole (2019).
    with h5py.File(f"outputs/concentrationDistributionLudlow2016{suffixCapital}.hdf5", "r") as model:
        nodeData      = model["Outputs/Output1/nodeData"]
        nodeIsIsolated= nodeData["nodeIsIsolated"][:]
        massHalo      = nodeData["massHaloEnclosedCurrent"][:]
        concentration = nodeData["concentration"][:]

    # Compute N-body measurement uncertainties.
    massParticle             = 1.6e5  # COCO simulation particle mass.
    countParticles           = massHalo / massParticle
    alpha                    = -0.20 + 1.46 * np.log10(concentration) - 0.25 * np.log10(concentration)**2
    b                        = -0.54
    concentrationUncertainty = 10.0**(alpha + b * np.log10(countParticles))
    massUncertainty          = 0.135 * (1000.0 / countParticles)**(1.0/3.0)

    # Construct perturbed halo masses.
    massMeasured = massHalo * np.exp(massUncertainty * np.random.randn(len(massHalo)))

    # Select halos for inclusion.
    selection = np.where(
        (nodeIsIsolated == 1)
        & (np.log10(massMeasured) >= 9.402)
        & (np.log10(massMeasured)  < 9.902)
    )[0]

    # Construct perturbed concentrations.
    concentrationsMeasured = np.log10(concentration[selection]) + concentrationUncertainty[selection] * np.random.randn(len(selection))

    # Evaluate the mean and scatter.
    mean    = np.mean(concentrationsMeasured)
    scatter = np.std(concentrationsMeasured)
    print(f"Mean - {mean} (target = {test['mean']} \u00b1 {test['meanTolerance']})")
    print(f"Scatter - {scatter} (target = {test['scatter']} \u00b1 {test['scatterTolerance']})")
    status_str = "SUCCESS" if (
        abs(mean    - test["mean"])    < test["meanTolerance"] and
        abs(scatter - test["scatter"]) < test["scatterTolerance"]
    ) else "FAILED"
    print(f"{status_str}: Ludlow2016 {suffix} concentration")
