#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Run models that check that the results of Johnson et al. (2021) can be reproduced.
# Andrew Benson (ported to Python)

np.random.seed(42)

# Make output directory.
subprocess.run("mkdir -p outputs/", shell=True)

# Specify the tests to run. Mean and tolerance targets are taken from Table 1 of Johnson et al. (2021).
tests = [
    {"suffix": "",              "mean": 1.084, "meanTolerance": 0.070, "scatter": 0.134, "scatterTolerance": 0.040},
    {"suffix": "Subsample1e8",  "mean": 1.084, "meanTolerance": 0.070, "scatter": 0.134, "scatterTolerance": 0.040},
    {"suffix": "Subsample1e7",  "mean": 1.084, "meanTolerance": 0.070, "scatter": 0.134, "scatterTolerance": 0.040},
]

for test in tests:
    suffix        = test["suffix"]
    suffixCapital = suffix[0].upper() + suffix[1:] if suffix else ""

    # Run the model.
    status = subprocess.run(
        f"cd ..; ./Galacticus.exe testSuite/parameters/concentrationDistributionJohnson2021{suffixCapital}.xml",
        shell=True
    )
    if status.returncode != 0:
        print(f"FAILED: Johnson2021 {suffix} concentration model failed to run")
        sys.exit(0)

    with h5py.File(f"outputs/concentrationDistributionJohnson2021{suffixCapital}.hdf5", "r") as model:
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

    mean    = np.mean(concentrationsMeasured)
    scatter = np.std(concentrationsMeasured)
    print(f"Mean   : {mean:5.3f} (target = {test['mean']} \u00b1 {test['meanTolerance']})")
    print(f"Scatter: {scatter:5.3f} (target = {test['scatter']} \u00b1 {test['scatterTolerance']})")
    status_str = "SUCCESS" if (
        abs(mean    - test["mean"])    < test["meanTolerance"] and
        abs(scatter - test["scatter"]) < test["scatterTolerance"]
    ) else "FAILED"
    print(f"{status_str}: Johnson2021 {suffix} concentration")
