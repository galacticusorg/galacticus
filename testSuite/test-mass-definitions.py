#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Test mass definition conversions.
# Andrew Benson (ported to Python)

PI = np.pi

# Make output directory.
subprocess.run("mkdir -p outputs/", shell=True)

# Iterate over time options.
for optionTime in ("current", "infall"):
    # Run the model.
    optionTimeCapitalized = optionTime[0].upper() + optionTime[1:]
    status = subprocess.run(
        f"cd ..; ./Galacticus.exe testSuite/parameters/massDefinitionsTime{optionTimeCapitalized}.xml",
        shell=True
    )
    if status.returncode != 0:
        print("FAIL: failed to run mass definitions model")
        continue

    # Read all data.
    with h5py.File(f"outputs/massDefinitionsTime{optionTimeCapitalized}.hdf5", "r") as model:
        cosmology       = model["Parameters/cosmologyParameters"]
        OmegaMatter     = cosmology.attrs["OmegaMatter"]
        HubbleConstant  = cosmology.attrs["HubbleConstant"]
        nodes           = model["Outputs/Output1/nodeData"]
        halos           = {}
        for key in ("redshift", "redshiftLastIsolated", "basicMass", "massHaloEnclosedCurrent", "darkMatterOnlyRadiusVirial"):
            halos[key] = nodes[key][:]

    # Find mean density of the universe at the present time.
    gravitationalConstant = 4.301e-9
    densityMean           = OmegaMatter * 3.0 * HubbleConstant**2 / 8.0 / PI / gravitationalConstant
    if optionTime == "current":
        halos["densityMean"] = densityMean
    elif optionTime == "infall":
        halos["densityMean"] = densityMean * ((1.0 + halos["redshiftLastIsolated"]) / (1.0 + halos["redshift"]))**3
    else:
        raise ValueError(f"unrecognized time option: {optionTime}")

    # Compute the target density, which is 200 times the mean density.
    densityTarget = 200.0 * halos["densityMean"]
    # Find the density of halos under the internal mass definition.
    halos["density"]       = halos["basicMass"] / (4.0 * PI * halos["darkMatterOnlyRadiusVirial"]**3 / 3.0)
    # For isothermal halos, the enclosed mass is inversely proportional to the square root of the enclosed density.
    halos["massConverted"] = halos["basicMass"] * np.sqrt(halos["density"] / densityTarget)
    # Compute the difference between what the mass should be, and what Galacticus computed.
    error = np.abs(halos["massHaloEnclosedCurrent"] - halos["massConverted"]) / halos["massConverted"]
    if np.all(error < 1.0e-4):
        print(f"SUCCESS: mass definitions for '{optionTime}' time")
    else:
        print(f"FAIL: mass definitions for '{optionTime}' time")
