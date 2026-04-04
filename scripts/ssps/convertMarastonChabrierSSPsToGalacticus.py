#!/usr/bin/env python3
# Convert Maraston (Chabrier IMF) stellar population data to Galacticus' format.
# Andrew Benson (18-October-2010); ported to Python by Copilot (2026).

import os
import re
import numpy as np
import h5py

# Create a data directory.
dataDirectory = "SSP_Maraston"
os.makedirs(dataDirectory, exist_ok=True)

# Specify list of metallicities (exclude those with only crude time grids).
metallicityList = ["m0.0010", "m0.0100", "m0.0200", "m0.0400"]
metallicityValues = {
    "m0.0400": +0.35,
    "m0.0200": +0.00,
    "m0.0100": -0.33,
    "m0.0010": -1.35,
}

# Specify number of ages in the files.
ageCount = 220

# Loop over all metallicities.
metallicityData = []
fluxData  = None
ages      = None
lambdas   = None

for iMetal, metallicity in enumerate(metallicityList):
    metallicityData.append(metallicityValues[metallicity])

    fileName = os.path.join(dataDirectory, f"M05_Chabrier_FullSED_{metallicity}.dat")
    print(fileName)

    lastAge  = -1.0
    iAge     = -1
    ages_tmp = []
    wl_tmp   = []
    flux_tmp = []

    with open(fileName, 'r') as fh:
        for line in fh:
            if re.match(r'^\s*$', line):
                continue
            line = line.strip()
            cols = line.split()
            age  = float(cols[0])
            lam  = float(cols[2])
            flux = float(cols[3])

            if age != lastAge:
                if len(wl_tmp) > 0:
                    iAge += 1
                    nWl = len(flux_tmp)
                    if fluxData is None:
                        fluxData = np.zeros((nWl, ageCount, len(metallicityList)))
                    fluxData[:, iAge, iMetal] = flux_tmp
                wl_tmp   = []
                flux_tmp = []
                ages_tmp.append(age)
                lastAge = age

            wl_tmp.append(lam)
            flux_tmp.append(flux * lam**2)

    # Handle last age block.
    if len(wl_tmp) > 0:
        iAge += 1
        if fluxData is None:
            fluxData = np.zeros((len(flux_tmp), ageCount, len(metallicityList)))
        fluxData[:, iAge, iMetal] = flux_tmp

    if ages is None:
        ages    = np.array(ages_tmp)
        lambdas = np.array(wl_tmp)

# Convert fluxes to Lsolar/Hz.
solarLuminosity = 3.826e33
angstroms       = 1.0e-10
speedLight      = 2.998e8
fluxData *= angstroms / speedLight / solarLuminosity

# Convert ages to Gyr.
ages /= 1.0e9

# Create the HDF5 output file.
outFile = "data/SSP_Spectra_Maraston_imfChabrier.hdf5"
os.makedirs(os.path.dirname(outFile), exist_ok=True)
with h5py.File(outFile, 'w') as hdf:
    ds = hdf.create_dataset("ages", data=ages)
    ds.attrs["units"]     = "Gigayears"
    ds.attrs["unitsInSI"] = 3.15576e+16

    ds = hdf.create_dataset("wavelengths", data=lambdas)
    ds.attrs["units"]     = "Angstroms"
    ds.attrs["unitsInSI"] = 1.0e-10

    ds = hdf.create_dataset("metallicities", data=np.array(metallicityData))
    ds.attrs["units"] = "dex"
    ds.attrs["type"]  = "logarithmic, relative to Solar"

    ds = hdf.create_dataset("spectra", data=fluxData)
    ds.attrs["units"]     = "Lsolar/Hz"
    ds.attrs["unitsInSI"] = 3.827e+26

    hdf.attrs["fileFormat"] = np.int64(1)

    # Note: baseURL is not defined in this script; using the known URL.
    baseURL = "http://www-astro.physics.ox.ac.uk/~maraston/SSPn/SED/"
    src = hdf.require_group("source")
    src.attrs["source"]    = "Maraston, C. 2005, MNRAS, 362, 799; Maraston, C. 1998, MNRAS, 300, 872"
    src.attrs["sourceURL"] = ("http://adsabs.harvard.edu/abs/2005MNRAS.362..799M "
                              "http://adsabs.harvard.edu/abs/1998MNRAS.300..872M")
    src.attrs["URL"]       = baseURL + "Claudia%27s_Stellar_Population_Models.html"
    src.attrs["provenance"] = "Provided by Bruno Henriques <bhenriques@mpa-garching.mpg.de>; April 8, 2015"
