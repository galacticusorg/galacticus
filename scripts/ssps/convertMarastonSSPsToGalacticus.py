#!/usr/bin/env python3
# Convert Maraston stellar population data to Galacticus' format.
# Andrew Benson (18-October-2010); ported to Python by Copilot (2026).

import os
import re
import subprocess
import sys
import numpy as np
import h5py

# Define the base URL for the data.
baseURL = "http://www-astro.physics.ox.ac.uk/~maraston/SSPn/SED/"

# Specify the files to be downloaded.
downloadFiles = [
    "AgegridSSP_Mar05",
    "Sed_Mar05_SSP_Kroupa.tar.gz",
    "Sed_Mar05_SSP_Salpeter.tar.gz",
]

# Create a data directory.
dataDirectory = "SSP_Maraston"
os.makedirs(dataDirectory, exist_ok=True)

# Download data.
for fname in downloadFiles:
    subprocess.run(['wget', baseURL + fname, '-O', os.path.join(dataDirectory, fname)], check=True)
    if fname.endswith('.tar.gz'):
        subprocess.run(f"cd {dataDirectory}; tar xvfz {fname}", shell=True, check=True)

# Specify list of IMFs to convert.
IMFs = {
    "Kroupa":   {"label": "kr"},
    "Salpeter": {"label": "ss"},
}

# Specify horizontal branch morphologies to convert.
hbMorphologies = {
    # "Blue": {"label": "bhb"},  # Ignore blue horizontal branch files (only two ages).
    "Red": {"label": "rhb"},
}

# Specify list of metallicities (exclude those with only crude time grids).
metallicityList = ["z0001", "z001", "z002", "z004"]
metallicityValues = {
    "z007":  +0.67,
    "z004":  +0.35,
    "z002":  +0.00,
    "z001":  -0.33,
    "z0001": -1.35,
    "z10m4": -2.25,
}

# Count lines in the age grid file.
ageCount = 0
with open(os.path.join(dataDirectory, "AgegridSSP_Mar05"), 'r') as fh:
    for _ in fh:
        ageCount += 1

galDataPath = os.environ.get('GALACTICUS_DATA_PATH', '.')

# Loop over all IMFs.
for IMF, imfData in IMFs.items():

    # Loop over all horizontal branch morphologies.
    for hbMorphology, hbData in hbMorphologies.items():

        fluxData = None

        # Loop over all metallicities.
        metallicityData = []
        ages    = None
        lambdas = None

        for iMetal, metallicity in enumerate(metallicityList):
            metallicityData.append(metallicityValues[metallicity])

            # Construct file name.
            fileName = os.path.join(
                dataDirectory,
                f"sed.{imfData['label']}{metallicity}.{hbData['label']}"
            )

            # Open the file and read data.
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
                    age    = float(cols[0])
                    lam    = float(cols[2])
                    flux   = float(cols[3])

                    if age != lastAge:
                        if len(wl_tmp) > 0:
                            iAge += 1
                            nWl = len(flux_tmp)
                            if fluxData is None:
                                fluxData = np.zeros((len(metallicityList), ageCount, nWl))
                            fluxData[iMetal, iAge, :] = flux_tmp
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
                    fluxData = np.zeros((len(metallicityList), ageCount, len(flux_tmp)))
                fluxData[iMetal, iAge, :] = flux_tmp

            if ages is None:
                ages    = np.array(ages_tmp)
                lambdas = np.array(wl_tmp)

        # Convert fluxes to Lsun/Hz.
        solarLuminosity   = 3.826e33
        angstromsToMeters = 1.0e-10
        speedOfLight      = 2.998e8
        fluxData *= angstromsToMeters / solarLuminosity / speedOfLight

        # Create the HDF5 output file.
        outPath = os.path.join(
            galDataPath, 'static', 'stellarPopulations',
            f"SSP_Spectra_Maraston_hbMorphology{hbMorphology}_imf{IMF}.hdf5"
        )
        os.makedirs(os.path.dirname(outPath), exist_ok=True)
        with h5py.File(outPath, 'w') as hdf:
            hdf.create_dataset("ages",          data=ages)
            hdf.create_dataset("wavelengths",   data=lambdas)
            hdf.create_dataset("metallicities", data=np.array(metallicityData))
            hdf.create_dataset("spectra",       data=fluxData)

            src = hdf.require_group("source")
            src.attrs["source"]    = "Maraston, C. 2005, MNRAS, 362, 799; Maraston, C. 1998, MNRAS, 300, 872"
            src.attrs["sourceURL"] = ("http://adsabs.harvard.edu/abs/2005MNRAS.362..799M "
                                      "http://adsabs.harvard.edu/abs/1998MNRAS.300..872M")
            src.attrs["URL"] = baseURL + "Claudia%27s_Stellar_Population_Models.html"

# Remove the temporary data directory.
subprocess.run(['rm', '-rf', dataDirectory])
