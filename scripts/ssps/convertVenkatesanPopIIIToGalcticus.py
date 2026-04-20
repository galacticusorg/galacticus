#!/usr/bin/env python3
# Convert Aparna Venkatesan's PopIII SEDs to Galacticus' format.
# Andrew Benson (25-June-2012); ported to Python by Copilot (2026).

import os
import datetime
import numpy as np
import h5py

# Specify ages and metallicities.
ages          = np.array([0.0, 1.0e-2])
metallicities = np.array([-4.0])

# Specify files to convert.
files = [
    {
        'file':        'data/pop3-inst-10-140.txt',
        'out':         'SSP_Spectra_Venkatesan_PopIII_m10:140.hdf5',
        'massMinimum': '10',
        'massMaximum': '140',
    },
    {
        'file':        'data/pop3-inst.txt',
        'out':         'SSP_Spectra_Venkatesan_PopIII_m1:100.hdf5',
        'massMinimum': '1',
        'massMaximum': '100',
    },
]

# Loop over files.
for fileSpec in files:

    # Read the text file.
    wavelength = []
    youngSED   = []
    oldSED     = []

    with open(fileSpec['file'], 'r') as fh:
        next(fh)  # skip header line
        for line in fh:
            line = line.strip()
            cols = line.split()
            wavelength.append(float(cols[0]))
            youngSED.append(float(cols[1]))
            oldSED.append(float(cols[2]))

    wavelength = np.array(wavelength)
    youngSED   = np.array(youngSED)
    oldSED     = np.array(oldSED)

    # Combine spectra into a single array: shape (nWavelength, nAge, nMetallicity).
    nWavelength = len(wavelength)
    spectra = np.zeros((nWavelength, 2, 1))
    spectra[:, 0, 0] = youngSED * wavelength**2
    spectra[:, 1, 0] = oldSED   * wavelength**2

    # Convert spectra to preferred units.
    solarLuminosity = 3.845e33  # ergs/s
    speedLight      = 2.998e18  # Angstroms/s
    spectra /= 1.0e6              # convert to 1 Msun population
    spectra /= solarLuminosity * speedLight

    # Write the data to file.
    outPath = fileSpec['out']
    if os.path.exists(outPath):
        os.remove(outPath)

    with h5py.File(outPath, 'w') as hdf:
        hdf.create_dataset("ages",          data=ages)
        hdf.create_dataset("metallicities", data=metallicities)
        hdf.create_dataset("wavelengths",   data=wavelength)
        hdf.create_dataset("spectra",       data=spectra)

        now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='milliseconds')
        hdf.attrs["fileFormat"]  = np.int64(1)
        hdf.attrs["createdBy"]   = "Galacticus"
        hdf.attrs["description"] = (
            f"Simple stellar population spectra for Population III stars from "
            f"Venkatesan et al. (2003; ApJ; 584; 621; "
            f"http://adsabs.harvard.edu/abs/2003ApJ...584..621V) using Salpeter IMF "
            f"from M={fileSpec['massMinimum']} to {fileSpec['massMaximum']} M\u2609"
        )
        hdf.attrs["timestamp"] = now

        ds = hdf["ages"]
        ds.attrs["unitsInSI"] = 3.15576e+16
        ds.attrs["units"]     = "Gigayears"

        ds = hdf["metallicities"]
        ds.attrs["units"] = "Solar"

        ds = hdf["wavelengths"]
        ds.attrs["unitsInSI"] = 1.0e-10
        ds.attrs["units"]     = "Angstroms"

        ds = hdf["spectra"]
        ds.attrs["unitsInSI"] = 3.827e+33
        ds.attrs["units"]     = "Lsolar/Hz"

        src = hdf.require_group("source")
        src.attrs["source"]    = "Venkatesan et al. (2003; ApJ; 584; 621)"
        src.attrs["sourceURL"] = "http://adsabs.harvard.edu/abs/2003ApJ...584..621V"

        imf = hdf.require_group("initialMassFunction")
        imf.attrs["shape"]       = "Salpeter"
        imf.attrs["massMinimum"] = fileSpec['massMinimum'] + " M\u2609"
        imf.attrs["massMaximum"] = fileSpec['massMaximum'] + " M\u2609"
