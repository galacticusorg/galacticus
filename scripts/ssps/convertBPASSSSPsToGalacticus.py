#!/usr/bin/env python3
# Convert BPASS stellar population data to Galacticus' format.
# Andrew Benson (03-June-2014); ported to Python by Copilot (2026).

import argparse
import os
import re
import subprocess
import sys
import datetime
import numpy as np
import h5py

# Parse command line arguments.
parser = argparse.ArgumentParser(
    prog='convertBPASSSSPsToGalacticus.py',
    description='Convert BPASS stellar population data to Galacticus format.'
)
parser.add_argument('dataDirectory', help='Directory in which to place/find the BPASS data')
args = parser.parse_args()
dataDirectory = args.dataDirectory

# Download the archived data.
fileURL = "http://flexiblelearning.auckland.ac.nz/bpass/2/files/"
files = [
    {'metallicity': 0.001, 'fileName': 'sed_bpass_z001_tar.gz'},
    {'metallicity': 0.004, 'fileName': 'sed_bpass_z004_tar.gz'},
    {'metallicity': 0.008, 'fileName': 'sed_bpass_z008_tar.gz'},
    {'metallicity': 0.020, 'fileName': 'sed_bpass_z020_tar.gz'},
    {'metallicity': 0.040, 'fileName': 'sed_bpass_z040_tar.gz'},
]

for f in files:
    destPath = os.path.join(dataDirectory, f['fileName'])
    if not os.path.exists(destPath):
        subprocess.run(['wget', fileURL + f['fileName'], '-O', destPath], check=True)
    if not os.path.exists(destPath):
        sys.exit(f"Convert_BPASS_SSPs_to_Galacticus.py: failed to download file {f['fileName']}")

# Extract the archives.
for f in files:
    m = re.search(r'(z\d+)', f['fileName'])
    if m:
        f['zLabel'] = m.group(1)
    sampleFile = os.path.join(dataDirectory, 'SEDS', f"sed.bpass.constant.cloudy.bin.{f['zLabel']}")
    if not os.path.exists(sampleFile):
        subprocess.run(f"cd {dataDirectory}; tar xvfz {f['fileName']}", shell=True, check=True)

# Extract population ages.
ages = []
readmePath = os.path.join(dataDirectory, 'SEDS', 'sed.bpass.readme.txt')
with open(readmePath, 'r') as fh:
    for line in fh:
        m = re.match(r'^\s*\d+\)\s*Flux\( log\(Age\/yrs\)=([\d\.]+)\s*\) \/  L\(Sun\)\/A\.', line)
        if m:
            ages.append(float(m.group(1)))
ageCount = len(ages)
# Ages are given in log10 years; convert to Gyr.
ages = np.array([10.0**(a - 9.0) for a in ages])

# Extract metallicities and convert to log10 relative to Solar.
metallicityValues = np.array([f['metallicity'] for f in files])
metallicitySolar  = 0.0188
metallicities     = np.log10(metallicityValues / metallicitySolar)

galDataPath = os.environ.get('GALACTICUS_DATA_PATH', '.')

# Iterate over single and binary populations.
for population in ('single', 'binary'):

    wavelengths = None
    spectra     = None

    for iFile, f in enumerate(files):
        fileName = os.path.join(
            dataDirectory, 'SEDS',
            f"sed.bpass.instant.nocont.{population[:3]}.{f['zLabel']}"
        )

        # First pass: get wavelength grid.
        wl_list = []
        with open(fileName, 'r') as fh:
            for line in fh:
                line = line.strip()
                cols = line.split()
                wl_list.append(float(cols[0]))
        wavelengths = np.array(wl_list)
        nWavelengths = len(wavelengths)

        if spectra is None:
            spectra = np.zeros((nWavelengths, ageCount, len(files)))

        # Second pass: read spectra.
        with open(fileName, 'r') as fh:
            for j, line in enumerate(fh):
                line = line.strip()
                cols = line.split()
                spectra[j, :, iFile] = [float(cols[k]) for k in range(1, ageCount + 1)]

        # Convert units (spectra are for 10^6 Msun burst).
        angstromsToMeters = 1.0e-10
        speedOfLight      = 2.998e8
        for i in range(ageCount):
            spectra[:, i, iFile] *= wavelengths**2 * angstromsToMeters / speedOfLight / 1.0e6

    # Write data to Galacticus' format file.
    outPath = os.path.join(galDataPath, 'dynamic', 'stellarPopulations', f"SSP_Spectra_BPASS_{population}.hdf5")
    os.makedirs(os.path.dirname(outPath), exist_ok=True)
    with h5py.File(outPath, 'w') as hdf:
        ds = hdf.create_dataset("ages", data=ages)
        ds.attrs["units"]     = "Gigayears"
        ds.attrs["unitsInSI"] = 3.15576e+16

        ds = hdf.create_dataset("wavelengths", data=wavelengths)
        ds.attrs["units"]     = "Angstroms"
        ds.attrs["unitsInSI"] = 1.0e-10

        ds = hdf.create_dataset("metallicities", data=metallicities)
        ds.attrs["units"] = "dex"
        ds.attrs["type"]  = "logarithmic, relative to Solar"

        ds = hdf.create_dataset("spectra", data=spectra)
        ds.attrs["units"]     = "Lsolar/Hz"
        ds.attrs["unitsInSI"] = 3.827e+26

        # Add metadata.
        references = [
            "Eldridge & Stanway, 2012, MNRAS, 419, 479",
            "Eldridge, Langer & Tout, 2011, MNRAS, 414, 3501",
            "Eldridge & Stanway, 2009, MNRAS, 400, 1019",
            "Eldridge, Izzard & Tout, 2008, MNRAS, 384, 1109",
        ]
        referenceURLs = [
            "http://adsabs.harvard.edu/abs/2012MNRAS.419..479E",
            "http://adsabs.harvard.edu/abs/2011MNRAS.414.3501E",
            "http://adsabs.harvard.edu/abs/2009MNRAS.400.1019E",
            "http://adsabs.harvard.edu/abs/2008MNRAS.384.1109E",
        ]
        now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='milliseconds')
        hdf.attrs["createdBy"]           = "Andrew Benson <abenson@obs.carnegiescience.edu>"
        hdf.attrs["description"]         = "Simple stellar population spectra from the BPASS library for the canonical BPASS IMF."
        hdf.attrs["url"]                 = "http://www.bpass.org.uk/"
        hdf.attrs["fileFormat"]          = np.int64(1)
        hdf.attrs["timeStamp"]           = now
        hdf.attrs["references"]          = "; ".join(references)
        hdf.attrs["referenceURLs"]       = "; ".join(referenceURLs)
        hdf.attrs["initialMassFunction"] = "Slope of -1.3 between 0.1 and 0.5M\u2609 and a typical Salpeter slope of -2.35 above 0.5M\u2609 and a maximum mass of 120M\u2609."
