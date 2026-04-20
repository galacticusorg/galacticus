#!/usr/bin/env python3
# Convert Bruzual & Charlot stellar population data to Galacticus' format.
# Andrew Benson (19-October-2010); ported to Python by Copilot (2026).

import argparse
import os
import re
import subprocess
import sys
import numpy as np
import h5py

# Parse command line arguments.
parser = argparse.ArgumentParser(
    prog='convertBC2003SSPsToGalacticus.py',
    description='Convert Bruzual & Charlot (2003) stellar population data to Galacticus format.'
)
parser.add_argument('dataDirectory', help='Directory containing the BC2003 data archives')
args = parser.parse_args()
dataDirectory = args.dataDirectory

# Specify models to convert.
modelsToConvert = ["padova_1994_salpeter_imf"]

# Ensure data files have been downloaded.
for model in modelsToConvert:
    archivePath = os.path.join(dataDirectory, f"bc03.models.{model}.tar.gz")
    if not os.path.exists(archivePath):
        sys.exit("Please download the model archives from the BC2003 website before running this script.")

# Lookup table for metallicity labels.
metallicityLookup = {
    '22': 0.0001,
    '32': 0.0004,
    '42': 0.0040,
    '52': 0.0080,
    '62': 0.0200,
    '72': 0.0500,
}

# Lookup table for resolutions.
resolutionLookup = {
    'lr': 'lowResolution',
    'hr': 'highResolution',
}

# Loop over models.
for model in modelsToConvert:
    # Extract the IMF name and other data from the model name.
    m = re.match(r'([a-z]+)_(\d+)_([a-z]+)_imf', model)
    if not m:
        sys.exit("Model name is in unrecognized format.")
    tracks = m.group(1).capitalize()
    year   = m.group(2)
    IMF    = m.group(3)

    # Unpack the data.
    subprocess.run(f"cd {dataDirectory}; tar xvfz bc03.models.{model}.tar.gz", shell=True, check=True)

    # Build the base directory for this model.
    baseDirectory = os.path.join(dataDirectory, "bc03", "models", f"{tracks}{year}", IMF)

    # Unpack the ASCII data files.
    subprocess.run(
        f'find {baseDirectory} -name "*_ASCII.gz" | xargs --no-run-if-empty gunzip -f',
        shell=True
    )

    # Loop over low and high-res spectra.
    for resolution in ('lr', 'hr'):

        # Build a list of files to process.
        files = {}
        for fname in os.listdir(baseDirectory):
            m2 = re.match(rf'bc2003_{resolution}_m(\d+)_[a-z]+_ssp\.ised_ASCII$', fname)
            if m2:
                metallicity = m2.group(1)
                files[metallicity] = fname

        # Skip if no files found for this resolution.
        if not files:
            continue

        # Reset datasets.
        metallicities     = []
        ageDataset        = None
        wavelengthDataset = None
        spectra           = None

        # Process files.
        sortedKeys = sorted(files.keys())
        nMetallicities = len(sortedKeys)

        for iMetallicity, fileKey in enumerate(sortedKeys):
            ages        = []
            wavelengths = []

            fileName = os.path.join(baseDirectory, files[fileKey])
            print(f"Processing file: {fileName}")

            # Extract the metallicity.
            metallicities.append(metallicityLookup[fileKey])

            def read_record(fh, target_count):
                """Read tokens from successive lines until target_count values collected."""
                values = []
                while len(values) < target_count:
                    line = fh.readline()
                    if not line:
                        break
                    line = line.strip()
                    cols = line.split()
                    values.extend(float(c) for c in cols)
                return values

            with open(fileName, 'r') as fh:
                # Read ages: first token on first line is count, then values.
                line = fh.readline().strip()
                cols = line.split()
                ageCount = int(cols[0])
                ages = [float(c) for c in cols[1:]]
                while len(ages) < ageCount:
                    extra = read_record(fh, ageCount - len(ages))
                    ages.extend(extra)
                ages = np.array(ages[:ageCount])

            # Store ages if not already done. (Convert to Gyr as we do so.)
            if ageDataset is None:
                ageDataset = ages / 1.0e9

            with open(fileName, 'r') as fh:
                # Skip the age lines (re-read past them).
                remaining_ages = ageCount
                first_age_line = True
                while remaining_ages > 0:
                    line = fh.readline().strip()
                    cols = line.split()
                    if first_age_line:
                        remaining_ages -= len(cols) - 1  # subtract count field
                        first_age_line = False
                    else:
                        remaining_ages -= len(cols)

                # Skip over lines prior to those beginning with "Padova".
                gotPadova = False
                lineBuffer = None
                while True:
                    line = fh.readline()
                    if not line:
                        break
                    if line.startswith('Padova'):
                        gotPadova = True
                    else:
                        if gotPadova:
                            # This is the first non-Padova line after Padova lines.
                            lineBuffer = line
                            break
                        # else: still scanning pre-Padova lines, discard

                # Read wavelengths: first token is count, then values.
                wavelengths = []
                if lineBuffer:
                    cols = lineBuffer.strip().split()
                    wavelengthCount = int(cols[0])
                    wavelengths.extend(float(c) for c in cols[1:])
                else:
                    line = fh.readline().strip()
                    cols = line.split()
                    wavelengthCount = int(cols[0])
                    wavelengths.extend(float(c) for c in cols[1:])
                while len(wavelengths) < wavelengthCount:
                    line = fh.readline().strip()
                    wavelengths.extend(float(c) for c in line.split())
                wavelengths = np.array(wavelengths[:wavelengthCount])
                wavelengthDataset = wavelengths

                # Create array to hold spectra.
                if spectra is None:
                    spectra = np.zeros((nMetallicities, ageCount, wavelengthCount))

                # Loop over ages.
                for iAge in range(ageCount):
                    # Read spectrum: first token is count, then values.
                    spectrum = []
                    line = fh.readline().strip()
                    cols = line.split()
                    specCount = int(cols[0])
                    spectrum.extend(float(c) for c in cols[1:])
                    while len(spectrum) < specCount:
                        line = fh.readline().strip()
                        spectrum.extend(float(c) for c in line.split())
                    spectrum = np.array(spectrum[:specCount])

                    # Store the spectrum (multiply by wavelength^2 for converting to LSolar/Hz).
                    spectra[iMetallicity, iAge, :] = spectrum * wavelengthDataset**2

                    # Skip the fitting function records.
                    fit = []
                    line = fh.readline().strip()
                    cols = line.split()
                    fitCount = int(cols[0])
                    fit.extend(float(c) for c in cols[1:])
                    while len(fit) < fitCount:
                        line = fh.readline().strip()
                        fit.extend(float(c) for c in line.split())

        # Convert metallicities to logarithmic relative to Solar.
        metallicitySolar = 0.0188
        metallicities = np.log10(np.array(metallicities) / metallicitySolar)

        # Convert spectra to correct units.
        angstromsToMeters = 1.0e-10
        speedOfLight      = 2.998e8
        spectra *= angstromsToMeters / speedOfLight

        # Create the HDF5 output file.
        IMF_cap = IMF.capitalize()
        outFile = f"SSP_Spectra_BC2003_{resolutionLookup[resolution]}_imf{IMF_cap}.hdf5"
        with h5py.File(outFile, 'w') as hdf:
            hdf.create_dataset("ages",          data=ageDataset)
            hdf.create_dataset("wavelengths",   data=wavelengthDataset)
            hdf.create_dataset("metallicities", data=metallicities)
            hdf.create_dataset("spectra",       data=spectra)

            src = hdf.require_group("source")
            src.attrs["source"]    = "Bruzual and Charlot, 2003, MNRAS, 344, 1000"
            src.attrs["sourceURL"] = "http://adsabs.harvard.edu/abs/2003MNRAS.344.1000B"
            src.attrs["URL"]       = "http://www2.iap.fr/users/charlot/bc2003/"

print(f"Conversion finished. You may delete the {dataDirectory} directory now.")
