#!/usr/bin/env python3
# Convert BPASS v2.3 stellar population data to Galacticus' format.
# Andrew Benson (24-March-2022); ported to Python by Copilot (2026).

import argparse
import os
import re
import sys
import datetime
import math
import numpy as np
import h5py

# Parse command line arguments.
parser = argparse.ArgumentParser(
    prog='convertBPASSv2.3SSPsToGalacticus.py',
    description='Convert BPASS v2.3 stellar population data to Galacticus format.'
)
parser.add_argument('dataDirectory', help='Directory containing the BPASSv2.3 spectra files')
args = parser.parse_args()
dataDirectory = args.dataDirectory

# Construct lists of ages and metallicities.
metallicityValues = np.array([1.0e-5, 1.0e-4, 1.0e-3, 2.0e-3, 3.0e-3, 4.0e-3, 6.0e-3,
                               8.0e-3, 1.0e-2, 1.4e-2, 2.0e-2, 3.0e-2, 4.0e-2])
metallicitiesSolar = np.log10(metallicityValues / 0.02)
ages               = 10.0 ** (np.arange(51) * 0.1 - 3.0)
countMetallicities = len(metallicityValues)
countAges          = len(ages)
countWavelengths   = 100000

galDataPath = os.environ.get('GALACTICUS_DATA_PATH', '.')


def metallicityZLabel(z):
    """Return the z-label string used in BPASS v2.3 filenames for metallicity z."""
    if z < 0.99e-3:
        return "em" + str(int(-round(math.log10(z))))
    else:
        return f"{int(round(1000.0 * z)):03d}"


# Iterate over all spectra files in the data directory.
filePattern = re.compile(r'^spectra-(bin|sin)-imf(135_300)\.a([\+\-\d]+)\.z([\dem]+)\.dat$')

for modelName in sorted(os.listdir(dataDirectory)):
    m = filePattern.match(modelName)
    if not m:
        continue

    binarity   = m.group(1)
    imfLabel   = m.group(2)
    alphaLabel = m.group(3)

    outPath = os.path.join(
        galDataPath, 'dynamic', 'stellarPopulations',
        f"SSP_Spectra_BPASSv2.3_{binarity}-imf{imfLabel}-alpha{alphaLabel}.hdf5"
    )
    if os.path.exists(outPath):
        continue

    spectra     = np.zeros((countWavelengths, countAges, countMetallicities))
    wavelengths = np.zeros(countWavelengths)

    for iMetallicity, z in enumerate(metallicityValues):
        zLabel = metallicityZLabel(z)
        modelFileName = os.path.join(
            dataDirectory,
            f"spectra-{binarity}-imf{imfLabel}.a{alphaLabel}.z{zLabel}.dat"
        )
        print(modelFileName)

        with open(modelFileName, 'r') as fh:
            for iWavelength, line in enumerate(fh):
                cols = line.split()
                wavelengths[iWavelength] = float(cols[0])
                spectra[iWavelength, :, iMetallicity] = [float(cols[k]) for k in range(1, countAges + 1)]

        # Convert units (spectra are for 10^6 Msun burst).
        angstromsToMeters = 1.0e-10
        speedOfLight      = 2.998e8
        for i in range(countAges):
            spectra[:, i, iMetallicity] *= wavelengths**2 * angstromsToMeters / speedOfLight / 1.0e6

    # Write data to Galacticus' format file.
    os.makedirs(os.path.dirname(outPath), exist_ok=True)
    with h5py.File(outPath, 'w') as hdf:
        ds = hdf.create_dataset("ages", data=ages)
        ds.attrs["units"]     = "Gigayears"
        ds.attrs["unitsInSI"] = 3.15576e+16

        ds = hdf.create_dataset("wavelengths", data=wavelengths)
        ds.attrs["units"]     = "Angstroms"
        ds.attrs["unitsInSI"] = 1.0e-10

        ds = hdf.create_dataset("metallicities", data=metallicitiesSolar)
        ds.attrs["units"] = "dex"
        ds.attrs["type"]  = "logarithmic, relative to Solar"

        ds = hdf.create_dataset("spectra", data=spectra)
        ds.attrs["units"]     = "L\u2609/Hz"
        ds.attrs["unitsInSI"] = 3.827e+26

        # Add metadata.
        references = [
            "Eldridge, Stanway et al. (2017; PASA, 34, 58)",
            "Stanway & Eldridge et al. (2018; MNRAS, 479, 75)",
            "Byrne et al. (2022, MNRAS in press)",
        ]
        referenceURLs = [
            "https://ui.adsabs.harvard.edu/abs/2017PASA...34...58E",
            "https://ui.adsabs.harvard.edu/abs/2018MNRAS.479...75S",
            "https://warwick.ac.uk/fac/sci/physics/research/astro/research/catalogues/bpass/bpassv2p3/bpassv2p3_release_paper.pdf",
        ]

        alphaValue = int(alphaLabel) / 10.0
        description = "Simple stellar population spectra from the BPASS v2.3 library"
        description += " including binary stars" if binarity == "bin" else " not including binary stars"
        description += f" with \u03b1-alement enhancement of \u0394log\u2081\u2080(\u03b1/Fe)={alphaValue:+4.1f}"
        alpha = f"\u03b1-alement enhancement of \u0394log\u2081\u2080(\u03b1/Fe)={alphaValue:+4.1f}."

        m_chab = re.match(r'_?chab(\d+)', imfLabel)
        m_all  = re.match(r'(\d+)all_(\d+)', imfLabel)
        m_pl   = re.match(r'(\d+)_(\d+)', imfLabel)
        if m_chab:
            description += f" for a Chabrier (2003) IMF with upper mass of {m_chab.group(1)}M\u2609."
            imf = f"Chabrier (2003) IMF with upper mass of {m_chab.group(1)}M\u2609."
        elif m_all:
            slope = -int(m_all.group(1)) / 100.0 - 1.0
            description += f" for a power law IMF with slope {slope:5.2f} from 0.1M\u2609 to {m_all.group(2)}M\u2609."
            imf = f"Power law IMF with slope {slope:5.2f} from 0.1M\u2609 to {m_all.group(2)}M\u2609."
        elif m_pl:
            slope1 = -1.30
            slope2 = -int(m_pl.group(1)) / 100.0 - 1.0
            description += (f" for a broken power law IMF with slope {slope1:5.2f} from 0.1M\u2609 to 0.5M\u2609"
                            f" and a slope of {slope2:5.2f} from 0.5M\u2609 to {m_pl.group(2)}M\u2609.")
            imf = (f"Broken power law IMF with slope {slope1:5.2f} from 0.1M\u2609 to 0.5M\u2609"
                   f" and a slope of {slope2:5.2f} from 0.5M\u2609 to {m_pl.group(2)}M\u2609.")
        else:
            sys.exit('failed to parse IMF label')

        now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='milliseconds')
        hdf.attrs["createdBy"]           = "Andrew Benson <abenson@obs.carnegiescience.edu>"
        hdf.attrs["description"]         = description
        hdf.attrs["url"]                 = "http://bpass.auckland.ac.nz"
        hdf.attrs["fileFormat"]          = np.int64(1)
        hdf.attrs["timeStamp"]           = now
        hdf.attrs["references"]          = "; ".join(references)
        hdf.attrs["referenceURLs"]       = "; ".join(referenceURLs)
        hdf.attrs["initialMassFunction"] = imf
        hdf.attrs["alphaEnhancement"]    = alpha
