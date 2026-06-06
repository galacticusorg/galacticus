#!/usr/bin/env python3
# Determine the relation between stellar mass and SDSS r-band absolute magnitude using semi-analytic models from the Millennium
# Database. Specifically, the DeLucia2006a models are used - these correspond to De Lucia & Blaizot (2007;
# http://adsabs.harvard.edu/abs/2007MNRAS.375....2D) which use a Chabrier IMF, consistent with the IMF used by Li & White (2009).
# Andrew Benson (10-July-2012)

import os
import sys
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.cosmology import LambdaCDM
import astropy.units as u

execPath = os.environ.get('GALACTICUS_EXEC_PATH', '').rstrip('/') + '/'

# Define a working directory.
workDirectory = execPath + 'constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07/massLuminosityWork/'
os.makedirs(workDirectory, exist_ok=True)

# Define little Hubble parameter for the Millennium Simulation.
hubble = 0.7

# Define the limiting apparent magnitude of the Li & White (2009) sample.
apparentMagnitudeLimit = 17.6

# Define Millennium Simulation cosmology.
cosmology = LambdaCDM(H0=70.0, Om0=0.3, Ode0=0.7)

# Define the database URL.
databaseURL = 'http://gavo.mpa-garching.mpg.de/Millennium/?action=doQuery&SQL='

# Get list of snapshot numbers and corresponding redshifts.
snapshotsFile = workDirectory + 'millenniumDB_snapshots.csv'
if not os.path.exists(snapshotsFile):
    sqlQuery = 'select snapnum, z from millimil..Snapshots'
    subprocess.run(['wget', databaseURL + sqlQuery, '-O', snapshotsFile], check=True)

snapshotNumbers = []
redshifts       = []
with open(snapshotsFile, 'r') as f:
    for line in f:
        if line.startswith('#') or 'snapnum' in line:
            continue
        parts = line.strip().split(',')
        if len(parts) >= 2:
            try:
                snapshotNumbers.append(int  (parts[0]))
                redshifts      .append(float(parts[1]))
            except ValueError:
                continue
snapshotNumbers = np.array(snapshotNumbers)
redshifts       = np.array(redshifts      )

# Initialize vectors to store results.
redshiftTable     = [0.0]
limitingMassTable = [6.0]

# Loop over snapshots in reverse order.
for iSnapshot in range(len(snapshotNumbers) - 1, -1, -1):
    redshift = redshifts[iSnapshot]

    # Only consider reasonably low redshifts.
    if 0.0 < redshift <= 0.5:

        # Compute the limiting absolute magnitude at this redshift.
        dL = cosmology.luminosity_distance(redshift).to(u.pc).value
        limitingMagnitude = apparentMagnitudeLimit - 5.0 * np.log10(dL / 10.0)

        # Define the SQL query.
        snapnum  = snapshotNumbers[iSnapshot]
        sqlQuery = (
            'select mass.stellarMass, light.r_sdss'
            ' from millimil..DeLucia2006a as mass, millimil..DeLucia2006a_sdss2mass as light'
            f' where mass.snapnum = {snapnum}'
            ' and mass.galaxyID = light.galaxyID'
        )

        # Download the data.
        csvFile = workDirectory + f'millenniumDB_{snapnum}.csv'
        if not os.path.exists(csvFile):
            subprocess.run(['wget', databaseURL + sqlQuery, '-O', csvFile], check=True)

        # Read the data.
        mass      = []
        magnitude = []
        with open(csvFile, 'r') as f:
            for line in f:
                if line.startswith('#') or 'stellarMass' in line:
                    continue
                parts = line.strip().split(',')
                if len(parts) >= 2:
                    try:
                        mass     .append(float(parts[0]))
                        magnitude.append(float(parts[1]))
                    except ValueError:
                        continue
        if not mass:
            continue
        mass      = np.array(mass     )
        magnitude = np.array(magnitude)

        # Convert mass to logarithmic scale and Solar units.
        logMass = np.log10(1.0e10 * mass / hubble)

        # Define mass bins.
        logMassMinimum = 8.0
        logMassMaximum = 13.0
        logMassCount   = 10
        logMassBins    = np.linspace(logMassMinimum, logMassMaximum, logMassCount)

        # Find the median magnitude as a function of mass.
        medianMag = np.zeros(logMassCount)
        for iBin in range(logMassCount):
            if iBin < logMassCount - 1:
                inBin = np.where((logMass >= logMassBins[iBin]) & (logMass < logMassBins[iBin + 1]))[0]
            else:
                inBin = np.where(logMass >= logMassBins[iBin])[0]
            if len(inBin) > 0:
                medianMag[iBin] = np.median(magnitude[inBin])

        # Fit a degree-3 polynomial to magnitude vs. mass (only for bins with negative magnitude).
        nonZero = np.where(medianMag < 0.0)[0]
        if len(nonZero) < 4:
            continue
        coeffsDescending = np.polyfit(medianMag[nonZero], logMassBins[nonZero], 3)
        coeffsAscending  = coeffsDescending[::-1]

        # Compute the limiting mass at this redshift.
        limitingMass = (
            coeffsAscending[0] +
            coeffsAscending[1] * limitingMagnitude +
            coeffsAscending[2] * limitingMagnitude**2
        )

        redshiftTable    .append(redshift    )
        limitingMassTable.append(limitingMass)

redshiftTable     = np.array(redshiftTable    )
limitingMassTable = np.array(limitingMassTable)

# Generate a fit to the data (degree 5).
coeffsDescending = np.polyfit(limitingMassTable, redshiftTable, 5)
coeffsAscending  = coeffsDescending[::-1]
fitMass = np.linspace(8.0, 13.0, 1000)
fitRedshift = np.zeros(len(fitMass))
closing = ''
for i, c in enumerate(coeffsAscending):
    fitRedshift += c * fitMass ** i
    print(f'{c:.10e}d0', end='')
    if i == len(coeffsAscending) - 1:
        print(closing)
    else:
        print('+logarithmicMass*(', end='')
        closing += ')'

# Create a plot.
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_title(r'Redshift vs. Limiting Stellar Mass for Li \& White (2009) Sample')
ax.set_xlabel(r'Limiting stellar mass; $M_\star\;[M_\odot]$')
ax.set_ylabel(r'Redshift; $z$')
ax.set_xscale('log')
ax.set_xlim(1.0e8, 1.0e13)
ax.set_ylim(0.0, 0.5)
ax.plot(10.0 ** fitMass, fitRedshift,
        color='mediumseagreen', linewidth=2.5, label='Fit')
ax.plot(10.0 ** limitingMassTable, redshiftTable,
        'o', color='goldenrod', markersize=8, label='De Lucia (2006) SAM')
ax.legend(loc='lower left')
plt.tight_layout()
plotFile = workDirectory + 'massLuminosityRelation.pdf'
plt.savefig(plotFile)
plt.close()
