#!/usr/bin/env python3
# Determine the relation between stellar mass and Spitzer IRAC 4.5um apparent magnitude using semi-analytic models from the
# Millennium Database. Specifically, the Henriques2012a models (http://adsabs.harvard.edu/abs/2012MNRAS.421.2904H) are used -
# these correspond to Guo et al. (2011; http://adsabs.harvard.edu/abs/2011MNRAS.413..101G) which use a Chabrier IMF. We therefore
# include a correction to the Salpeter IMF assumed by Caputi et al. (2011).
# Andrew Benson (21-August-2012)

import os
import sys
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

execPath = os.environ.get('GALACTICUS_EXEC_PATH', '').rstrip('/') + '/'

myMillenniumUser   = os.environ.get('MYMILLENNIUM_USER'  , '')
myMillenniumPasswd = os.environ.get('MYMILLENNIUM_PASSWD', '')

# Define a working directory.
workDirectory = execPath + 'constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/massLuminosityWork/'
os.makedirs(workDirectory, exist_ok=True)

# Define little Hubble parameter for the Millennium Simulation.
hubble = 0.7

# Define the database URL.
databaseURL = 'http://gavo.mpa-garching.mpg.de/MyMillennium?action=doQuery&SQL='

# Define magnitude range for sample.
chabrierToSalpeter = 1.8
magnitudeLimit     = 24.0
magnitudeWidth     =  0.2
magnitudeMinimum   = magnitudeLimit - 0.5 * magnitudeWidth - 2.5 * np.log10(chabrierToSalpeter)
magnitudeMaximum   = magnitudeLimit + 0.5 * magnitudeWidth - 2.5 * np.log10(chabrierToSalpeter)

# Define the SQL query.
sqlQuery = (
    'select light.z_app, mass.stellarmass, light.i2'
    ' from Guo2010a..MR as mass, Henriques2012a.wmap1.BC03_001 as light'
    ' where mass.galaxyid = light.galid'
    f' and light.i2 < {magnitudeMaximum} and light.i2 > {magnitudeMinimum}'
    ' and light.z_app > 2.5 and light.z_app < 5.5'
)

# Download the data.
csvFile = workDirectory + 'millenniumDB.csv'
if not os.path.exists(csvFile):
    subprocess.run([
        'wget', '--http-user='+myMillenniumUser, '--http-passwd='+myMillenniumPasswd,
        databaseURL + sqlQuery, '-O', csvFile
    ], check=True)

# Read the data (skip header lines starting with # or containing 'stellarMass').
redshift  = []
mass      = []
magnitude = []
with open(csvFile, 'r') as f:
    for line in f:
        if line.startswith('#') or 'stellarMass' in line:
            continue
        parts = line.strip().split(',')
        if len(parts) >= 3:
            try:
                redshift .append(float(parts[0]))
                mass     .append(float(parts[1]))
                magnitude.append(float(parts[2]))
            except ValueError:
                continue
redshift  = np.array(redshift )
mass      = np.array(mass     )
magnitude = np.array(magnitude)

# Convert mass to logarithmic and scale to Solar units.
logMass = np.log10(1.0e10 * mass / hubble)

# Define redshift bins.
redshiftMinimum = 2.5
redshiftMaximum = 5.5
redshiftCount   = 10
redshiftBins    = np.linspace(redshiftMinimum, redshiftMaximum, redshiftCount)

# Find the median stellar mass as a function of redshift (equal weights).
medianMass = np.zeros(redshiftCount)
for iBin in range(redshiftCount):
    if iBin < redshiftCount - 1:
        inBin = np.where((redshift >= redshiftBins[iBin]) & (redshift < redshiftBins[iBin + 1]))[0]
    else:
        inBin = np.where(redshift >= redshiftBins[iBin])[0]
    if len(inBin) > 0:
        medianMass[iBin] = np.median(logMass[inBin])

# Generate a fit to the data (fit redshift as function of medianMass, degree 2).
nonZero = np.where(medianMass > 0.0)[0]
coeffsDescending = np.polyfit(medianMass[nonZero], redshiftBins[nonZero], 2)
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
ax.set_title('Redshift vs. Limiting Stellar Mass for Caputi et al. (2011) Sample')
ax.set_xlabel(r'Limiting stellar mass; $M_\star\;[M_\odot]$')
ax.set_ylabel(r'Redshift; $z$')
ax.set_xscale('log')
ax.set_xlim(6.0e9, 6.0e10)
ax.set_ylim(2.5, 5.5)
ax.plot(10.0 ** fitMass, fitRedshift,
        color='mediumseagreen', linewidth=2.5, label='Fit')
ax.plot(10.0 ** medianMass, redshiftBins,
        'o', color='goldenrod', markersize=8, label='Guo et al. (2011) SAM')
ax.legend(loc='lower left')
plt.tight_layout()
plotFile = workDirectory + 'massLuminosityRelation.pdf'
plt.savefig(plotFile)
plt.close()
