#!/usr/bin/env python3
# Fit stellar mass completeness limits from the ULTRAVISTA survey (Muzzin et al. 2013). Fits the tabulated results in the data
# file downloaded from the ULTRAVISTA web site.
# Andrew Benson (13-August-2014)

import os
import sys
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

execPath = os.environ.get('GALACTICUS_EXEC_PATH', '')

def galacticusPath():
    return execPath.rstrip('/') + '/'

# Specify work directory.
workDirectory = galacticusPath() + 'constraints/dataAnalysis/stellarMassFunctions_ULTRAVISTA_z0.2_4.0/'

# Get the data file if we do not already have it.
dataFile = workDirectory + 'Mstar_redshift_completeness_emp_uvista_v4.1_95.dat'
if not os.path.exists(dataFile):
    subprocess.run([
        'wget',
        'http://www.strw.leidenuniv.nl/galaxyevolution/ULTRAVISTA/Mstar_redshift_completeness_emp_uvista_v4.1_95.dat',
        '-O', dataFile
    ], check=True)

# Read the data file.
data     = np.loadtxt(dataFile, usecols=(0, 1))
redshift = data[:, 0]
mass     = data[:, 1]

# Use only points with redshift < 3.
unique = np.where(redshift < 3.0)[0]

# Create a plot.
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_title(r"Redshift vs. Limiting Stellar Mass for Muzzin et al. (2013) Sample")
ax.set_xlabel(r'Limiting stellar mass; $M_\star\;[M_\odot]$')
ax.set_ylabel(r'Redshift; $z$')
ax.set_xscale('log')
ax.set_xlim(3.0e8, 3.0e11)
ax.set_ylim(0.2, 4.0)

# Fit a degree-6 polynomial.
coeffsDescending = np.polyfit(mass[unique], redshift[unique], 6)
coeffsAscending  = coeffsDescending[::-1]
fitMass = np.linspace(7.0, 12.0, 1000)
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

# Apply suppression factor for high masses.
fitRedshift *= 1.0 / (1.0 - np.exp((fitMass - 11.24) / 0.02))

ax.plot(10.0 ** fitMass, fitRedshift,
        color='steelblue', linewidth=1.5, label='fit')
ax.plot(10.0 ** mass, redshift, '--',
        color='steelblue', linewidth=2.5, label='observed')

ax.legend(loc='lower right')
plt.tight_layout()
plotFile = workDirectory + 'massRedshiftRelation.pdf'
plt.savefig(plotFile)
plt.close()
