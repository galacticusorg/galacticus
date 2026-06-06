#!/usr/bin/env python3
# Fit stellar mass completeness limits from the ZFOURGE survey (Tomczak et al. 2014). Fits the tabulated results in the data file
# given by R. Quadri.
# Andrew Benson (11-August-2014)

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

execPath = os.environ.get('GALACTICUS_EXEC_PATH', '')

# Read data.
fields = [
    {'name': 'NMBS'    },
    {'name': 'ZFOURGE' },
]
data = np.loadtxt('constraints/dataAnalysis/stellarMassFunctions_ZFOURGE_z0.2_2.5/zfourge-SMF-supplemental/masslimits.dat', usecols=(0, 1, 2))
redshift     = data[:, 0]
massNMBS     = data[:, 1]
massZFOURGE  = data[:, 2]
fields[0]['mass'    ] = massNMBS
fields[1]['mass'    ] = massZFOURGE
fields[0]['redshift'] = redshift
fields[1]['redshift'] = redshift

# Create a plot.
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_title(r"Redshift vs. Limiting Stellar Mass for Tomczak et al. (2014) Sample")
ax.set_xlabel(r'Limiting stellar mass; $M_\star\;[M_\odot]$')
ax.set_ylabel(r'Redshift; $z$')
ax.set_xscale('log')
ax.set_xlim(3.0e7, 3.0e11)
ax.set_ylim(0.2, 3.1)

colors = ['steelblue', 'darkorange']
for iField, field in enumerate(fields):
    # Fit a degree-4 polynomial.
    coeffsDescending = np.polyfit(field['mass'], field['redshift'], 4)
    coeffsAscending  = coeffsDescending[::-1]
    fitMass = np.linspace(7.0, 12.0, 1000)
    fitRedshift = np.zeros(len(fitMass))
    closing = ''
    print(field['name'] + ' : ', end='')
    for i, c in enumerate(coeffsAscending):
        fitRedshift += c * fitMass ** i
        print(f'{c:.10e}d0', end='')
        if i == len(coeffsAscending) - 1:
            print(closing)
        else:
            print('+logarithmicMass*(', end='')
            closing += ')'
    ax.plot(10.0 ** fitMass, fitRedshift,
            color=colors[iField], linewidth=1.5, label=field['name'] + ' [fit]')
    ax.plot(10.0 ** field['mass'], field['redshift'], '--',
            color=colors[iField], linewidth=2.5, label=field['name'] + ' [observed]')

ax.legend(loc='lower right')
plt.tight_layout()
plotFile = 'constraints/dataAnalysis/stellarMassFunctions_ZFOURGE_z0.2_2.5/massRedshiftRelation.pdf'
plt.savefig(plotFile)
plt.close()
