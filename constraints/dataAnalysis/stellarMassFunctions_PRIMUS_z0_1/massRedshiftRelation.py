#!/usr/bin/env python3
# Fit stellar mass completeness limits from the PRIMUS survey (Moustakas et al. 2013).
# Andrew Benson (24-April-2014)

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Completeness limit data for "All" galaxies taken from Table 2 of Moustakas et al. (2013).
fields = [
    {
        'name'    : 'COSMOS',
        'redshift': np.array([0.250, 0.350, 0.450, 0.575, 0.725, 0.900]),
        'mass'    : np.array([8.730, 9.140, 9.510, 9.920, 10.330, 10.710])
    },
    {
        'name'    : 'XMM-SXDS',
        'redshift': np.array([0.250, 0.350, 0.450, 0.575, 0.725, 0.900]),
        'mass'    : np.array([8.860, 9.230, 9.580, 9.970, 10.380, 10.780])
    },
    {
        'name'    : 'XMM-CFHTLS',
        'redshift': np.array([0.250, 0.350, 0.450, 0.575, 0.725, 0.900]),
        'mass'    : np.array([8.950, 9.230, 9.510, 9.870, 10.310, 10.830])
    },
    {
        'name'    : 'CDFS',
        'redshift': np.array([0.250, 0.350, 0.450, 0.575, 0.725, 0.900]),
        'mass'    : np.array([9.620, 9.870, 10.100, 10.370, 10.650, 10.940])
    },
    {
        'name'    : 'ELAIS-S1',
        'redshift': np.array([0.250, 0.350, 0.450, 0.575, 0.725, 0.900]),
        'mass'    : np.array([9.700, 9.990, 10.260, 10.560, 10.870, 11.170])
    },
]

# Create a plot.
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_title(r"Redshift vs. Limiting Stellar Mass for Moustakas et al. (2013) Sample")
ax.set_xlabel(r'Limiting stellar mass; $M_\star\;[M_\odot]$')
ax.set_ylabel(r'Redshift; $z$')
ax.set_xscale('log')
ax.set_xlim(3.0e8, 1.0e11)
ax.set_ylim(0.0, 1.0)

colors = ['steelblue', 'darkorange', 'green', 'red', 'purple']
for iField, field in enumerate(fields):
    # Fit a degree-3 polynomial.
    coeffsDescending = np.polyfit(field['mass'], field['redshift'], 3)
    coeffsAscending  = coeffsDescending[::-1]
    fitMass = np.linspace(8.0, 13.0, 1000)
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
            color=colors[iField], linewidth=2.5, label=field['name'] + ' [fit]')
    ax.plot(10.0 ** field['mass'], field['redshift'],
            'o', color=colors[iField], markersize=6, label=field['name'] + ' [PRIMUS]')

ax.legend(loc='upper left')
plt.tight_layout()
plotFile = 'constraints/dataAnalysis/stellarMassFunctions_PRIMUS_z0_1/massRedshiftRelation.pdf'
plt.savefig(plotFile)
plt.close()
