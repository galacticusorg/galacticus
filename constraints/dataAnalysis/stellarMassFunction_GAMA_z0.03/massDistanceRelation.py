#!/usr/bin/env python3
# Determine the relation between stellar mass and limiting distance for the GAMA stellar mass function.
# Andrew Benson (04-June-2014)

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

execPath = os.environ.get('GALACTICUS_EXEC_PATH', '').rstrip('/') + '/'
dataPath = os.environ.get('GALACTICUS_DATA_PATH', '').rstrip('/') + '/'

# Define a working directory.
workDirectory = execPath + 'constraints/dataAnalysis/stellarMassFunction_GAMA_z0.03/'

# Define survey solid angle.
squareDegreesPerSteradian = 3282.80635
solidAngle = 143.0 / squareDegreesPerSteradian

# Read the stellar mass function.
xmlFile = dataPath + 'static/observations/massFunctionsStellar/Stellar_Mass_Function_GAMA_2012.xml'
tree = ET.parse(xmlFile)
root = tree.getroot()
columns = root.find('massFunction').find('columns')
logarithmicMass = np.array([float(x.text) for x in columns.find('mass'        ).findall('datum')], dtype=float)
massFunction    = np.array([float(x.text) for x in columns.find('massFunction').findall('datum')], dtype=float)
number          = np.array([float(x.text) for x in columns.find('number'      ).findall('datum')], dtype=float)
binWidth        = np.full(len(massFunction), 0.2)
binWidth[0:2]   = 0.5

# Compute the effective volume.
volume = number / massFunction / binWidth

# Define the fields.
fields = [
    {'label': 'G09/G15', 'depth': 19.4, 'fieldCount': 2},
    {'label': 'G12',     'depth': 19.8, 'fieldCount': 1},
]

# Find relative volumes of the fields.
volumeTotal = 0.0
for field in fields:
    field['distance'] = 10.0 ** (0.4 * (field['depth'] - 19.0))
    field['volume'  ] = field['distance'] ** 3
    volumeTotal += field['fieldCount'] * field['volume']
for field in fields:
    field['volume'] /= volumeTotal

# Create a plot.
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_title('Limiting Distance vs. Stellar Mass for GAMA Survey')
ax.set_xlabel(r'Stellar mass; $M_\star\;[M_\odot]$')
ax.set_ylabel(r'Limiting distance; $D_{\rm max}$ [Mpc]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1.0e6, 1.0e12)
ax.set_ylim(20.0, 500.0)

colors = ['steelblue', 'darkorange']
for iField, field in enumerate(fields):
    # Find the maximum distance.
    distance            = (3.0 * volume * field['volume'] / (solidAngle / 3.0)) ** (1.0 / 3.0)
    logarithmicDistance = np.log10(distance)

    # Fit a polynomial to the results (degree 2, only for logMass < 9).
    fitTo   = np.where(logarithmicMass < 9.0)[0]
    coeffsAscending = np.polyfit(logarithmicMass[fitTo], logarithmicDistance[fitTo], 2)[::-1]

    # Generate a fit to the data.
    fitMass = np.linspace(6.0, 12.0, 1000)
    fitDistance = np.zeros(len(fitMass))
    closing = ''
    print(field['label'] + ':\t', end='')
    for i, c in enumerate(coeffsAscending):
        fitDistance += c * fitMass ** i
        print(f'{c:.10e}d0', end='')
        if i == len(coeffsAscending) - 1:
            print(closing)
        else:
            print('+logarithmicMass*(', end='')
            closing += ')'

    ax.plot(10.0 ** fitMass, 10.0 ** fitDistance,
            color=colors[iField], linewidth=2, label=f'Fit [{field["label"]}]')
    ax.plot(10.0 ** logarithmicMass, 10.0 ** logarithmicDistance,
            'o', color=colors[iField], markersize=6, label=f'GAMA [{field["label"]}]')

ax.legend(loc='lower right')
plt.tight_layout()
plotFile = workDirectory + 'massDistanceRelation.pdf'
plt.savefig(plotFile)
plt.close()
