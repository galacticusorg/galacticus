#!/usr/bin/env python3
# Determine the relation between stellar mass and limiting distance for the VIPERS stellar mass functions.
# Andrew Benson (04-June-2014)

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from astropy.cosmology import LambdaCDM
import astropy.units as u

execPath = os.environ.get('GALACTICUS_EXEC_PATH', '').rstrip('/') + '/'
dataPath = os.environ.get('GALACTICUS_DATA_PATH', '').rstrip('/') + '/'

# Define a working directory.
workDirectory = execPath + 'constraints/dataAnalysis/stellarMassFunctions_VIPERS_z0_1/'

# Define survey solid angle (computed from mangle polygons).
solidAngle = 0.003137

# Read the stellar mass function.
xmlFile = dataPath + 'static/observations/massFunctionsStellar/Stellar_Mass_Functions_VIPERS_2013.xml'
tree = ET.parse(xmlFile)
root = tree.getroot()

cosmoElem    = root.find('cosmology')
omegaMatter  = float(cosmoElem.find('omegaMatter'    ).text)
omegaDE      = float(cosmoElem.find('omegaDarkEnergy').text)
hubble       = float(cosmoElem.find('hubble'         ).text)

cosmologyObserved = LambdaCDM(H0=hubble, Om0=omegaMatter, Ode0=omegaDE)

# Estimated sampling rates for each redshift.
samplingRate = np.array([0.44, 0.44, 0.4])

# Create a plot.
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_title('Limiting Distance vs. Stellar Mass for VIPERS Survey')
ax.set_xlabel(r'Stellar mass; $M_\star\;[M_\odot]$')
ax.set_ylabel(r'Limiting distance; $D_{\rm max}$ [Mpc]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(3.0e9, 3.0e12)
ax.set_ylim(2000.0, 4000.0)

colors1 = ['steelblue', 'darkorange', 'green']
colors2 = ['navy',      'saddlebrown', 'darkgreen']

massFunction = root.find('massFunction')
columns = massFunction.findall('columns')
for iMF, column in enumerate(columns):
    redshiftBin = iMF

    logarithmicMass = np.array([float(d.text) for d in column.find('mass'        ).findall('datum')])
    massFunction    = np.array([float(d.text) for d in column.find('massFunction').findall('datum')])
    number          = np.array([float(d.text) for d in column.find('number'      ).findall('datum')])
    binWidth        = np.full(len(massFunction), 0.2)
    redshiftMinimum = float(column.find('redshiftLow' ).text)
    redshiftMaximum = float(column.find('redshiftHigh').text)
    redshiftLabel   = f'${redshiftMinimum:.1f}<z<{redshiftMaximum:.1f}$'

    # Find distance to minimum and maximum redshifts.
    distanceMinimum = cosmologyObserved.comoving_distance(redshiftMinimum).to(u.Mpc).value
    distanceMaximum = cosmologyObserved.comoving_distance(redshiftMaximum).to(u.Mpc).value

    # Compute the effective volume.
    volume = number / samplingRate[redshiftBin] / massFunction / binWidth

    # Find the maximum distance.
    distance            = (3.0 * volume / solidAngle + distanceMinimum**3) ** (1.0 / 3.0)
    logarithmicDistance = np.log10(distance)

    # Fit a polynomial to the results.
    fitTo = np.where(
        (massFunction > 0.0) &
        (distance < 0.98 * distanceMaximum) &
        (logarithmicMass < 12)
    )[0]
    coeffsDescending = np.polyfit(logarithmicMass[fitTo], logarithmicDistance[fitTo], 2)
    coeffsAscending  = coeffsDescending[::-1]

    fitMass = np.linspace(4.0, 13.0, 1000)
    fitDistance = np.zeros(len(fitMass))
    closing = ''
    print(f'{redshiftMinimum} < z < {redshiftMaximum}: ', end='')
    for i, c in enumerate(coeffsAscending):
        fitDistance += c * fitMass ** i
        print(f'{c:.10e}d0', end='')
        if i == len(coeffsAscending) - 1:
            print(closing)
        else:
            print('+logarithmicMass*(', end='')
            closing += ')'

    ax.plot(10.0 ** fitMass, 10.0 ** fitDistance,
            color=colors2[redshiftBin], linewidth=2.5, label=f'Fit [{redshiftLabel}]')
    ax.plot(10.0 ** logarithmicMass, 10.0 ** logarithmicDistance,
            'o', color=colors1[redshiftBin], markersize=6, label=f'VIPERS [{redshiftLabel}]')

ax.legend(loc='upper left')
plt.tight_layout()
plotFile = workDirectory + 'massDistanceRelation.pdf'
plt.savefig(plotFile)
plt.close()
