#!/usr/bin/env python3
# Estimate completeness as a function of mass for the UKIDSS UDS survey of Caputi et al. (2011).
# Andrew Benson (28-April-2014)

import os
import sys
import numpy as np
from scipy.special import erf, erfinv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from xml.dom import minidom

execPath = os.environ.get('GALACTICUS_EXEC_PATH', '').rstrip('/') + '/'
dataPath = os.environ.get('GALACTICUS_DATA_PATH', '').rstrip('/') + '/'

# Number of sigma above sky for detection (arbitrary - will scale out of the results).
sigma = 3.0

# Argument in error function required for 80% completeness.
# completeness = 0.5*(1 - erf(x/sqrt(2))); at 80%: erf(x/sqrt(2)) = -0.6, x = sqrt(2)*erfinv(-0.6)
y80 = np.sqrt(2.0) * erfinv(2.0 * (0.5 - 0.8))

# Specify mass completeness limits.
bins = [
    {'redshiftMinimum':  3.00, 'redshiftMaximum':  3.50, 'completeness50': 10.30, 'completeness80': 10.93},
    {'redshiftMinimum':  3.50, 'redshiftMaximum':  4.25, 'completeness50': 10.41, 'completeness80': 11.06},
    {'redshiftMinimum':  4.25, 'redshiftMaximum':  5.00, 'completeness50': 10.51, 'completeness80': 11.18},
]

# Load the existing data file.
dataFile = dataPath + 'static/observations/massFunctionsStellar/Stellar_Mass_Functions_UKIDSS_UDS_2011.xml'
tree = ET.parse(dataFile)
root = tree.getroot()

# Begin constructing the plot.
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xlabel(r'Galaxy stellar mass; $M_\star\;[M_\odot]$')
ax.set_ylabel('Completeness; []')
ax.set_xscale('log')
ax.set_xlim(1.0e10, 1.0e12)
ax.set_ylim(-0.05, 1.05)

colors = ['steelblue', 'darkorange', 'green']
massFunctions = root.findall('massFunction')

for iBin, binData in enumerate(bins):
    m50 = 10.0 ** binData['completeness50']
    m80 = 10.0 ** binData['completeness80']

    # Solve for required parameters in the completeness model.
    flux = y80**2 * m80 / ((m50 - m80)**2 - y80**2 * m50 / sigma)
    sky  = (flux * m50 / sigma)**2

    # Construct masses.
    logMass = np.linspace(10.0, 12.0, 1000)
    mass    = 10.0 ** logMass

    # Determine signal and noise as a function of mass.
    signal = flux * mass
    noise  = np.sqrt(sky)

    # Evaluate the completeness model.
    x            = (sigma * noise - signal) / np.sqrt(noise**2 + signal)
    completeness = 0.5 * (1.0 - erf(x / np.sqrt(2.0)))

    # Construct label.
    label = f'${binData["redshiftMinimum"]:.2f}<z<{binData["redshiftMaximum"]:.2f}$'

    ax.plot(mass, completeness, color=colors[iBin], linewidth=2.5, label=label)

    # Compute completeness in each observed mass bin.
    if iBin < len(massFunctions):
        mf = massFunctions[iBin]
        massCol = mf.find('columns').find('mass')
        if massCol is not None:
            observedMass         = np.array([float(d.text) for d in massCol.findall('datum')])
            observedSignal       = flux * 10.0 ** observedMass
            observedX            = (sigma * noise - observedSignal) / np.sqrt(noise**2 + observedSignal)
            observedCompleteness = 0.5 * (1.0 - erf(observedX / np.sqrt(2.0)))
            # Add/update completeness column.
            columns = mf.find('columns')
            compCol = columns.find('completeness')
            if compCol is None:
                compCol = ET.SubElement(columns, 'completeness')
            # Remove existing datum elements.
            for d in compCol.findall('datum'):
                compCol.remove(d)
            for val in observedCompleteness:
                d = ET.SubElement(compCol, 'datum')
                d.text = str(val)
            compCol.set('scaling',     'linear')
            compCol.set('description', 'Completeness in this mass bin.')

ax.legend(loc='lower left')
plt.tight_layout()
plotFile = execPath + 'constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/completeness.pdf'
plt.savefig(plotFile)
plt.close()

# Write the augmented data back to file.
xmlStr = minidom.parseString(ET.tostring(root, encoding='unicode')).toprettyxml(indent='  ')
with open(dataFile, 'w') as f:
    f.write(xmlStr)
