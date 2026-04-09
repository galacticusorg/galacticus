#!/usr/bin/env python3
# Construct a mass error model for the ALFALFA survey.
# Andrew Benson (10-July-2013)

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

# Define working directory.
workDirectory = 'constraints/dataAnalysis/hiMassFunction_ALFALFA_z0.00/'

# Load the data.
tree = ET.parse(workDirectory + 'alfalfaHIMassErrorModel.xml')
root = tree.getroot()
logarithmicMass      = np.array([float(d.text) for d in root.find('mass'     ).findall('datum')])
logarithmicMassError = np.array([float(d.text) for d in root.find('massError').findall('datum')])

# Generate a fit to the data. Functional form: a + exp(-(logMass - b) / c)
aBest          = 0.100
bBest          = 5.885
cBest          = 0.505
fitMeasureBest = 1.0e30

aVals = np.arange(0.05, 0.155, 0.005)
bVals = np.arange(5.00, 7.005, 0.005)
cVals = np.arange(0.25, 0.755, 0.005)

for a in aVals:
    for b in bVals:
        for c in cVals:
            fit        = a + np.exp(-(logarithmicMass - b) / c)
            fitMeasure = np.sum((fit - logarithmicMassError)**2)
            if fitMeasure < fitMeasureBest:
                aBest          = a
                bBest          = b
                cBest          = c
                fitMeasureBest = fitMeasure

print('Parameters of best fitting model:')
print(f'  a = {aBest}')
print(f'  b = {bBest}')
print(f'  c = {cBest}')

massFine = np.linspace(5.0, 12.0, 1000)
fitFine  = aBest + np.exp(-(massFine - bBest) / cBest)

# Create a plot.
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xlabel(r'HI gas mass; $M_{\rm HI}\,[{\rm M}_\odot]$')
ax.set_ylabel(r'Error in $\log_{10}$ of HI gas mass; $\sigma$')
ax.set_xscale('log')
ax.set_xlim(3.1e5, 3.1e11)
ax.set_ylim(0.0, 0.8)
ax.plot(10.0 ** logarithmicMass, logarithmicMassError,
        'o', color='peachpuff', markersize=8, markeredgecolor='gray',
        label='Haynes et al. (2011)')
ax.plot(10.0 ** massFine, fitFine,
        color='mediumseagreen', linewidth=2.5, label='Fit')
ax.legend(loc='upper right')
plt.tight_layout()
plotFile = workDirectory + 'alfalfaHIMassErrorModel.pdf'
plt.savefig(plotFile)
plt.close()
