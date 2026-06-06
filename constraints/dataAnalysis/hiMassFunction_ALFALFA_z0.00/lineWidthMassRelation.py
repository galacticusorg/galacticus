#!/usr/bin/env python3
# Construct the median line-width vs. HI mass relation for the ALFALFA survey.
# Andrew Benson (9-July-2013)

import os
import sys
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Define working directory.
workDirectory = 'constraints/dataAnalysis/hiMassFunction_ALFALFA_z0.00/'

# Retrieve the ALFALFA 40% dataset if we don't already have it.
csvRaw   = workDirectory + 'alfalfaSurveyData.csv'
csvClean = workDirectory + 'alfalfaSurveyDataClean.csv'
if not os.path.exists(csvRaw):
    subprocess.run([
        'wget',
        'http://egg.astro.cornell.edu/alfalfa/data/a40files/a40.datafile1.csv',
        '-O', csvRaw
    ], check=True)

# Remove troublesome double-double quotes.
if not os.path.exists(csvClean):
    with open(csvRaw, 'r') as inFile, open(csvClean, 'w') as outFile:
        for line in inFile:
            line = line.replace('""', '')
            outFile.write(line)

# Load the data (skip header line starting with 'AGCNr').
W50                = []
errorW50           = []
logarithmicMassHI  = []
with open(csvClean, 'r') as f:
    for line in f:
        if line.startswith('AGCNr'):
            continue
        parts = line.strip().split(',')
        if len(parts) > 14:
            try:
                W50_               = float(parts[ 7])
            except ValueError:
                continue
            try:
                errorW50_          = float(parts[ 8])
            except ValueError:
                continue
            try:
                logarithmicMassHI_ = float(parts[14])
            except ValueError:
                continue
            W50              .append(W50_              )
            errorW50         .append(errorW50_         )
            logarithmicMassHI.append(logarithmicMassHI_)

W50               = np.array(W50              )
errorW50          = np.array(errorW50         )
logarithmicMassHI = np.array(logarithmicMassHI)

# Find usable galaxies.
usable = np.where((logarithmicMassHI > 0.0) & (errorW50 > 0.0))[0]

# Bin the data and find the weighted median line width vs. logarithmic mass.
logMassUsable = logarithmicMassHI[usable]
w50Usable     = np.log10(W50[usable])
errUsable     = errorW50[usable]
weight        = 1.0 / errUsable**2

sortIdx           = np.argsort(logMassUsable)
logMassUsable     = logMassUsable[sortIdx]
w50Usable         = w50Usable    [sortIdx]
weight            = weight       [sortIdx]

logMassMin  = logMassUsable[ 0]
logMassMax  = logMassUsable[-1]
nBins       = len(usable) // 500
binEdges    = np.linspace(logMassMin, logMassMax, nBins)

# Compute weighted median in each bin.
medianW50 = np.zeros(nBins)
for iBin in range(nBins):
    if iBin < nBins - 1:
        inBin = np.where((logMassUsable >= binEdges[iBin]) & (logMassUsable < binEdges[iBin + 1]))[0]
    else:
        inBin = np.where(logMassUsable >= binEdges[iBin])[0]
    if len(inBin) == 0:
        continue
    w   = weight   [inBin]
    val = w50Usable[inBin]
    idx = np.argsort(val)
    cumW = np.cumsum(w[idx])
    totalW = cumW[-1]
    medIdx = np.searchsorted(cumW, 0.5 * totalW)
    medianW50[iBin] = val[idx[medIdx]]

# Fit a degree-2 polynomial to the median relation.
nonZero = np.where(medianW50 > 0.0)[0]
coeffsDescending = np.polyfit(binEdges[nonZero], medianW50[nonZero], 2)
coeffsAscending  = coeffsDescending[::-1]
fitW50 = np.zeros(len(binEdges))
closing = ''
for i, c in enumerate(coeffsAscending):
    fitW50 += c * binEdges ** i
    print(f'{c:.10e}d0', end='')
    if i == len(coeffsAscending) - 1:
        print(closing)
    else:
        print('+logarithmicMass*(', end='')
        closing += ')'

# Create a plot.
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xlabel(r'HI gas mass; $M_{\rm HI}\;[{\rm M}_\odot]$')
ax.set_ylabel(r'Line width; $W_{50}$ [km/s]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1.0e8, 1.0e11)
ax.set_ylim(10.0, 600.0)

# Plot a random 1% subsample of the data.
rng      = np.random.default_rng(seed=0)
plotSel  = rng.random(len(usable)) < 0.01
ax.errorbar(
    10.0 ** logarithmicMassHI[usable][plotSel],
    W50[usable][plotSel],
    yerr=errorW50[usable][plotSel],
    fmt='o', color='goldenrod', markersize=1.5, linewidth=0.5,
    label=r'$\alpha.40$'
)
ax.plot(10.0 ** binEdges, 10.0 ** medianW50,
        'o', color='peachpuff', markersize=6, markeredgecolor='gray', label='Median')
ax.plot(10.0 ** binEdges, 10.0 ** fitW50,
        color='mediumseagreen', linewidth=2.5, label='Fit')

ax.legend(loc='upper left')
plt.tight_layout()
plotFile = workDirectory + 'lineWidthMassRelation.pdf'
plt.savefig(plotFile)
plt.close()
