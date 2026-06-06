#!/usr/bin/env python3
# Find relation between mass and maximum distance for the sample of Bernardi et al. (2013;
# http://adsabs.harvard.edu/abs/2013MNRAS.436..697B). Andrew Benson (24-April-2014)

import os
import sys
import subprocess
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

execPath = os.environ.get('GALACTICUS_EXEC_PATH', '').rstrip('/') + '/'
dataPath = os.environ.get('GALACTICUS_DATA_PATH', '').rstrip('/') + '/'

# Tabulated data provided by Mariangela Bernardi.
log10MassStellar = np.array([  # Stellar masses are in units of M_Solar.
     8.655397415,
     8.758553505,
     8.856801987,
     8.954696655,
     9.052925110,
     9.151350975,
     9.255033493,
     9.351667404,
     9.454163551,
     9.553510666,
     9.651016235,
     9.750605583,
     9.852743149,
     9.954152107,
    10.053371429,
    10.153795242,
    10.253108025,
    10.352433205,
    10.452750206,
    10.551837921,
    10.651668549,
    10.751879692,
    10.849768639,
    10.950572014,
    11.049318314,
    11.147963524,
    11.247805595,
    11.346734047,
    11.446274757,
    11.545268059,
    11.644490242,
    11.742430687,
    11.840696335,
    11.941007614,
    12.043179512,
    12.135819435,
])
volumeMaximum = np.array([  # Comoving volumes are in units of 10^9 Mpc^3.
     0.002602330,
     0.002865336,
     0.003689888,
     0.005141599,
     0.006540476,
     0.008913630,
     0.011576199,
     0.015452851,
     0.020556541,
     0.027066843,
     0.035638396,
     0.047748027,
     0.063546862,
     0.082955151,
     0.106538044,
     0.134253938,
     0.161168912,
     0.196770098,
     0.232025141,
     0.281275746,
     0.349784664,
     0.444777516,
     0.562816249,
     0.717034707,
     0.912019113,
     1.143364514,
     1.425395784,
     1.767075533,
     2.164316387,
     2.631458218,
     3.202964305,
     3.810579128,
     4.508823099,
     5.098926730,
     5.946447509,
     6.690833525,
])

# Solid angle of the sample.
solidAngle = None
polygonFile = execPath + 'constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07_Bernardi/sdss_dr72safe0_res6d.pol'
if not os.path.exists(polygonFile):
    subprocess.run([
        'wget',
        'https://zenodo.org/records/10998446/files/sdss_dr72safe0_res6d.pol.gz',
        '-O',  polygonFile + '.gz'
    ], check=True)
    subprocess.run([
        'gunzip',
        polygonFile + '.gz'
    ], check=True)
result = subprocess.run([
    dataPath + 'dynamic/mangle-2.3.3/bin/harmonize',
    execPath + 'constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07_Bernardi/sdss_dr72safe0_res6d.pol',
    '/dev/null'
], check=True, capture_output=True)
for line in result.stdout.splitlines():
    match = re.search(r'^area of \(weighted\) region is ([0-9\.]+) str/',line)
    if match:
        solidAngle = match.group(1)

# Convert volumes to maximum distance.
log10DistanceMaximum = np.log10((3.0 * volumeMaximum * 1.0e9 / solidAngle) ** (1.0 / 3.0))

# Fit a polynomial (degree 6).
coeffsDescending = np.polyfit(log10MassStellar, log10DistanceMaximum, 6)
coeffsAscending  = coeffsDescending[::-1]
fitMass = np.linspace(8.6, 12.2, 1000)
fitDistance = np.zeros(len(fitMass))
closing = ''
for i, c in enumerate(coeffsAscending):
    fitDistance += c * fitMass ** i
    print(f'{c:.10e}d0', end='')
    if i == len(coeffsAscending) - 1:
        print(closing)
    else:
        print('+logarithmicMass*(', end='')
        closing += ')'

# Create a plot.
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_title(r"Maximum Distance vs. Limiting Stellar Mass for Bernardi et al. (2013) Sample")
ax.set_xlabel(r'Limiting stellar mass; $M_\star\;[M_\odot]$')
ax.set_ylabel(r'Maximum distance; $D_{\rm max}$ [Mpc]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(3.0e8, 2.0e12)
ax.set_ylim(100.0, 2500.0)
ax.plot(10.0 ** log10MassStellar, 10.0 ** log10DistanceMaximum,
        'o', color='mediumseagreen', markersize=6, label='Bernardi et al. (2013)')
ax.plot(10.0 ** fitMass, 10.0 ** fitDistance,
        color='goldenrod', linewidth=2.5, label='fit')
ax.legend(loc='upper left')
plt.tight_layout()
plotFile = execPath + 'constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07_Bernardi/massDistanceRelation.pdf'
plt.savefig(plotFile)
plt.close()
