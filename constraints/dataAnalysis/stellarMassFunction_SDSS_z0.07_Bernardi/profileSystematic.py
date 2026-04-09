#!/usr/bin/env python3
# Compute mass systematic model coefficients to describe the mass systematic in the SDSS stellar
# mass function arising from the choice of profile fitting. Based on the results of Bernardi et
# al. (2013; http://adsabs.harvard.edu/abs/2013MNRAS.436..697B).
# Andrew Benson (07-April-2014)

import os
import sys
import numpy as np
from scipy.special import gamma as gamma_func
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Specify the various models and mass function fits used by Bernardi et al. (2013).
models = [
    {
        'name'    : 'cmodel',
        'phiStar' : 0.766e-2,
        'mStar'   : 0.4103e9,
        'alpha'   : 1.764   ,
        'beta'    : 0.384   ,
        'phiGamma': 0.557e-2,
        'mGamma'  : 4.7802e9,
        'gamma'   : 0.053
    },
    {
        'name'    : "Sérsic",
        'phiStar' : 1.040e-2,
        'mStar'   : 0.0094e9,
        'alpha'   : 1.665   ,
        'beta'    : 0.255   ,
        'phiGamma': 0.675e-2,
        'mGamma'  : 2.7031e9,
        'gamma'   : 0.296
    },
    {
        'name'    : "SérExp",
        'phiStar' : 0.892e-2,
        'mStar'   : 0.0014e9,
        'alpha'   : 2.330   ,
        'beta'    : 0.239   ,
        'phiGamma': 0.738e-2,
        'mGamma'  : 3.2324e9,
        'gamma'   : 0.305
    },
    {
        'name'    : "Sérsic (Simard)",
        'phiStar' : 0.820e-2,
        'mStar'   : 0.0847e9,
        'alpha'   : 1.755   ,
        'beta'    : 0.310   ,
        'phiGamma': 0.539e-2,
        'mGamma'  : 5.2204e9,
        'gamma'   : 0.072
    },
]

# Initialize suitable ranges of masses.
logMassLimited    = np.linspace(9.0, 12.5, 36)
massLimited       = 10.0 ** logMassLimited
logMassNormalized = logMassLimited - 10.8
logMass           = np.linspace(8.0, 13.5, 56)
mass              = 10.0 ** logMass

# Compute mass functions for each model.
for model in models:
    gam_val = gamma_func(model['alpha'] / model['beta'])
    model['phi'] = (
        model['phiStar'] * model['beta'] *
        (mass / model['mStar'])**model['alpha'] *
        np.exp(-(mass / model['mStar'])**model['beta']) / gam_val +
        model['phiGamma'] * (mass / model['mGamma'])**model['gamma'] *
        np.exp(-mass / model['mGamma'])
    )
    if model['name'] == 'cmodel':
        model['phiLimited'] = (
            model['phiStar'] * model['beta'] *
            (massLimited / model['mStar'])**model['alpha'] *
            np.exp(-(massLimited / model['mStar'])**model['beta']) / gam_val +
            model['phiGamma'] * (massLimited / model['mGamma'])**model['gamma'] *
            np.exp(-massLimited / model['mGamma'])
        )
        model['logMassLimited'] = logMassLimited
    else:
        # Abundance matching: interpolate to find mass in this model corresponding to cmodel phi.
        model['logMassLimited'] = np.interp(models[0]['phiLimited'], model['phi'], logMass)
    model['offset'] = models[0]['logMassLimited'] - model['logMassLimited']

# Create a plot.
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xlim(9.0, 12.5)
ax.set_ylim(-1.0, 1.0)
ax.set_xlabel(r'$\log_{10} (M_{\star,{\tt cmodel}}/{\rm M}_\odot)$')
ax.set_ylabel(r'$\log_{10} M_{\star,{\tt cmodel}}-\log_{10} M_{\star}$')

colors = ['steelblue', 'darkorange', 'green']
iColor = -1
for model in models:
    if model['name'] != 'cmodel':
        iColor += 1
        ax.plot(logMassLimited, model['offset'],
                'o', color=colors[iColor], markersize=6, label=model['name'])

# Construct a simple fit to the offsets.
fits = [
    {'mu0': -0.10, 'kappa0': -0.00, 'mu1': -0.00, 'kappa1': -0.33, 'beta': +0.50},
    {'mu0': -0.10, 'kappa0': -0.00, 'mu1': -0.00, 'kappa1': -0.25, 'beta': +0.50},
]
fitColors = ['red', 'purple']
for iFit, fit in enumerate(fits):
    systematic = (
        (fit['mu0'] + fit['kappa0'] * logMassNormalized) * (1.0 - 1.0 / (1.0 + np.exp(-logMassNormalized / fit['beta']))) +
        (fit['mu1'] + fit['kappa1'] * logMassNormalized) * (       1.0 / (1.0 + np.exp(-logMassNormalized / fit['beta'])))
    )
    ax.plot(logMassLimited, systematic,
            color=fitColors[iFit], linewidth=2.5, label='model')

ax.legend(loc='lower left')
plt.tight_layout()
plotFile = 'constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07_Bernardi/profileSystematic.pdf'
plt.savefig(plotFile, bbox_inches='tight')
plt.close()
