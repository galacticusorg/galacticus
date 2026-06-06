#!/usr/bin/env python3
# Fit polynomials to the mass and redshift dependence of stellar mass errors in the PRIMUS survey fields.
# Andrew Benson (22-May-2014)

import os
import sys
import numpy as np
import h5py
from astropy.io import fits

execPath = os.environ.get('GALACTICUS_EXEC_PATH', '')

# Get field solid angles.
with h5py.File(dataPath + '/static/surveyGeometry/PRIMUS/solidAngles.hdf5', 'r') as f:
    solidAngles = f['solidAngle'][:]
weights = solidAngles / solidAngles.sum()

# Define files for each field.
fileNames = [
    'cosmos.fits',
    'xmm_swire.fits',
    'cfhtls_xmm.fits',
    'cdfs.fits',
    'es1.fits',
]

fortranLines = []
latexLines   = []

# Sequence of exponents on mass and redshift.
sequence = [
    (0, 0),
    (1, 0),
    (2, 0),
    (0, 1),
    (0, 2),
    (1, 1),
]

for iField, fileName in enumerate(fileNames):
    fieldNum  = iField + 1
    filePath  = execPath + '/constraints/dataAnalysis/stellarMassFunctions_PRIMUS_z0_1/' + fileName
    if not os.path.exists(filePath):
        raise FileNotFoundError(
            f'massErrors.py: this script requires the PRIMUS data tables be available ({filePath})'
        )

    with fits.open(filePath) as hdul:
        table     = hdul[1].data
        zPrimus   = table['ZPRIMUS' ].astype(float)
        mass50    = table['MASS_50' ].astype(float)
        massErr   = table['MASS_ERR'].astype(float)

    # Construct matrices required to optimise the polynomial coefficients.
    nSeq = len(sequence)
    b = np.zeros(nSeq)
    a = np.zeros((nSeq, nSeq))
    for i, (mi, zi) in enumerate(sequence):
        b[i] = np.sum(massErr * mass50**mi * zPrimus**zi)
        for j, (mj, zj) in enumerate(sequence):
            a[i, j] = np.sum(mass50**mi * zPrimus**zi * mass50**mj * zPrimus**zj)

    # Solve for polynomial coefficients.
    x = np.linalg.solve(a, b)

    # Format field name.
    fieldName = fileName.upper().replace('_BENSON.FITS', '').replace('_', ' ')

    fortranLines.append(f'   ! Error fit for field: {fieldName}\n')
    fortranLines.append(f'   error({fieldNum})=                                             &\n')
    for i, (mi, zi) in enumerate(sequence):
        coeff        = x[i]
        pCoeff       = f'{coeff:16.10e}'.replace('e', 'd')
        prefix       = '+' if coeff >= 0.0 else ''
        massTerm     = f'*logMass**{mi}' if mi > 1 else ('*logMass   ' if mi == 1 else '           ')
        redshiftTerm = f'*redshift**{zi}' if zi > 1 else ('*redshift   ' if zi == 1 else '            ')
        fortranLines.append(f'      &         {prefix}{pCoeff}{massTerm}{redshiftTerm}')
        if i < len(sequence) - 1:
            fortranLines.append(' & \n')
        else:
            fortranLines.append('\n')

    latexLine = fieldName
    for i in range(nSeq):
        coeff  = x[i]
        pCoeff = f'{coeff:6.4f}'
        prefix = '+' if coeff >= 0.0 else ''
        latexLine += f' & {prefix}{pCoeff}'
    latexLines.append(latexLine + ' \\\\\n')

with open('constraints/dataAnalysis/stellarMassFunctions_PRIMUS_z0_1/massErrors.F90', 'w') as f:
    f.writelines(fortranLines)

with open('constraints/dataAnalysis/stellarMassFunctions_PRIMUS_z0_1/massErrors.tex', 'w') as f:
    f.writelines(latexLines)
