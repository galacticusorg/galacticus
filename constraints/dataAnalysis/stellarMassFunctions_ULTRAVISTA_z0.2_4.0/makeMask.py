#!/usr/bin/env python3
# Construct a mangle mask file for the UltraVISTA survey of Muzzin et al. (2013; http://adsabs.harvard.edu/abs/2013ApJ...777...18M)
# Andrew Benson (27-August-2014)

import os
import sys
import subprocess
import numpy as np

execPath = os.environ.get('GALACTICUS_EXEC_PATH', '').rstrip('/') + '/'
dataPath = os.environ.get('GALACTICUS_DATA_PATH', '').rstrip('/') + '/'

# Define data directory.
dataDirectoryName = execPath + 'constraints/dataAnalysis/stellarMassFunctions_ULTRAVISTA_z0.2_4.0/'
# Define work directory.
workDirectoryName = dataDirectoryName + 'work/'
# Define mangle directory.
mangle = dataPath + 'dynamic/mangle-2.3.3/bin/'

os.makedirs(workDirectoryName, exist_ok=True)

# Generate raw field geometry.
fieldFile = dataDirectoryName + 'field.ply'
if not os.path.exists(fieldFile):
    with open(fieldFile, 'w') as f:
        f.write('149.373 150.779 1.604 2.81\n')

# Specify radii for bright and medium stars.
brightRadius = 75.0 / 3600.0
mediumRadius = 30.0 / 3600.0

# Generate star mask holes.
starsFile = dataDirectoryName + 'stars.ply'
if not os.path.exists(starsFile):
    with open(starsFile, 'w') as starFile:
        # Read USNO star list.
        usnoData = np.loadtxt(dataDirectoryName + 'USNO_star_list',
                              usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
        id_usno   = usnoData[:, 0]
        ra_usno   = usnoData[:, 1]
        dec_usno  = usnoData[:, 2]
        b1        = usnoData[:, 3]

        # Find USNO bright and medium stars.
        usnoBright = np.where((b1 >  0.0) & (b1 <= 10.0))[0]
        usnoMedium = np.where((b1 > 10.0) & (b1 <= 13.0))[0]

        for i in usnoBright:
            starFile.write(f'{ra_usno[i]}\t{dec_usno[i]}\t{brightRadius}\n')
        for i in usnoMedium:
            starFile.write(f'{ra_usno[i]}\t{dec_usno[i]}\t{mediumRadius}\n')

        # Read 2MASS stars list (skip lines beginning with | or \).
        twoMassData = []
        with open(dataDirectoryName + '2mass_psc_list.dat', 'r') as tmf:
            for line in tmf:
                if line.startswith('|') or line.startswith('\\'):
                    continue
                parts = line.split()
                if len(parts) >= 12:
                    twoMassData.append([float(p) for p in parts[:12]])
        if twoMassData:
            twoMassArr   = np.array(twoMassData)
            ra_2mass     = twoMassArr[:, 0]
            dec_2mass    = twoMassArr[:, 1]
            k_2mass      = twoMassArr[:, 8]

            twomassBright = np.where((k_2mass >  0.0) & (k_2mass <=  8.0))[0]
            twomassMedium = np.where((k_2mass >  8.0) & (k_2mass <= 10.5))[0]

            for i in twomassBright:
                starFile.write(f'{ra_2mass[i]}\t{dec_2mass[i]}\t{brightRadius}\n')
            for i in twomassMedium:
                starFile.write(f'{ra_2mass[i]}\t{dec_2mass[i]}\t{mediumRadius}\n')

# Generate bad pixel map.
uvistaBadPly = workDirectoryName + 'UVISTA_Ks_15_12_10_bad.ply'
uvistaWeightFits = workDirectoryName + 'UVISTA_Ks_15_12_10_skysub_015_v1.weight.fits'
if not os.path.exists(uvistaWeightFits):
    subprocess.run([
        'wget',
        'http://ultravista.org/release1/data/UVISTA_Ks_15_12_10_skysub_015_v1.weight.fits',
        '-O', uvistaWeightFits
    ], check=True)
if not os.path.exists(uvistaBadPly):
    subprocess.run(['python3', dataDirectoryName + 'badPixels.py', workDirectoryName], check=True)

def run_mangle(cmd):
    subprocess.run(cmd, check=True)

# Pixelize the survey field.
if not os.path.exists(dataDirectoryName + 'fieldP.ply'):
    run_mangle([mangle + 'pixelize', '-Ps0,11', '-ir',
                dataDirectoryName + 'field.ply',
                dataDirectoryName + 'fieldP.ply'])

# Pixelize the stars.
if not os.path.exists(dataDirectoryName + 'starsP.ply'):
    run_mangle([mangle + 'pixelize', '-Ps0,11', '-ic',
                dataDirectoryName + 'stars.ply',
                dataDirectoryName + 'starsP.ply'])

# Pixelize the bad pixels.
if not os.path.exists(workDirectoryName + 'UVISTA_Ks_15_12_10_badP.ply'):
    run_mangle([mangle + 'pixelize', '-Ps0,11', '-ir',
                workDirectoryName + 'UVISTA_Ks_15_12_10_bad.ply',
                workDirectoryName + 'UVISTA_Ks_15_12_10_badP.ply'])

# Snap the survey field.
if not os.path.exists(dataDirectoryName + 'fieldPS.ply'):
    run_mangle([mangle + 'snap',
                dataDirectoryName + 'fieldP.ply',
                dataDirectoryName + 'fieldPS.ply'])

# Snap the stars.
if not os.path.exists(dataDirectoryName + 'starsPS.ply'):
    run_mangle([mangle + 'snap',
                dataDirectoryName + 'starsP.ply',
                dataDirectoryName + 'starsPS.ply'])

# Snap the bad pixels.
if not os.path.exists(workDirectoryName + 'UVISTA_Ks_15_12_10_badPS.ply'):
    run_mangle([mangle + 'snap',
                '-a1.0e-7d', '-b1.0e-7d', '-t1.0e-7d',
                workDirectoryName + 'UVISTA_Ks_15_12_10_badP.ply',
                workDirectoryName + 'UVISTA_Ks_15_12_10_badPS.ply'])

# Holeize the stars.
holesFile = dataDirectoryName + 'holesPS.ply'
if not os.path.exists(holesFile):
    with open(dataDirectoryName + 'starsPS.ply', 'r') as starsIn, \
         open(holesFile, 'w') as holesOut:
        for line in starsIn:
            line = line.replace('1 weight', '0 weight')
            holesOut.write(line)

# Holeize the bad pixels.
badHolesFile = workDirectoryName + 'badHolesPS.ply'
if not os.path.exists(badHolesFile):
    with open(workDirectoryName + 'UVISTA_Ks_15_12_10_badPS.ply', 'r') as badPixelsIn, \
         open(badHolesFile, 'w') as holesOut:
        for line in badPixelsIn:
            line = line.replace('1 weight', '0 weight')
            holesOut.write(line)

# Balkanize the survey and holes.
maskFile = workDirectoryName + 'mask.ply'
if not os.path.exists(maskFile):
    run_mangle([mangle + 'balkanize',
                dataDirectoryName + 'fieldPS.ply',
                dataDirectoryName + 'holesPS.ply',
                badHolesFile,
                maskFile])

# Unify the mask.
surveyMaskFile = dataDirectoryName + 'surveyMask.ply'
if not os.path.exists(surveyMaskFile):
    run_mangle([mangle + 'unify', maskFile, surveyMaskFile])
