#!/usr/bin/env python3
"""
Test the log-likelihoods reported by the galaxy property vs. halo (or stellar) mass relation output analyses.

Each of these analyses evaluates a multivariate normal likelihood over a selected set of bins, via the shared
helpers in `source/output/analyses/utilities.F90`. Here a small model is run against synthetic target datasets
which exercise every branch of those helpers (the mean and the scatter of each relation; bins selected
explicitly, automatically, or not at all; and normalized and unnormalized likelihoods), and each reported
log-likelihood is compared with an independent calculation from the model and target datasets written to the
output file.
"""
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

import h5py
import numpy as np

# Paths below are relative to the root of the Galacticus source tree, so work from there: this test is run
# both from that root and, by `test-all.py` and the CI workflows, from the `testSuite` directory.
os.chdir(os.environ.get('GALACTICUS_EXEC_PATH', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')))

# Judged improbable by Galacticus (see `logImprobable` in `source/models/likelihoods/constants.F90`).
logImprobable  = -1.0e-16*np.finfo(np.float64).max
# Model values at or below these thresholds (log₁₀ mass, and log₁₀ radius in Mpc) indicate empty bins.
valueMinimum   = {'sizeVsStellarMassRelation': -6.0}
valueMinimumAll=  1.0e-3
toleranceRelative = 1.0e-9

executable    = 'Galacticus.exe'
parameterFile = 'testSuite/parameters/outputAnalysisLikelihoods.xml'
outputPath    = 'testSuite/outputs/outputAnalysisLikelihoods'
logFile       = outputPath+'/galacticus.log'

def targetsWrite():
    """Write synthetic target datasets for the relation analyses."""
    def cosmology(file, label):
        group = file.create_group('cosmology')
        group.attrs['OmegaMatter'    ] = 0.3153
        group.attrs['OmegaDarkEnergy'] = 0.6847
        group.attrs['OmegaBaryon'    ] = 0.0493
        group.attrs['HubbleConstant' ] = 67.36
        file.attrs['label'    ] = np.bytes_(label)
        file.attrs['reference'] = np.bytes_('Synthetic test data')
    def relationHalo(fileName, label, property, normalization, slope):
        massHalo = np.logspace(11.0, 14.0, 10)
        value    = 10.0**(normalization+slope*(np.log10(massHalo)-12.0))
        with h5py.File(fileName, 'w') as file:
            cosmology(file, label)
            file.attrs['haloMassDefinition'] = np.bytes_('virial')
            group = file.create_group('redshiftInterval1')
            group.attrs['redshiftMinimum'] = 0.22
            group.attrs['redshiftMaximum'] = 0.48
            group['massHalo'                 ] = massHalo
            group[property                   ] = value
            group[property+'Error'           ] = 0.3*value
            group[property+'Scatter'         ] = np.full(massHalo.size, 0.40)
            group[property+'ScatterError'    ] = np.full(massHalo.size, 0.05)
    def relationSize(fileName, label):
        massStellar = np.logspace(9.0, 11.5, 8)
        radius      = 3.0e-3*(massStellar/1.0e10)**0.25
        with h5py.File(fileName, 'w') as file:
            cosmology(file, label)
            group = file.create_group('sample1')
            group.attrs['redshiftMinimum'] = 0.22
            group.attrs['redshiftMaximum'] = 0.48
            group.attrs['selection'      ] = np.bytes_('none')
            group['massStellar'                ] = massStellar
            group['radiusEffective'            ] = radius
            group['radiusEffectiveError'       ] = 0.1*radius
            group['radiusEffectiveScatter'     ] = np.full(massStellar.size, 0.20)
            group['radiusEffectiveScatterError'] = np.full(massStellar.size, 0.03)
    for suffix in ('A', 'B'):
        relationHalo(outputPath+'/blackHole'+suffix+'.hdf5', 'SyntheticBH'  +suffix, 'massBlackHole',  7.0, 1.5)
        relationHalo(outputPath+'/stellar'  +suffix+'.hdf5', 'SyntheticSMHM'+suffix, 'massStellar'  , 10.3, 0.8)
        relationSize(outputPath+'/size'     +suffix+'.hdf5', 'SyntheticSize'+suffix                              )

def analysesRead():
    """Return the analyses requested in the parameter file, each with its target label."""
    analyses = []
    for analysis in ET.parse(parameterFile).getroot().find('outputAnalysis').findall('outputAnalysis'):
        options    = {child.tag: child.get('value') for child in analysis}
        targetFile = options.get('fileNameTarget')
        label      = None
        if targetFile is not None:
            with h5py.File(targetFile, 'r') as file:
                label = file.attrs['label'].decode()
        analyses.append((analysis.get('value'), label, options))
    return analyses

def analysisMatch(groupName, analyses):
    """Return the analysis (class name and options) whose results are stored in the named group."""
    isScatter = 'Scatter' in groupName
    for className, label, options in analyses:
        if (options.get('computeScatter') == 'true') != isScatter:
            continue
        if label is None:
            if className.endswith('Leauthaud2012') and 'Leauthaud' in groupName:
                return className, options
        elif label in groupName:
            return className, options
    return None, None

def logLikelihoodCompute(group, className, options):
    """Compute the log-likelihood of an analysis from the model and target datasets in its output group."""
    value            = group[group.attrs['yDataset'         ].decode()][:]
    covariance       = group[group.attrs['yCovariance'      ].decode()][:]
    valueTarget      = group[group.attrs['yDatasetTarget'   ].decode()][:]
    covarianceTarget = group[group.attrs['yCovarianceTarget'].decode()][:]
    # The Leauthaud2012 analysis compares the scatter to a fixed target, not to one read from file.
    if className.endswith('Leauthaud2012') and options.get('computeScatter') == 'true':
        valueTarget      = np.full(value.size, 0.16      )
        covarianceTarget = np.diag(np.full(value.size, 0.04**2))
    minimum  = valueMinimum.get(className, valueMinimumAll)
    binsText = options.get('likelihoodBins')
    if   binsText == 'auto':
        bins = np.nonzero(value != 0.0)[0]
    elif binsText is None:
        bins = np.array([], dtype=int)
    else:
        bins = np.array([int(bin)-1 for bin in binsText.split()], dtype=int)
    if (bins.size == 0 and np.any(value <= minimum)) or np.any(value[bins] <= minimum):
        return logImprobable
    selection  = bins if bins.size > 0 else np.arange(value.size)
    difference = (value-valueTarget)[selection]
    combined   = (covariance+covarianceTarget)[np.ix_(selection, selection)]
    try:
        cholesky = np.linalg.cholesky(combined)
    except np.linalg.LinAlgError:
        return logImprobable
    solution      = np.linalg.solve(cholesky, difference)
    logLikelihood = -0.5*np.dot(solution, solution)
    if options.get('likelihoodNormalize') == 'true':
        logLikelihood += -np.sum(np.log(np.diag(cholesky)))-0.5*difference.size*np.log(2.0*np.pi)
    return logLikelihood

# Check that prerequisites are met.
if not os.path.exists(executable):
    print("SKIPPED: "+executable+" does not exist - build it before running this test")
    sys.exit(0)
if 'GALACTICUS_DATA_PATH' not in os.environ:
    print("SKIPPED: GALACTICUS_DATA_PATH is not set - the Leauthaud2012 target dataset is unavailable")
    sys.exit(0)

# Run the model.
os.makedirs(outputPath, exist_ok=True)
targetsWrite()
with open(logFile, 'w') as log:
    status = subprocess.run(['./'+executable, parameterFile], stdout=log, stderr=subprocess.STDOUT).returncode
if status != 0:
    print("FAILED: model did not run - see "+logFile)
    sys.exit(0)

# Compare each reported log-likelihood with an independent calculation.
analyses   = analysesRead()
discrepant = []
compared   = 0
with h5py.File(outputPath+'/galacticus.hdf5', 'r') as file:
    for groupName, group in sorted(file['analyses'].items()):
        if 'logLikelihood' not in group.attrs:
            continue
        className, options = analysisMatch(groupName, analyses)
        if className is None:
            continue
        reported = float(group.attrs['logLikelihood'])
        expected = logLikelihoodCompute(group, className, options)
        compared += 1
        if not np.isclose(reported, expected, rtol=toleranceRelative, atol=0.0):
            discrepant.append(f'{groupName}: reported {reported:.15e}, expected {expected:.15e}')
if compared < len(analyses):
    print(f"FAILED: only {compared} of {len(analyses)} analyses reported a log-likelihood")
elif discrepant:
    print("FAILED: log-likelihoods differ from an independent calculation:")
    for message in discrepant:
        print("   "+message)
else:
    print(f"SUCCESS: all {compared} relation analysis log-likelihoods match an independent calculation")
sys.exit(0)
