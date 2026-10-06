#!/usr/bin/env python3
"""
Test the `morphologicalFractionGAMAMoffett2016` output analysis.

The analysis reads early-type and total galaxy counts in bins of stellar mass from the GAMA dataset of Moffett et
al. (2016), and from them constructs its target early-type fraction, the confidence interval on that fraction
(Clopper-Pearson in bins of fewer than 100 galaxies, the Wilson score interval with continuity correction
otherwise), and a binomial log-likelihood. Here a small model is run with this analysis, and the target dataset,
its error bars, and the log-likelihood written to the output file are each compared with an independent
calculation from the counts in the dataset.

This guards against regressions of issue #1565, in which the total counts were not read and constructing the
analysis caused a segmentation fault, and of an error in the Wilson score interval, which used the quantile
probability, 1-α/2, in place of the corresponding standard normal deviate.
"""
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

import h5py
import numpy as np
from scipy.stats import beta, norm

# Paths below are relative to the root of the Galacticus source tree, so work from there: this test is run
# both from that root and, by `test-all.py` and the CI workflows, from the `testSuite` directory.
os.chdir(os.environ.get('GALACTICUS_EXEC_PATH', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')))

executable        = 'Galacticus.exe'
outputPath        = 'testSuite/outputs/morphologicalFractionGAMAMoffett2016'
parameterFile     = outputPath+'/parameters.xml'
logFile           = outputPath+'/galacticus.log'
outputFile        = outputPath+'/galacticus.hdf5'
analysisName      = 'morphologicalFractionGAMAMoffett2016'
# Confidence level of the target intervals, and the count below which the Clopper-Pearson interval is used (both as
# set in `source/output/analyses/morphological_fraction/GAMA_Moffett2016.F90`).
confidenceLevel   = 0.683
countExact        = 100.0
# Judged impossible by Galacticus (see `logImpossible` in `source/models/likelihoods/constants.F90`).
logImpossible     = -1.0e-6*np.finfo(np.float64).max
toleranceRelative = 1.0e-6

def parametersWrite():
    """Write a parameter file which adds the analysis to the quick test model."""
    tree       = ET.parse('parameters/quickTest.xml')
    parameters = tree.getroot()
    # Populate every bin of the analysis - an empty bin makes the log-likelihood impossible, and so leaves most of
    # its calculation untested. The quick test model forms no galaxies above about 10¹⁰M☉, so sample more trees
    # over a wider range of mass, and then stretch and shift the model stellar masses (as log₁₀M →
    # log₁₀M+1.5+0.5(log₁₀M-11.3)) using the analysis' own systematic error polynomial to span all bins.
    masses     = parameters.find('mergerTreeBuildMasses')
    masses.find('massTreeMinimum').set('value', '1.0e10')
    masses.find('massTreeMaximum').set('value', '1.0e15')
    masses.find('treesPerDecade' ).set('value', '10'    )
    outputter  = parameters.find('mergerTreeOutputter')
    parameters.remove(outputter)
    multi      = ET.SubElement(parameters, 'mergerTreeOutputter', value='multi'   )
    multi.append(outputter)
    ET.SubElement(multi, 'mergerTreeOutputter', value='analyzer')
    outputTimes = ET.SubElement(parameters, 'outputTimes', value='list')
    ET.SubElement(outputTimes, 'redshifts', value='0.00 0.02 0.04 0.06 0.08')
    analysis    = ET.SubElement(parameters, 'outputAnalysis', value=analysisName)
    ET.SubElement(analysis, 'systematicErrorPolynomialCoefficient', value='1.5 0.5')
    ET.SubElement(parameters, 'outputFileName', value=outputFile)
    tree.write(parameterFile)

def intervalCompute(countEarly, countAll):
    """Return the lower and upper limits of the confidence interval on the early-type fraction in each bin."""
    alpha    = 1.0-confidenceLevel
    fraction = countEarly/countAll
    lower    = np.empty(countAll.size)
    upper    = np.empty(countAll.size)
    for i, (k, n, p) in enumerate(zip(countEarly, countAll, fraction)):
        if n < countExact:
            # Clopper-Pearson interval.
            lower[i] = 0.0 if k <= 0.0 else beta.ppf(      0.5*alpha, k      , n-k+1.0)
            upper[i] = 1.0 if k >= n   else beta.ppf(1.0-0.5*alpha, k+1.0, n-k    )
        else:
            # Wilson score interval with continuity correction.
            z             = norm.ppf(1.0-0.5*alpha)
            argumentLower = z**2-1.0/n+4.0*n*p*(1.0-p)+(4.0*p-2.0)
            argumentUpper = z**2-1.0/n+4.0*n*p*(1.0-p)-(4.0*p-2.0)
            lower[i] = 0.0 if argumentLower < 0.0 else max(0.0, (2.0*n*p+z**2-(z*np.sqrt(argumentLower)+1.0))/2.0/(n+z**2))
            upper[i] = 1.0 if argumentUpper < 0.0 else min(1.0, (2.0*n*p+z**2+(z*np.sqrt(argumentUpper)+1.0))/2.0/(n+z**2))
    return lower, upper

def logLikelihoodCompute(countEarly, countAll, fractionModel):
    """Return the binomial log-likelihood of the target counts given the model early-type fraction."""
    logLikelihood = 0.0
    for k, n, f in zip(countEarly, countAll, fractionModel):
        for count, probability in ((k, f), (n-k, 1.0-f)):
            if count > 0.0:
                if probability <= 0.0:
                    return logImpossible
                logLikelihood += count*np.log(probability)
    return logLikelihood

# Check that prerequisites are met.
if not os.path.exists(executable):
    print("SKIPPED: "+executable+" does not exist - build it before running this test")
    sys.exit(0)
if 'GALACTICUS_DATA_PATH' not in os.environ:
    print("SKIPPED: GALACTICUS_DATA_PATH is not set - the GAMA target dataset is unavailable")
    sys.exit(0)

# Run the model.
os.makedirs(outputPath, exist_ok=True)
parametersWrite()
with open(logFile, 'w') as log:
    status = subprocess.run(['./'+executable, parameterFile], stdout=log, stderr=subprocess.STDOUT).returncode
if status != 0:
    print(f"FAILED: model did not run (exit status {status}) - see "+logFile)
    sys.exit(0)

# Read the counts from the target dataset.
with h5py.File(os.path.join(os.environ['GALACTICUS_DATA_PATH'], 'static', 'observations', 'morphology', 'earlyTypeFractionGAMA.hdf5'), 'r') as file:
    countEarly = file['countEarly'][:]
    countAll   = file['countAll'  ][:]

# Compare the target dataset, its covariance, and the log-likelihood with independent calculations.
failures = []
with h5py.File(outputFile, 'r') as file:
    if 'analyses/'+analysisName not in file:
        print(f"FAILED: no output group 'analyses/{analysisName}' - see "+logFile)
        sys.exit(0)
    group            = file['analyses/'+analysisName]
    fractionModel    = group[group.attrs['yDataset'         ].decode()][:]
    fractionTarget   = group[group.attrs['yDatasetTarget'   ].decode()][:]
    errorLower       = group[group.attrs['yErrorLowerTarget'].decode()][:]
    errorUpper       = group[group.attrs['yErrorUpperTarget'].decode()][:]
    logLikelihood    = float(group.attrs['logLikelihood'])
if np.any(fractionModel <= 0.0) or np.any(fractionModel >= 1.0):
    print(f"FAILED: the model early-type fraction is not in (0,1) in every bin, so the log-likelihood is not fully tested: {fractionModel}")
    sys.exit(0)
fraction     = countEarly/countAll
lower, upper = intervalCompute(countEarly, countAll)
expected     = {
    'target fraction'   : (fractionTarget           , fraction                                                            ),
    'target lower error': (errorLower               , fraction-lower                                                      ),
    'target upper error': (errorUpper               , upper-fraction                                                      ),
    'log-likelihood'    : (np.array([logLikelihood]), np.array([logLikelihoodCompute(countEarly, countAll, fractionModel)]))
}
for name, (reported, independent) in expected.items():
    if reported.shape != independent.shape:
        failures.append(f"{name}: reported shape {reported.shape}, expected {independent.shape}")
    elif not np.allclose(reported, independent, rtol=toleranceRelative, atol=0.0):
        failures.append(f"{name}: reported {reported}, expected {independent}")
if failures:
    print("FAILED: the GAMA morphological fraction analysis differs from an independent calculation:")
    for message in failures:
        print("   "+message)
else:
    print(f"SUCCESS: the GAMA morphological fraction target ({countAll.size} bins), its error bars, and the log-likelihood match an independent calculation")
sys.exit(0)
