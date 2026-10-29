#!/usr/bin/env python3
"""Test the counts likelihood of the progenitor mass function analysis, with and without the variance of the model.

With `likelihoodInCounts`, the analysis evaluates a negative binomial likelihood of the target counts in each bin, with an
over-dispersion c² = C² (C being `rootVarianceTargetFractional`), or, if `covarianceTargetOnly` is false, c² = C² + v + C² v,
where v is the fractional variance of the model in the bin (from the diagonal of the model covariance). If c² = 0 it is a
Poisson likelihood. Here a small model is run with four progenitor mass function analyses, which differ only in C (0.1, or 0)
and in whether the model variance is included, and the log-likelihood reported by each is compared with an independent
calculation from the model, its covariance, and the target counts written to the output file.

Andrew Benson, Claude (29-October-2026).
"""

import os
import subprocess
import sys

import h5py
import lxml.etree as ET
import numpy as np
from scipy.special import gammaln

root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
os.environ['GALACTICUS_EXEC_PATH'] = root
os.chdir(root)

directory = 'testSuite/outputs/progenitorMassFunctionLikelihood'
massRatioBins = np.logspace(-2.0, 0.0, 8)
massRatioLikelihood = (0.015, 1.0)
toleranceRelative = 1.0e-9
# The analyses: label, C, and whether the model variance is included.
analyses = [('TargetOnly', 0.1, False), ('ModelVariance', 0.1, True), ('PoissonTargetOnly', 0.0, False),
            ('PoissonModelVariance', 0.0, True)]


def values(array):
    return ' '.join(f'{value:.17g}' for value in np.ravel(array))


def writeParameters():
    """The quick test model, with its outputs replaced by the progenitor mass function analyses."""
    tree = ET.parse('parameters/quickTest.xml')
    for name in ('mergerTreeOutputter', 'outputTimes', 'outputFileName'):
        for element in tree.getroot().findall(name):
            tree.getroot().remove(element)
    ET.SubElement(tree.getroot(), 'mergerTreeOutputter', value='analyzer')
    ET.SubElement(tree.getroot(), 'outputFileName', value=f'{directory}/galacticus.hdf5')
    tree.getroot().find('mergerTreeBuildMasses').find('treesPerDecade').set('value', '20')
    outputTimes = ET.SubElement(tree.getroot(), 'outputTimes', value='list')
    ET.SubElement(outputTimes, 'redshifts', value='0.0 0.5')
    errors = ET.SubElement(tree.getroot(), 'nbodyHaloMassError', value='powerLaw')
    for name, value in (('normalization', '0.135'), ('exponent', '-0.333333'), ('fractionalErrorHighMass', '0.01'),
                        ('correlationModelTrivial', 'false'), ('correlationNormalization', '1.0'),
                        ('correlationMassExponent', '1.0'), ('correlationRedshiftExponent', '0.0')):
        ET.SubElement(errors, name, value=value)
    # Target counts, with some empty bins.
    value = 0.3 * massRatioBins ** -0.9
    count = np.round(40.0 * value * massRatioBins ** 0.6)
    count[3] = 0
    multi = ET.SubElement(tree.getroot(), 'outputAnalysis', value='multi')
    for label, rootVariance, modelVariance in analyses:
        analysis = ET.SubElement(multi, 'outputAnalysis', value='progenitorMassFunction')
        for name, setting in (('label', label), ('comment', 'Progenitor mass function'), ('targetLabel', 'synthetic'),
                              ('massRatioMinimum', f'{massRatioBins[0]}'), ('massRatioMaximum', f'{massRatioBins[-1]}'),
                              ('countMassRatio', f'{massRatioBins.size}'), ('massParentMinimum', '1.0e12'),
                              ('massParentMaximum', '1.0e13'), ('redshiftProgenitor', '0.5'), ('redshiftParent', '0.0'),
                              ('massRatioLikelihoodMinimum', f'{massRatioLikelihood[0]}'),
                              ('massRatioLikelihoodMaximum', f'{massRatioLikelihood[1]}'), ('likelihoodInCounts', 'true'),
                              ('covarianceTargetOnly', 'false' if modelVariance else 'true'), ('covarianceDiagonalize', 'true'),
                              ('rootVarianceTargetFractional', f'{rootVariance}'), ('alwaysIsolatedOnly', 'false'),
                              ('functionValueTarget', values(value)), ('functionCovarianceTarget', values(np.diag((0.2 * value) ** 2))),
                              ('functionCountTarget', ' '.join(f'{int(n)}' for n in count)),
                              ('virialDensityContrastDefinition', 'sphericalCollapseClsnlssMttrCsmlgclCnstnt')):
            ET.SubElement(analysis, name, value=setting)
    path = f'{directory}/parameters.xml'
    tree.write(path, xml_declaration=True, encoding='UTF-8')
    return path


def logLikelihood(group, rootVariance, modelVariance):
    """An independent calculation of the log-likelihood of an analysis."""
    x = group[group.attrs['xDataset'].decode()][:]
    y = group[group.attrs['yDataset'].decode()][:]
    covariance = group[group.attrs['yCovariance'].decode()][:]
    countTarget = group['countTarget'][:]
    countModel = y * group['countPerValue'][:]
    included = (x > group.attrs['xLikelihoodMinimum']) & (x < group.attrs['xLikelihoodMaximum'])
    total = 0.0
    for i in np.flatnonzero(included):
        n, mu = countTarget[i], countModel[i]
        if mu <= 0.0:
            if n > 0:
                return -np.inf
            continue
        dispersion = rootVariance ** 2
        if modelVariance:
            dispersion += covariance[i, i] / y[i] ** 2 * (1.0 + dispersion)
        if dispersion > 0.0:
            k = 1.0 / dispersion
            total += n * np.log(mu) - gammaln(n + 1.0) + gammaln(k + n) - gammaln(k) - n * np.log(k + mu) - k * np.log1p(mu / k)
        else:
            total += n * np.log(mu) - gammaln(n + 1.0) - mu
    return total


def main():
    os.makedirs(directory, exist_ok=True)
    parameters = writeParameters()
    with open(f'{directory}/galacticus.log', 'w') as log:
        status = subprocess.run(['./Galacticus.exe', parameters], stdout=log, stderr=subprocess.STDOUT,
                                env=dict(os.environ, OMP_NUM_THREADS='1')).returncode
    if status != 0:
        print(f'FAILED: the model did not run (see {directory}/galacticus.log)')
        return
    failures = []
    reported = {}
    with h5py.File(f'{directory}/galacticus.hdf5', 'r') as file:
        for label, rootVariance, modelVariance in analyses:
            group = file[f'analyses/progenitorMassFunction{label}']
            reported[label] = float(group.attrs['logLikelihood'])
            expected = logLikelihood(group, rootVariance, modelVariance)
            if not (np.isfinite(reported[label]) and abs(reported[label] - expected) <= toleranceRelative * abs(expected)):
                failures.append(f'{label}: reported {reported[label]!r}, expected {expected!r}')
            print(f'{label}: log-likelihood {reported[label]:.10g} (independent calculation {expected:.10g})')
            if modelVariance:
                # The model variance must be large enough to matter, or this test could not detect its omission.
                y = group[group.attrs['yDataset'].decode()][:]
                covariance = group[group.attrs['yCovariance'].decode()][:]
                fractional = np.diag(covariance)[y > 0.0] / y[y > 0.0] ** 2
                if fractional.max() < 1.0e-3:
                    failures.append(f'{label}: the model variance is too small (at most {fractional.max():.2e}) to test')
    # Including the model variance must change the likelihood.
    for withVariance, without in (('ModelVariance', 'TargetOnly'), ('PoissonModelVariance', 'PoissonTargetOnly')):
        if reported[withVariance] == reported[without]:
            failures.append(f'{withVariance}: the log-likelihood is unchanged by including the model variance')
    if failures:
        print('FAILED: progenitor mass function counts likelihoods:\n  ' + '\n  '.join(failures))
    else:
        print(f'SUCCESS: all {len(analyses)} progenitor mass function counts likelihoods match an independent calculation')


if __name__ == '__main__':
    main()
    sys.exit(0)
