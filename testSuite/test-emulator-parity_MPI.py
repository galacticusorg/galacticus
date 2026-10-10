#!/usr/bin/env python3
"""Test that the emulated likelihood (`posteriorSampleLikelihoodEmulated`) reproduces the likelihoods of the analyses which it
emulates, given the exact model.

Runs a small galaxy formation model (the quick test model) at the points of a design, with two analyses: a
`volumeFunction1D` (a stellar mass function, with a full covariance matrix for its target data) and a `meanFunction1D` (the
mean halo mass as a function of stellar mass, with a diagonal covariance matrix for its target data). Each analysis writes
its own log-likelihood to the model output. The training sets are then collected from those outputs
(`Galacticus.Emulation.train.collect`), exactly as for training an emulator.

For each design point, an emulator is constructed which predicts exactly the model vector of that point, with a variance
in each bin equal to the variance of the model in that bin (the diagonal of the model covariance written by the analysis).
The emulated likelihood evaluated with that emulator (including the emulator variance) must then equal the log-likelihood of
the analysis, which uses the sum of the model and target covariances. This holds for the `gaussianCovariance` form, and
for the `gaussianDiagonal` form when the target covariance is diagonal, provided that the model covariance is diagonal,
which is checked. The test therefore checks that the target data, their covariance, and the model vector reach the emulated
likelihood intact through the analysis outputs, the collector, and the emulator file, and that the emulated likelihood has
the same form and normalization as the likelihoods of the analyses.

Posterior sampling requires an MPI build of Galacticus.

Andrew Benson, Claude (10-October-2026).
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys

import h5py
import lxml.etree as ET
import numpy as np

root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
os.environ['GALACTICUS_EXEC_PATH'] = root
sys.path.insert(0, os.path.join(root, 'python'))
os.chdir(root)

directory = 'testSuite/outputs/emulatorParity'
# The varied parameter, and its prior.
parameterName = 'hotHaloOutflowReincorporation/multiplier'
parameterPrior = (1.0, 10.0)
countDesign = 4
# The analyses: the label, the bins (in log10 stellar mass), and the target data.
massBins = np.arange(6.0, 11.01, 0.5)
labelVolume = 'parityVolume'
labelMean = 'parityMean'


def targets():
    """Construct the target data: for the mass function, with correlated errors; for the mean halo mass, independent errors."""
    yVolume = 10.0 ** (-1.5 - 0.3 * (massBins - 8.0))
    sigmaVolume = 0.3 * yVolume
    correlation = 0.4 ** np.abs(np.subtract.outer(np.arange(massBins.size), np.arange(massBins.size)))
    yMean = 11.0 + 0.5 * (massBins - 8.0)
    return {
        labelVolume: (yVolume, correlation * np.outer(sigmaVolume, sigmaVolume)),
        labelMean: (yMean, np.diag(np.full(massBins.size, 0.2 ** 2))),
    }


def values(array):
    return ' '.join(f'{value:.17g}' for value in np.ravel(array))


def addElements(parent, elements):
    """Add a list of (name, value) elements to an XML element."""
    for name, value in elements:
        ET.SubElement(parent, name, value=value)


def writeBase():
    """Write the base parameter file: the quick test model, with its outputs replaced by the two analyses."""
    tree = ET.parse('parameters/quickTest.xml')
    for outputter in tree.getroot().findall('mergerTreeOutputter'):
        tree.getroot().remove(outputter)
    ET.SubElement(tree.getroot(), 'mergerTreeOutputter', value='analyzer')
    # Use more trees, so that most bins of the analyses are populated.
    tree.getroot().find('mergerTreeBuildMasses').find('treesPerDecade').set('value', '20')
    target = targets()
    analyses = ET.SubElement(tree.getroot(), 'outputAnalysis', value='multi')
    # Settings common to both analyses: stellar mass bins at z=0, with no operators other than the logarithm of stellar mass.
    common = [('xAxisLabel', 'M_star'), ('xAxisIsLog', 'true'), ('propertyLabel', 'massStellar'),
              ('propertyComment', 'Stellar mass'), ('propertyUnits', 'Msun'), ('propertyQuantity', 'solMass'),
              ('propertyIsComoving', 'false'), ('propertyUnitsInSI', '1.98847e30'), ('binCenter', values(massBins)),
              ('bufferCount', '0'), ('outputWeight', values(np.ones(massBins.size))), ('covarianceModel', 'poisson'),
              ('targetLabel', 'synthetic')]

    def operators(analysis):
        ET.SubElement(analysis, 'nodePropertyExtractor', value='massStellar')
        ET.SubElement(analysis, 'outputAnalysisPropertyOperator', value='log10')
        ET.SubElement(analysis, 'outputAnalysisWeightOperator', value='identity')
        ET.SubElement(analysis, 'outputAnalysisDistributionOperator', value='identity')
        ET.SubElement(analysis, 'galacticFilter', value='always')
        outputTimes = ET.SubElement(analysis, 'outputTimes', value='list')
        ET.SubElement(outputTimes, 'redshifts', value='0.0')

    volume = ET.SubElement(analyses, 'outputAnalysis', value='volumeFunction1D')
    addElements(volume, [('label', labelVolume), ('comment', 'Stellar mass function'), ('yAxisLabel', 'phi'),
                         ('yAxisIsLog', 'true')] + common
                + [('distributionLabel', 'massFunction'), ('distributionComment', 'Mass function'),
                   ('distributionUnits', 'Mpc^-3'), ('distributionQuantity', 'Mpc**-3'), ('distributionIsComoving', 'true'),
                   ('distributionUnitsInSI', '3.4036771916e-68'), ('functionValueTarget', values(target[labelVolume][0])),
                   ('functionCovarianceTarget', values(target[labelVolume][1]))])
    operators(volume)
    ET.SubElement(volume, 'outputAnalysisDistributionNormalizer', value='identity')
    unoperator = ET.SubElement(volume, 'unoperatorParameters')
    ET.SubElement(unoperator, 'outputAnalysisPropertyOperator', value='antiLog10')

    mean = ET.SubElement(analyses, 'outputAnalysis', value='meanFunction1D')
    addElements(mean, [('label', labelMean), ('comment', 'Mean halo mass'), ('yAxisLabel', 'M_halo'),
                       ('yAxisIsLog', 'false')] + common
                + [('meanLabel', 'massHalo'), ('meanComment', 'Mean log10 halo mass'), ('meanUnits', 'dex'),
                   ('meanQuantity', 'dimensionless'), ('meanIsComoving', 'false'), ('meanUnitsInSI', '1.0'),
                   ('meanValueTarget', values(target[labelMean][0])), ('meanCovarianceTarget', values(target[labelMean][1]))])
    operators(mean)
    weight = ET.SubElement(mean, 'weight')
    ET.SubElement(weight, 'nodePropertyExtractor', value='massBasic')
    ET.SubElement(weight, 'outputAnalysisPropertyOperator', value='log10')
    unoperator = ET.SubElement(mean, 'unoperator')
    ET.SubElement(unoperator, 'outputAnalysisPropertyOperator', value='antiLog10')
    path = f'{directory}/base.xml'
    tree.write(path, xml_declaration=True, encoding='UTF-8')
    return path


def modelParameter():
    """The active model parameter, as XML."""
    return f'''    <modelParameter value="active">
      <name value="{parameterName}"/>
      <distributionFunction1DPrior value="uniform">
        <limitLower value="{parameterPrior[0]}"/>
        <limitUpper value="{parameterPrior[1]}"/>
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0"/>
        <scale  value="1.0e-3"/>
      </distributionFunction1DPerturber>
    </modelParameter>
'''


def writeDesign():
    path = f'{directory}/design.xml'
    with open(path, 'w') as file:
        file.write(f'''<parameters>
  <formatVersion>2</formatVersion>
  <task value="emulatorDesign">
    <designFileName     value="{directory}/design.hdf5"/>
    <changeFilesRoot    value="{directory}/changes/point"/>
    <outputFileNameRoot value="{directory}/models/point"/>
{modelParameter()}    <posteriorSamples value="sobol">
      <countSamples value="{countDesign}"/>
      <randomShift  value="true"/>
    </posteriorSamples>
  </task>
  <randomNumberGenerator value="GSL">
    <seed value="5"/>
  </randomNumberGenerator>
</parameters>
''')
    return path


def exactEmulator(ef, design, trainingSet, row):
    """Construct an emulator which predicts exactly the model vector of one row of a training set, with a variance in each bin
    equal to the variance of the model. The emulator has one principal component per bin (the identity), and a Gaussian process
    for each with zero targets, unit amplitude, and length scales so short that its covariance with any evaluation point away
    from the training points underflows to zero: its mean is then zero and its variance is the amplitude, so the prediction is
    the bin mean (the model vector) with variance the square of the bin scale (the model variance)."""
    countBins = trainingSet.y.shape[1]
    inputs = design.quantiles[trainingSet.point_index, :]
    logLengthScales = np.full(inputs.shape[1], np.log(1.0e-6))
    components = []
    for _ in range(countBins):
        noiseVariance = np.full(inputs.shape[0], 1.0e-3)
        covariance = ef.kernel_matern52(inputs, inputs, 0.0, logLengthScales)
        covariance[np.diag_indices_from(covariance)] += noiseVariance + 1.0e-10
        components.append(ef.Component(
            log_amplitude=0.0, log_length_scales=logLengthScales, targets=np.zeros(inputs.shape[0]),
            noise_variance=noiseVariance, alpha=np.zeros(inputs.shape[0]), cholesky_factor=np.linalg.cholesky(covariance),
        ))
    return ef.Emulator(
        kernel='matern52ARD', jitter=1.0e-10, input_names=list(design.names), inputs=inputs,
        bin_mean=trainingSet.y[row].copy(), bin_scale=trainingSet.root_variance[row].copy(),
        pca_components=np.eye(countBins), coefficient_mean=np.zeros(countBins), coefficient_scale=np.ones(countBins),
        components=components, pca_variance_retained=1.0,
    )


def writeSampler(name, emulatorFileName, labels, forms):
    """Write a parameter file for a grid simulation of the emulated likelihood."""
    emulators = ''.join(f'''    <emulator value="gaussianProcess">
      <fileName value="{emulatorFileName}"/>
      <label    value="{label}"/>
    </emulator>
''' for label in labels)
    path = f'{directory}/{name}.xml'
    with open(path, 'w') as file:
        file.write(f'''<parameters>
  <formatVersion>2</formatVersion>
  <outputFileName value="{directory}/{name}.hdf5"/>
  <task value="posteriorSample">
    <initializeNodeClassHierarchy value="false"/>
  </task>
  <posteriorSampleLikelihood value="emulated">
{emulators}    <likelihoodForms         value="{' '.join(forms)}"/>
    <includeEmulatorVariance value="true"/>
  </posteriorSampleLikelihood>
  <posteriorSampleSimulation value="grid">
    <logFileRoot       value="{directory}/{name}"/>
    <logFlushCount     value="1"/>
    <appendLogs        value="false"/>
    <outputLikelihoods value="false"/>
{modelParameter()}    <posteriorSamples value="sobol">
      <countSamples value="2"/>
      <randomShift  value="true"/>
    </posteriorSamples>
  </posteriorSampleSimulation>
  <randomNumberGenerator value="GSL">
    <seed value="13"/>
  </randomNumberGenerator>
</parameters>
''')
    return path


def run(command, logName):
    """Run a command, logging its output, and return its exit status."""
    with open(f'{directory}/{logName}', 'w') as log:
        return subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=dict(os.environ, OMP_NUM_THREADS='1')).returncode


def directLikelihoods(path):
    """Read the log-likelihood of each analysis, and the largest off-diagonal element of its model covariance relative to the
    diagonal, from a model output."""
    results = {}
    with h5py.File(path, 'r') as file:
        for label in (labelVolume, labelMean):
            group = file['analyses'][label]
            covariance = group[group.attrs['yCovariance'].decode() if isinstance(group.attrs['yCovariance'], bytes)
                               else group.attrs['yCovariance']][()]
            diagonal = np.sqrt(np.outer(np.diag(covariance), np.diag(covariance)))
            offDiagonal = covariance - np.diag(np.diag(covariance))
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = np.where(diagonal > 0.0, np.abs(offDiagonal) / diagonal, np.abs(offDiagonal))
            results[label] = (float(group.attrs['logLikelihood']), float(np.max(ratio)))
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--processesPerNode', type=int, default=1)
    parser.add_argument('--allowRunAsRoot', type=str, default='no')
    args, _ = parser.parse_known_args()
    from Galacticus.Emulation import emulatorFile as ef
    from Galacticus.Emulation.campaign import Campaign, output_file_path
    from Galacticus.Emulation.train import Observable, collect
    processes = max(args.processesPerNode, 1)
    allowRunAsRoot = ['--allow-run-as-root'] if args.allowRunAsRoot == 'yes' else []
    if args.allowRunAsRoot == 'yes':
        # Allow the MPI build of Galacticus to run as root without `mpirun` (as for the design and the campaign).
        os.environ['OMPI_ALLOW_RUN_AS_ROOT'] = '1'
        os.environ['OMPI_ALLOW_RUN_AS_ROOT_CONFIRM'] = '1'
    shutil.rmtree(directory, ignore_errors=True)
    os.makedirs(directory)

    # Generate the design, and run the models.
    base = writeBase()
    if run(['./Galacticus.exe', writeDesign()], 'design.log') != 0:
        print(f'FAILED: emulatorDesign task failed (see {directory}/design.log)')
        return
    campaign = Campaign.create(f'{directory}/campaign.json', f'{directory}/design.hdf5', base)
    summary = campaign.run_local(parallel=processes, threads=1)
    if summary.get('complete', 0) != countDesign:
        failed = [run_ for run_ in campaign.runs if run_.status != 'complete']
        print(f'FAILED: {len(failed)} of {countDesign} model runs did not complete (e.g. run {failed[0].index}: {failed[0].message}; '
              f'see {failed[0].log_file})')
        return
    # Collect the training sets, exactly as for training an emulator.
    design, trainingSets, _ = collect(f'{directory}/design.hdf5', [Observable(labelVolume), Observable(labelMean)])
    # The checks: the emulator(s) to use, and their likelihood forms.
    checks = [
        ('volumeCovariance', [labelVolume], ['gaussianCovariance']),
        ('meanDiagonal', [labelMean], ['gaussianDiagonal']),
        ('meanCovariance', [labelMean], ['gaussianCovariance']),
        ('both', [labelVolume, labelMean], ['gaussianCovariance', 'gaussianDiagonal']),
    ]
    errors = {name: [] for name, _, _ in checks}
    for run_ in campaign.runs:
        point = run_.point
        direct = directLikelihoods(output_file_path(run_.output_file))
        for label, (_, offDiagonal) in direct.items():
            if offDiagonal > 1.0e-12:
                print(f'FAILED: the model covariance of {label} at point {point} is not diagonal (off-diagonal elements up to '
                      f'{offDiagonal:.1e} of the diagonal), so the parity test is not valid')
                return
        for label, (logLikelihood, _) in direct.items():
            print(f'point {point}: direct log-likelihood of {label} is {logLikelihood:.10g}')
        # Construct the exact emulators of this point.
        emulators = {}
        for label, trainingSet in trainingSets.items():
            row = int(np.nonzero(trainingSet.point_index == point)[0][0])
            emulators[label] = exactEmulator(ef, design, trainingSet, row)
            mean, variance = ef.predict(emulators[label], np.full(len(design.names), 0.5))
            if (not np.array_equal(mean, trainingSet.y[row])
                    or not np.allclose(variance, trainingSet.root_variance[row] ** 2, rtol=1.0e-14, atol=0.0)):
                print(f'FAILED: the exact emulator of {label} at point {point} does not reproduce the model')
                return
        emulatorFileName = f'{directory}/emulator{point:04d}.hdf5'
        ef.write(emulatorFileName, ef.EmulatorFile(design=design, training_sets=trainingSets, emulators=emulators), overwrite=True)
        # Evaluate the emulated likelihood, and compare with the direct likelihood.
        for name, labels, forms in checks:
            nameRun = f'{name}{point:04d}'
            for stale in glob.glob(f'{directory}/{nameRun}_*.log'):
                os.remove(stale)
            status = run(['mpirun', '--oversubscribe', '-np', '1'] + allowRunAsRoot
                         + ['./Galacticus.exe', writeSampler(nameRun, emulatorFileName, labels, forms)], f'{nameRun}.out')
            if status != 0:
                print(f'FAILED: grid simulation ({nameRun}) exited with status {status} (see {directory}/{nameRun}.out)')
                return
            logLikelihoods = []
            for path in sorted(glob.glob(f'{directory}/{nameRun}_*.log')):
                with open(path) as file:
                    logLikelihoods += [float(line.split()[5]) for line in file if not line.startswith('#')]
            if len(logLikelihoods) != 2:
                print(f'FAILED: grid simulation ({nameRun}) evaluated {len(logLikelihoods)} states, expected 2')
                return
            expected = sum(direct[label][0] for label in labels)
            errors[name] += [abs(logLikelihood - expected) / max(abs(expected), 1.0) for logLikelihood in logLikelihoods]
    for name, labels, forms in checks:
        error = max(errors[name])
        message = (f"emulated log-likelihood of {' and '.join(labels)} ({' and '.join(forms)}) equals that of the "
                   f'analysis (maximum relative difference {error:.1e})')
        print(f'SUCCESS: {message}' if error < 1.0e-9 else f'FAILED: {message}')


if __name__ == '__main__':
    try:
        main()
    except Exception as error:  # Report any unexpected error as a failure, rather than a bare traceback.
        print(f'FAILED: unexpected error: {type(error).__name__}: {error}')
    sys.exit(0)
