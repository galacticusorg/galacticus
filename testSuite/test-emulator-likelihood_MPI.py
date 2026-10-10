#!/usr/bin/env python3
"""Test the emulated likelihood (`posteriorSampleLikelihoodEmulated`) against an independent computation in Python.

Builds an emulator file (with numpy only, as in `test-emulator-predict.py`) holding emulators of two observables over a
design in three parameters with uniform, log-uniform, and normal priors. The emulators take their inputs in an order
different from that of the design and of the active parameters, and one uses only a subset of the parameters, so that any
confusion in mapping parameters to inputs shows. The parameter names differ in length, so that they are padded when stored
as fixed-length strings in the emulator file. The target data include bins which must be excluded: a missing datum, a
datum which is minus infinity (as for an empty bin after a logarithmic transformation), and a datum with zero variance
(which is included only if the emulator variance is). A `grid` posterior simulation then evaluates the likelihood at a set
of points, and the log-likelihoods it reports are compared with those computed in Python from
`Galacticus.Emulation.emulatorFile.predict`, for each likelihood form and with and without the emulator variance. Finally,
a run with a prior differing from that under which the emulators were trained must be rejected.

Posterior sampling requires an MPI build of Galacticus, so the simulations are run under `mpirun`, with the grid points
divided between the processes.

Andrew Benson, Claude (6-October-2026).
"""

import argparse
import glob
import math
import os
import subprocess
import sys
from statistics import NormalDist

import numpy as np

root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
os.environ['GALACTICUS_EXEC_PATH'] = root
sys.path.insert(0, os.path.join(root, 'python'))
os.chdir(root)

from Galacticus.Emulation import emulatorFile as ef  # noqa: E402

directory = 'testSuite/outputs/emulatorLikelihood'
emulatorFileName = f'{directory}/emulator.hdf5'
countPoints = 30
processes = 1
allowRunAsRoot = []

# The parameters, in the order of the design (and of the active parameters), with their priors (as Galacticus parameters, and
# as cumulative distribution and inverse functions).
normal = NormalDist(mu=0.5, sigma=math.sqrt(0.04))
parameters = [
    {'name': 'p/a', 'prior': 'uniform', 'options': {'limitLower': 0.0, 'limitUpper': 2.0}, 'mapper': 'identity',
     'cumulative': lambda x: x / 2.0, 'invert': lambda f: 2.0 * f},
    {'name': 'p/b', 'prior': 'logUniform', 'options': {'limitLower': 1.0e-2, 'limitUpper': 1.0e+1}, 'mapper': 'logarithm',
     'cumulative': lambda x: math.log(x / 1.0e-2) / math.log(1.0e3), 'invert': lambda f: 1.0e-2 * 1.0e3 ** f},
    {'name': 'p/cLonger', 'prior': 'normal', 'options': {'mean': 0.5, 'variance': 0.04}, 'mapper': 'identity',
     'cumulative': normal.cdf, 'invert': normal.inv_cdf},
]
names = [parameter['name'] for parameter in parameters]

# The observables: the names of their inputs (in the emulator's order), and their target data.
countBins = {'observableDiagonal': 6, 'observableCovariance': 5}
inputNames = {'observableDiagonal': ['p/cLonger', 'p/a', 'p/b'], 'observableCovariance': ['p/b', 'p/a']}


def targets():
    """Construct the target data of each observable."""
    yDiagonal = np.array([-0.3, 0.2, np.nan, 0.5, -np.inf, 0.1])
    varianceDiagonal = np.array([0.02, 0.03, 0.02, 0.0, 0.01, 0.04])
    yCovariance = np.array([0.4, -0.2, 0.0, np.nan, 0.3])
    sigma = np.array([0.15, 0.2, 0.1, 0.2, 0.25])
    correlation = 0.5 ** np.abs(np.subtract.outer(np.arange(5), np.arange(5)))
    return {
        'observableDiagonal': (yDiagonal, np.diag(varianceDiagonal)),
        'observableCovariance': (yCovariance, correlation * np.outer(sigma, sigma)),
    }


def buildEmulator(names_, inputs, y, seed):
    """Build an emulator directly: standardize the bins, keep every principal component, and construct a Gaussian process for
    each standardized coefficient with fixed, unequal, hyperparameters."""
    rng = np.random.default_rng(seed)
    binMean = y.mean(axis=0)
    binScale = y.std(axis=0)
    standardized = (y - binMean) / binScale
    _, _, components = np.linalg.svd(standardized, full_matrices=False)
    coefficients = standardized @ components.T
    coefficientMean = coefficients.mean(axis=0)
    coefficientScale = coefficients.std(axis=0)
    gpComponents = []
    for k in range(components.shape[0]):
        target = (coefficients[:, k] - coefficientMean[k]) / coefficientScale[k]
        logAmplitude = np.log(0.5 + 0.3 * k)
        logLengthScales = np.log(np.linspace(0.3, 1.2, inputs.shape[1])) + 0.1 * k
        noiseVariance = 1.0e-3 * (1.0 + rng.uniform(size=inputs.shape[0]))
        covariance = ef.kernel_matern52(inputs, inputs, logAmplitude, logLengthScales)
        covariance[np.diag_indices_from(covariance)] += noiseVariance + 1.0e-10
        gpComponents.append(ef.Component(
            log_amplitude=logAmplitude, log_length_scales=logLengthScales, targets=target, noise_variance=noiseVariance,
            alpha=np.linalg.solve(covariance, target), cholesky_factor=np.linalg.cholesky(covariance),
        ))
    return ef.Emulator(
        kernel='matern52ARD', jitter=1.0e-10, input_names=names_, inputs=inputs, bin_mean=binMean, bin_scale=binScale,
        pca_components=components, coefficient_mean=coefficientMean, coefficient_scale=coefficientScale,
        components=gpComponents, pca_variance_retained=1.0,
    )


def buildFile():
    rng = np.random.default_rng(11)
    quantiles = rng.uniform(size=(countPoints, len(parameters)))
    values = np.array([[parameter['invert'](f) for parameter, f in zip(parameters, row)] for row in quantiles])
    design = ef.Design(
        names=names, priors=[ef.Prior(parameter['prior'], parameter['options']) for parameter in parameters],
        mappers=[parameter['mapper'] for parameter in parameters], quantiles=quantiles, values=values,
    )
    trainingSets, emulators = {}, {}
    for seed, (label, (yTarget, covarianceTarget)) in enumerate(targets().items()):
        columns = [names.index(name) for name in inputNames[label]]
        inputs = quantiles[:, columns]
        u = quantiles
        y = np.array([[np.sin(2.0 * u_[0] + 0.4 * b) + 0.5 * u_[1] * b / countBins[label] - 0.3 * u_[2] ** 2 * (seed + 1)
                       for b in range(countBins[label])] for u_ in u])
        trainingSets[label] = ef.TrainingSet(
            x=np.arange(countBins[label], dtype=float), y=y, root_variance=np.full(y.shape, 0.03),
            mask=np.zeros(y.shape, dtype=np.int8), point_index=np.arange(countPoints), y_target=yTarget,
            covariance_target=covarianceTarget,
        )
        emulators[label] = buildEmulator(inputNames[label], inputs, y, seed)
    content = ef.EmulatorFile(design=design, training_sets=trainingSets, emulators=emulators)
    ef.write(emulatorFileName, content, overwrite=True)
    return content


def logLikelihoodExpected(content, values, forms, includeEmulatorVariance):
    """Compute the log-likelihood at the given parameter values."""
    quantilesAll = {parameter['name']: parameter['cumulative'](value) for parameter, value in zip(parameters, values)}
    logLikelihood = 0.0
    for label, form in zip(countBins, forms):
        emulator = content.emulators[label]
        trainingSet = content.training_sets[label]
        mean, variance = ef.predict(emulator, np.array([quantilesAll[name] for name in emulator.input_names]))
        mean, variance = mean.ravel(), variance.ravel()
        y, covariance = trainingSet.y_target, trainingSet.covariance_target
        varianceTotal = np.diag(covariance) + (variance if includeEmulatorVariance else 0.0)
        with np.errstate(invalid='ignore'):
            included = np.isfinite(y) & np.isfinite(np.diag(covariance)) & (varianceTotal > 0.0)
        residual = (y - mean)[included]
        if form == 'gaussianDiagonal':
            logLikelihood += -0.5 * np.sum(residual ** 2 / varianceTotal[included] + np.log(2.0 * np.pi * varianceTotal[included]))
        else:
            covarianceTotal = covariance[np.ix_(included, included)].copy()
            if includeEmulatorVariance:
                covarianceTotal[np.diag_indices_from(covarianceTotal)] += variance[included]
            _, logDeterminant = np.linalg.slogdet(covarianceTotal)
            logLikelihood += -0.5 * (residual @ np.linalg.solve(covarianceTotal, residual) + logDeterminant
                                     + residual.size * np.log(2.0 * np.pi))
    return logLikelihood


def parameterFile(name, forms, includeEmulatorVariance, limitUpperA=2.0):
    """Write a parameter file for a grid simulation of the emulated likelihood."""
    modelParameters = ''
    for parameter in parameters:
        options = dict(parameter['options'])
        if parameter['name'] == 'p/a':
            options['limitUpper'] = limitUpperA
        prior = '\n'.join(f'        <{key} value="{value}"/>' for key, value in options.items())
        modelParameters += f'''    <modelParameter value="active">
      <name value="{parameter['name']}"/>
      <distributionFunction1DPrior value="{parameter['prior']}">
{prior}
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="{parameter['mapper']}"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0"/>
        <scale  value="1.0e-3"/>
      </distributionFunction1DPerturber>
    </modelParameter>
'''
    emulators = ''.join(f'''    <emulator value="gaussianProcess">
      <fileName value="{emulatorFileName}"/>
      <label    value="{label}"/>
    </emulator>
''' for label in countBins)
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
    <includeEmulatorVariance value="{'true' if includeEmulatorVariance else 'false'}"/>
  </posteriorSampleLikelihood>
  <posteriorSampleSimulation value="grid">
    <logFileRoot       value="{directory}/{name}"/>
    <logFlushCount     value="1"/>
    <appendLogs        value="false"/>
    <outputLikelihoods value="false"/>
{modelParameters}    <posteriorSamples value="sobol">
      <countSamples value="16"/>
      <randomShift  value="true"/>
    </posteriorSamples>
  </posteriorSampleSimulation>
  <randomNumberGenerator value="GSL">
    <seed value="7"/>
  </randomNumberGenerator>
</parameters>
''')
    return path


def run(name, forms, includeEmulatorVariance, limitUpperA=2.0):
    """Run a grid simulation, returning its exit status, and the path of its log."""
    path = parameterFile(name, forms, includeEmulatorVariance, limitUpperA)
    logPath = f'{directory}/{name}.out'
    for stale in glob.glob(f'{directory}/{name}_*.log'):
        os.remove(stale)
    with open(logPath, 'w') as log:
        status = subprocess.run(['mpirun', '--oversubscribe', '-np', str(processes)] + allowRunAsRoot + ['./Galacticus.exe', path],
                                stdout=log, stderr=subprocess.STDOUT, env=dict(os.environ, OMP_NUM_THREADS='1')).returncode
    return status, logPath


def check(content, name, forms, includeEmulatorVariance):
    status, logPath = run(name, forms, includeEmulatorVariance)
    if status != 0:
        print(f'FAILED: grid simulation ({name}) exited with status {status} (see {logPath})')
        return
    rows = []
    for path in sorted(glob.glob(f'{directory}/{name}_*.log')):
        with open(path) as file:
            for line in file:
                columns = line.split()
                rows.append((float(columns[5]), [float(value) for value in columns[6:]]))
    if len(rows) != 16:
        print(f'FAILED: grid simulation ({name}) evaluated {len(rows)} states, expected 16')
        return
    formsAll = forms if len(forms) == len(countBins) else forms * len(countBins)
    errors = [abs(logLikelihood - logLikelihoodExpected(content, values, formsAll, includeEmulatorVariance))
              / max(abs(logLikelihood), 1.0) for logLikelihood, values in rows]
    error = max(errors)
    message = f'log-likelihoods ({name}; maximum relative difference {error:.1e})'
    print(f'SUCCESS: {message}' if error < 1.0e-9 else f'FAILED: {message}')


def main():
    global processes, allowRunAsRoot
    parser = argparse.ArgumentParser()
    parser.add_argument('--processesPerNode', type=int, default=1)
    parser.add_argument('--allowRunAsRoot', type=str, default='no')
    args, _ = parser.parse_known_args()
    # The number of grid points (16) must be divisible by the number of processes.
    processes = 4 if args.processesPerNode >= 4 else (2 if args.processesPerNode >= 2 else 1)
    allowRunAsRoot = ['--allow-run-as-root'] if args.allowRunAsRoot == 'yes' else []
    os.makedirs(directory, exist_ok=True)
    content = buildFile()
    check(content, 'diagonalCovariance', ['gaussianDiagonal', 'gaussianCovariance'], True)
    check(content, 'diagonalCovarianceNoVariance', ['gaussianDiagonal', 'gaussianCovariance'], False)
    check(content, 'covarianceAll', ['gaussianCovariance'], True)
    # A prior differing from that under which the emulators were trained must be rejected.
    status, logPath = run('priorMismatch', ['gaussianDiagonal'], True, limitUpperA=2.5)
    with open(logPath) as file:
        rejected = "the prior of parameter 'p/a' differs" in file.read()
    print('SUCCESS: a mismatched prior is rejected' if status != 0 and rejected
          else f'FAILED: a mismatched prior was not rejected (see {logPath})')


if __name__ == '__main__':
    try:
        main()
    except Exception as error:  # Report any unexpected error as a failure, rather than a bare traceback.
        print(f'FAILED: unexpected error: {type(error).__name__}: {error}')
    sys.exit(0)
