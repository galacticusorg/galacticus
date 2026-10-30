#!/usr/bin/env python3
"""Test the Galacticus Gaussian process emulator against the reference implementation of the emulator file format.

Builds an emulator file in Python (with numpy only: the principal components and Gaussian processes are constructed
directly, with deliberately unequal hyperparameters---a different amplitude for each component, a different length
scale for each input, and noise varying from point to point---so that any confusion of inputs, components, or array
orientation shows), then evaluates it at a set of points with the `emulatorPredict` task, and checks that the predicted
means and variances agree with `Galacticus.Emulation.emulatorFile.predict`. This is done both with the Cholesky factors
of the Gaussian processes stored in the file, and without them (so that Galacticus must compute them).

Andrew Benson, Claude (5-October-2026).
"""

import os
import subprocess
import sys

import h5py
import numpy as np

root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
os.environ['GALACTICUS_EXEC_PATH'] = root
sys.path.insert(0, os.path.join(root, 'python'))
os.chdir(root)

from Galacticus.Emulation import emulatorFile as ef  # noqa: E402

directory = 'testSuite/outputs/emulatorPredict'
countPoints, countInputs, countBins, countEvaluate = 25, 3, 6, 12
label = 'testObservable'


def buildEmulator(inputs, y, storeFactor):
    """Build an emulator directly: standardize the bins, keep every principal component, and construct a Gaussian process for
    each standardized coefficient with fixed, unequal, hyperparameters."""
    rng = np.random.default_rng(5)
    binMean = y.mean(axis=0)
    binScale = y.std(axis=0)
    standardized = (y - binMean) / binScale
    _, _, components = np.linalg.svd(standardized, full_matrices=False)
    coefficients = standardized @ components.T
    coefficientMean = coefficients.mean(axis=0)
    coefficientScale = coefficients.std(axis=0)
    gpComponents = []
    for k in range(components.shape[0]):
        targets = (coefficients[:, k] - coefficientMean[k]) / coefficientScale[k]
        logAmplitude = np.log(0.5 + 0.3 * k)
        logLengthScales = np.log([0.2, 0.5, 1.3]) + 0.1 * k
        noiseVariance = 1.0e-4 * (1.0 + rng.uniform(size=countPoints))
        covariance = ef.kernel_matern52(inputs, inputs, logAmplitude, logLengthScales)
        covariance[np.diag_indices_from(covariance)] += noiseVariance + 1.0e-10
        gpComponents.append(ef.Component(
            log_amplitude=logAmplitude,
            log_length_scales=logLengthScales,
            targets=targets,
            noise_variance=noiseVariance,
            alpha=np.linalg.solve(covariance, targets),
            cholesky_factor=np.linalg.cholesky(covariance) if storeFactor else None,
        ))
    return ef.Emulator(
        kernel='matern52ARD', jitter=1.0e-10, input_names=['p/a', 'p/b', 'p/c'], inputs=inputs,
        bin_mean=binMean, bin_scale=binScale, pca_components=components,
        coefficient_mean=coefficientMean, coefficient_scale=coefficientScale, components=gpComponents,
        pca_variance_retained=1.0,
    )


def buildFile(path, storeFactor):
    rng = np.random.default_rng(1)
    inputs = rng.uniform(size=(countPoints, countInputs))
    x = np.linspace(9.0, 11.5, countBins)
    y = np.array([[np.sin(3.0 * u[0] + 0.3 * b) + u[1] * b / countBins - 0.5 * u[2] ** 2 for b in range(countBins)] for u in inputs])
    design = ef.Design(
        names=['p/a', 'p/b', 'p/c'], priors=[ef.Prior('uniform', {'limitLower': 0.0, 'limitUpper': 1.0})] * countInputs,
        mappers=['identity'] * countInputs, quantiles=inputs, values=inputs,
    )
    trainingSet = ef.TrainingSet(
        x=x, y=y, root_variance=np.full(y.shape, 0.01), mask=np.zeros(y.shape, dtype=np.int8), point_index=np.arange(countPoints),
        y_target=np.linspace(-1.0, 1.0, countBins), covariance_target=np.diag(np.full(countBins, 0.04)),
    )
    content = ef.EmulatorFile(design=design, training_sets={label: trainingSet}, emulators={label: buildEmulator(inputs, y, storeFactor)})
    ef.write(path, content, overwrite=True)
    return content


def check(name, storeFactor):
    emulatorFileName = f'{directory}/emulator_{name}.hdf5'
    content = buildFile(emulatorFileName, storeFactor)
    points = np.random.default_rng(2).uniform(size=(countEvaluate, countInputs))
    quantilesFileName = f'{directory}/quantiles.hdf5'
    with h5py.File(quantilesFileName, 'w') as file:
        file.create_dataset('quantiles', data=points)
    outputFileName = f'{directory}/predictions_{name}.hdf5'
    parameterFileName = f'{directory}/predict_{name}.xml'
    with open(parameterFileName, 'w') as file:
        file.write(f'''<parameters>
  <formatVersion>2</formatVersion>
  <task value="emulatorPredict">
    <quantilesFileName value="{quantilesFileName}"/>
    <outputFileName    value="{outputFileName}"/>
    <emulator value="gaussianProcess">
      <fileName value="{emulatorFileName}"/>
      <label    value="{label}"/>
    </emulator>
  </task>
</parameters>
''')
    with open(f'{directory}/predict_{name}.log', 'w') as log:
        status = subprocess.run(['./Galacticus.exe', parameterFileName], stdout=log, stderr=subprocess.STDOUT).returncode
    if status != 0:
        print(f'FAILED: emulatorPredict task ({name}) exited with status {status} (see {directory}/predict_{name}.log)')
        return
    with h5py.File(outputFileName, 'r') as file:
        mean, variance, x = file['mean'][:], file['variance'][:], file['x'][:]
    meanExpected, varianceExpected = ef.predict(content.emulators[label], points)
    if mean.shape != meanExpected.shape:
        print(f'FAILED: predictions ({name}) have shape {mean.shape}, expected {meanExpected.shape}')
        return
    errorMean = np.max(np.abs(mean - meanExpected) / np.maximum(np.abs(meanExpected), 1.0))
    errorVariance = np.max(np.abs(variance - varianceExpected) / np.maximum(varianceExpected, 1.0e-12))
    print(f'SUCCESS: abscissae of the outputs ({name})' if np.array_equal(x, content.training_sets[label].x)
          else f'FAILED: abscissae of the outputs ({name}) differ')
    print(f'SUCCESS: predicted means agree ({name}; maximum relative difference {errorMean:.1e})' if errorMean < 1.0e-10
          else f'FAILED: predicted means disagree ({name}; maximum relative difference {errorMean:.1e})')
    print(f'SUCCESS: predicted variances agree ({name}; maximum relative difference {errorVariance:.1e})' if errorVariance < 1.0e-7
          else f'FAILED: predicted variances disagree ({name}; maximum relative difference {errorVariance:.1e})')


def main():
    os.makedirs(directory, exist_ok=True)
    check('storedFactor', storeFactor=True)
    check('computedFactor', storeFactor=False)


if __name__ == '__main__':
    try:
        main()
    except Exception as error:  # Report any unexpected error as a failure, rather than a bare traceback.
        print(f'FAILED: unexpected error: {type(error).__name__}: {error}')
    sys.exit(0)
